#!/usr/bin/env python3
"""
fixedwing_radar_sim.py

Realistic slow/maneuvering fixed-wing drone trajectories + scanning-radar detections,
for testing track velocity estimation and a "turning" indicator flag.

Layers (each can be checked on its own):
  1. Coordinated-turn dynamics (bank -> turn rate), roll-rate lag, stall floor that rises with bank,
     speed bleed in tight turns.
  2. Autopilot guidance (L1-style): waypoint legs, loiter orbits, break turns, straight legs.
  3. Wind: steady wind + Gauss-Markov gusts (radar sees GROUND velocity).
  4. Radar: rotating antenna (hit time = when beam crosses target azimuth), range/azimuth noise,
     Pd < 1, optional false alarms.
  5. Monte Carlo runner with randomized airframe, wind, scenario, range.

Outputs (CSV, open in Excel or MATLAB readtable):
  truth.csv       - true state every --truth-every steps, incl. 'turning' label (|bank| > 10 deg)
  detections.csv  - radar plots with truth attached for scoring
  runs.csv        - parameters used for each run
  run_XX.png      - optional plots (--plot)

Conventions: x = East, y = North (m). Heading/course/azimuth clockwise from North.
Positive bank = right turn.

Example:
  python fixedwing_radar_sim.py --runs 50 --rpm 30 --plot
  python fixedwing_radar_sim.py --runs 10 --scenario loiter --seed 7
"""
import argparse
import csv
import math
import os
from dataclasses import dataclass, asdict

import numpy as np

G = 9.81
D2R = math.pi / 180.0
TWO_PI = 2.0 * math.pi
TURN_LABEL_BANK_DEG = 10.0


# ----------------------------------------------------------------------------- parameters
@dataclass
class Airframe:
    v_stall: float = 12.0          # m/s, wings level
    v_cruise: float = 20.0         # m/s
    v_max: float = 30.0            # m/s
    bank_max_deg: float = 45.0     # autopilot bank limit
    roll_tau: float = 0.5          # s, roll response time constant
    roll_rate_max_deg: float = 60.0
    speed_tau: float = 3.0         # s, throttle/speed response
    turn_bleed: float = 0.3        # fractional speed loss per (n-1) of load factor
    stall_margin: float = 1.2      # autopilot keeps V >= margin * stall(bank)


@dataclass
class Wind:
    speed: float = 5.0             # m/s steady
    dir_from_deg: float = 270.0    # meteorological: direction wind comes FROM
    gust_sigma: float = 1.0        # m/s per axis
    gust_tau: float = 3.0          # s correlation time


@dataclass
class Radar:
    x: float = 0.0
    y: float = 0.0
    rpm: float = 30.0              # 30 rpm -> 2 s scan, 20 rpm -> 3 s
    range_sigma: float = 15.0      # m
    az_sigma_deg: float = 0.3      # deg
    pd: float = 0.9                # probability of detection per scan
    false_per_scan: float = 0.0    # mean false alarms per scan near the target
    false_box: float = 500.0       # m, half-width of box false alarms appear in


# ----------------------------------------------------------------------------- helpers
def clip(v, lo, hi):
    return lo if v < lo else hi if v > hi else v


def wrap_pi(a):
    return (a + math.pi) % TWO_PI - math.pi


def bearing(dx, dy):
    return math.atan2(dx, dy)


def l1_bank(course, gs, brg_to_point, dist_to_point, bank_max):
    """L1 guidance: lateral accel = 2 V^2 sin(eta) / L1 -> bank angle."""
    eta = clip(wrap_pi(brg_to_point - course), -math.pi / 2, math.pi / 2)
    a_lat = 2.0 * gs * gs * math.sin(eta) / max(dist_to_point, 1.0)
    return clip(math.atan2(a_lat, G), -bank_max, bank_max)


# ----------------------------------------------------------------------------- mission segments
def seg_wp(pts, v=None, accept=None):
    return {"type": "wp", "pts": [tuple(p) for p in pts], "v": v, "accept": accept}


def seg_straight(dist, v=None):
    return {"type": "straight", "dist": dist, "v": v}


def seg_loiter(radius, turns, direction=1, center=None, v=None):
    """center=None -> orbit starts right where the aircraft is (circle beside current course)."""
    return {"type": "loiter", "radius": radius, "turns": turns, "dir": direction,
            "center": center, "v": v}


def seg_break(direction, angle_deg, v=None):
    return {"type": "break", "dir": direction, "angle": angle_deg * D2R, "v": v}


class Guidance:
    def __init__(self, mission, af, L1=60.0):
        self.mission, self.af, self.L1 = mission, af, L1
        self.i, self.s = 0, {}
        self.bank_max = af.bank_max_deg * D2R

    def done(self):
        return self.i >= len(self.mission)

    def current_type(self):
        return "done" if self.done() else self.mission[self.i]["type"]

    def command(self, x, y, course, gs, heading):
        while not self.done():
            seg = self.mission[self.i]
            out = getattr(self, "_" + seg["type"])(seg, x, y, course, gs, heading)
            if out is not None:
                return out
            self.i += 1
            self.s = {}
        return None

    def _v(self, seg):
        return seg["v"] if seg["v"] else self.af.v_cruise

    def _wp(self, seg, x, y, course, gs, heading):
        s, pts = self.s, seg["pts"]
        s.setdefault("prev", (x, y))
        k = s.setdefault("k", 0)
        while k < len(pts):
            (ax, ay), (bx, by) = s["prev"], pts[k]
            L = max(math.hypot(bx - ax, by - ay), 1e-6)
            ux, uy = (bx - ax) / L, (by - ay) / L
            along = (x - ax) * ux + (y - ay) * uy
            r_turn = gs * gs / (G * math.tan(self.bank_max))
            accept = seg["accept"] or max(30.0, r_turn)
            if math.hypot(bx - x, by - y) < accept or along >= L:
                s["prev"], k = (bx, by), k + 1
                s["k"] = k
                continue
            pa = min(along + self.L1, L)
            px, py = ax + ux * pa, ay + uy * pa
            bank = l1_bank(course, gs, bearing(px - x, py - y),
                           math.hypot(px - x, py - y), self.bank_max)
            return bank, self._v(seg)
        return None

    def _straight(self, seg, x, y, course, gs, heading):
        if "pts" not in self.s:
            self.s["pts_seg"] = seg_wp([(x + seg["dist"] * math.sin(course),
                                         y + seg["dist"] * math.cos(course))], v=seg["v"])
            self.s["pts"] = True
        return self._wp(self.s["pts_seg"], x, y, course, gs, heading)

    def _loiter(self, seg, x, y, course, gs, heading):
        s, R, d = self.s, seg["radius"], seg["dir"]
        if "c" not in s:
            if seg["center"] is not None:
                s["c"] = seg["center"]
            else:  # circle center to the side we will turn toward
                side = course + d * math.pi / 2
                s["c"] = (x + R * math.sin(side), y + R * math.cos(side))
            s["acc"], s["th"] = None, None
        cx, cy = s["c"]
        dist = math.hypot(x - cx, y - cy)
        th = bearing(x - cx, y - cy)
        if s["acc"] is None:
            if abs(dist - R) < 0.25 * R:
                s["acc"], s["th"] = 0.0, th
        else:
            s["acc"] += d * wrap_pi(th - s["th"])
            s["th"] = th
            if s["acc"] >= TWO_PI * seg["turns"]:
                return None
        phi = clip(self.L1 / R, 0.3, 1.5)
        tt = th + d * phi
        px, py = cx + R * math.sin(tt), cy + R * math.cos(tt)
        bank = l1_bank(course, gs, bearing(px - x, py - y),
                       math.hypot(px - x, py - y), self.bank_max)
        return bank, self._v(seg)

    def _break(self, seg, x, y, course, gs, heading):
        s = self.s
        if "acc" not in s:
            s["acc"], s["h"] = 0.0, heading
        s["acc"] += seg["dir"] * wrap_pi(heading - s["h"])
        s["h"] = heading
        if s["acc"] >= seg["angle"]:
            return None
        return seg["dir"] * self.bank_max, self._v(seg)


# ----------------------------------------------------------------------------- scenarios
def to_world(pts_fr, x0, y0, h0):
    """(forward, right) offsets in metres -> world (x, y)."""
    out = []
    for f, r in pts_fr:
        out.append((x0 + f * math.sin(h0) + r * math.sin(h0 + math.pi / 2),
                    y0 + f * math.cos(h0) + r * math.cos(h0 + math.pi / 2)))
    return out


def build_mission(kind, af, x0, y0, h0, rng):
    slow = af.stall_margin * af.v_stall * rng.uniform(1.05, 1.25)
    if kind == "lawnmower":
        L, sp, n = rng.uniform(600, 1500), rng.uniform(100, 250), rng.integers(4, 7)
        pts = []
        for i in range(n):
            ends = (0.0, L) if i % 2 == 0 else (L, 0.0)
            pts += [(ends[0], i * sp), (ends[1], i * sp)]
        return [seg_wp(to_world(pts, x0, y0, h0), v=af.v_cruise * 0.85)]
    if kind == "racetrack":
        L, W = rng.uniform(800, 2000), rng.uniform(150, 400)
        lap = [(L, 0), (L, W), (0, W), (0, 0)]
        return [seg_wp(to_world(lap * 2, x0, y0, h0))]
    if kind == "loiter":
        d = int(rng.choice([-1, 1]))
        return [seg_straight(rng.uniform(500, 1500)),
                seg_loiter(rng.uniform(40, 150), rng.uniform(2, 3.5), d, v=slow),
                seg_straight(rng.uniform(500, 1000))]
    if kind == "evasive":
        d = int(rng.choice([-1, 1]))
        return [seg_straight(rng.uniform(600, 1200)),
                seg_break(d, rng.uniform(90, 270)),
                seg_straight(rng.uniform(400, 800)),
                seg_break(-d, rng.uniform(60, 180)),
                seg_straight(rng.uniform(400, 800))]
    if kind == "mixed":
        segs = [seg_straight(rng.uniform(300, 800))]
        for _ in range(rng.integers(3, 6)):
            c = rng.choice(["break", "loiter", "sturn", "straight"])
            d = int(rng.choice([-1, 1]))
            if c == "break":
                segs.append(seg_break(d, rng.uniform(45, 200)))
            elif c == "loiter":
                segs.append(seg_loiter(rng.uniform(50, 150), rng.uniform(1, 2), d, v=slow))
            elif c == "sturn":
                segs += [seg_break(d, 90), seg_break(-d, 90)]
            segs.append(seg_straight(rng.uniform(200, 600)))
        return segs
    raise ValueError(f"unknown scenario {kind}")


SCENARIOS = ["lawnmower", "racetrack", "loiter", "evasive", "mixed"]


# ----------------------------------------------------------------------------- flight sim
def simulate(mission, af, wind, x0, y0, h0, alt, dt, t_max, rng):
    g = Guidance(mission, af)
    x, y, hdg, v, bank = x0, y0, h0, af.v_cruise, 0.0
    to = (wind.dir_from_deg + 180.0) * D2R           # direction wind blows TO
    w_steady = np.array([wind.speed * math.sin(to), wind.speed * math.cos(to)])
    gust = np.zeros(2)
    roll_rate_max = af.roll_rate_max_deg * D2R
    rec = {k: [] for k in ["t", "x", "y", "alt", "vx", "vy", "gs", "airspeed", "heading_deg",
                           "course_deg", "bank_deg", "turn_rate_dps", "segment"]}
    t = 0.0
    while t <= t_max:
        gust += (-gust / wind.gust_tau) * dt + wind.gust_sigma * math.sqrt(2 * dt / wind.gust_tau) \
            * rng.standard_normal(2)
        w = w_steady + gust
        vx, vy = v * math.sin(hdg) + w[0], v * math.cos(hdg) + w[1]
        gs, course = math.hypot(vx, vy), bearing(vx, vy)

        cmd = g.command(x, y, course, gs, hdg)
        if cmd is None:
            break
        bank_cmd, v_cmd = cmd

        # roll dynamics: first-order lag + rate limit
        bank += clip((bank_cmd - bank) / af.roll_tau, -roll_rate_max, roll_rate_max) * dt
        n = 1.0 / math.cos(bank)                          # load factor
        v_min = af.stall_margin * af.v_stall * math.sqrt(n)
        v_tgt = clip(v_cmd * (1.0 - af.turn_bleed * (n - 1.0)), v_min, af.v_max)
        v += (v_tgt - v) / af.speed_tau * dt
        turn_rate = G * math.tan(bank) / v                # coordinated turn
        hdg = (hdg + turn_rate * dt) % TWO_PI

        for k, val in zip(rec, [t, x, y, alt, vx, vy, gs, v, hdg / D2R % 360, course / D2R % 360,
                                bank / D2R, turn_rate / D2R, g.current_type()]):
            rec[k].append(val)
        x += vx * dt
        y += vy * dt
        t += dt

    out = {k: np.array(vals) for k, vals in rec.items()}
    c = np.unwrap(out["course_deg"] * D2R)
    out["course_rate_dps"] = np.gradient(c, out["t"]) / D2R
    out["lat_accel_g"] = out["gs"] * out["course_rate_dps"] * D2R / G
    out["turning"] = (np.abs(out["bank_deg"]) > TURN_LABEL_BANK_DEG).astype(int)
    return out


# ----------------------------------------------------------------------------- radar
def radar_detect(tr, radar, rng):
    t = tr["t"]
    T = 60.0 / radar.rpm
    ant = (rng.uniform(0, TWO_PI) + TWO_PI * t / T) % TWO_PI
    taz = np.arctan2(tr["x"] - radar.x, tr["y"] - radar.y)
    diff = (taz - ant + math.pi) % TWO_PI - math.pi
    cross = (diff[:-1] > 0) & (diff[1:] <= 0)           # beam sweeps past target
    rows = []
    for scan, i in enumerate(np.nonzero(cross)[0] + 1):
        truth = dict(truth_x=tr["x"][i], truth_y=tr["y"][i], truth_gs=tr["gs"][i],
                     truth_course_deg=tr["course_deg"][i], truth_turning=int(tr["turning"][i]),
                     truth_bank_deg=tr["bank_deg"][i])
        if rng.random() < radar.pd:
            r = math.hypot(tr["x"][i] - radar.x, tr["y"][i] - radar.y) \
                + rng.normal(0, radar.range_sigma)
            az = taz[i] + rng.normal(0, radar.az_sigma_deg * D2R)
            rows.append(dict(t=t[i], scan=scan, range_m=r, az_deg=az / D2R % 360,
                             x=radar.x + r * math.sin(az), y=radar.y + r * math.cos(az),
                             is_false=0, **truth))
        for _ in range(rng.poisson(radar.false_per_scan)):
            fx = tr["x"][i] + rng.uniform(-radar.false_box, radar.false_box)
            fy = tr["y"][i] + rng.uniform(-radar.false_box, radar.false_box)
            r, az = math.hypot(fx - radar.x, fy - radar.y), bearing(fx - radar.x, fy - radar.y)
            rows.append(dict(t=t[i], scan=scan, range_m=r, az_deg=az / D2R % 360, x=fx, y=fy,
                             is_false=1, **truth))
    return rows


# ----------------------------------------------------------------------------- monte carlo
def random_run(rng, args):
    v_stall = rng.uniform(10, 15)
    af = Airframe(v_stall=v_stall,
                  v_cruise=min(rng.uniform(1.4, 1.9) * v_stall, 28.0),
                  v_max=30.0,
                  bank_max_deg=rng.uniform(30, 55),
                  roll_tau=rng.uniform(0.3, 0.8),
                  roll_rate_max_deg=rng.uniform(40, 90))
    wind = Wind(speed=rng.uniform(0, args.max_wind), dir_from_deg=rng.uniform(0, 360),
                gust_sigma=rng.uniform(0.3, 2.0), gust_tau=rng.uniform(2, 5))
    radar = Radar(rpm=args.rpm if args.rpm else float(rng.choice([20, 30])),
                  pd=args.pd, false_per_scan=args.false_per_scan)
    rng_m, az0 = rng.uniform(2000, 12000), rng.uniform(0, TWO_PI)
    x0, y0 = rng_m * math.sin(az0), rng_m * math.cos(az0)
    h0 = rng.uniform(0, TWO_PI)
    kind = args.scenario if args.scenario != "random" else str(rng.choice(SCENARIOS))
    return af, wind, radar, kind, x0, y0, h0, rng.uniform(100, 600)


def write_csv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def plot_run(path, run, tr, dets):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.5))
    a1.plot(tr["x"], tr["y"], "k-", lw=1, label="truth")
    hits = [d for d in dets if not d["is_false"]]
    fa = [d for d in dets if d["is_false"]]
    a1.plot([d["x"] for d in hits], [d["y"] for d in hits], "o", ms=4, label="radar hits")
    if fa:
        a1.plot([d["x"] for d in fa], [d["y"] for d in fa], "x", ms=4, label="false alarms")
    a1.plot(tr["x"][0], tr["y"][0], "g^", ms=9, label="start")
    a1.set_aspect("equal")
    a1.set_xlabel("East (m)")
    a1.set_ylabel("North (m)")
    a1.legend()
    a1.set_title(f"Run {run}")
    a2.plot(tr["t"], tr["gs"], label="ground speed")
    a2.plot(tr["t"], tr["airspeed"], label="airspeed")
    a2.fill_between(tr["t"], 0, tr["gs"].max() * 1.1, where=tr["turning"] == 1,
                    alpha=0.15, color="orange", label="turning (truth)")
    a2.set_xlabel("t (s)")
    a2.set_ylabel("m/s")
    a2.legend()
    a2.set_title("Speed vs time")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", type=int, default=20)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--scenario", default="random", choices=["random"] + SCENARIOS)
    p.add_argument("--rpm", type=float, default=0, help="20 or 30; 0 = random per run")
    p.add_argument("--pd", type=float, default=0.9)
    p.add_argument("--false-per-scan", type=float, default=0.0)
    p.add_argument("--max-wind", type=float, default=10.0)
    p.add_argument("--dt", type=float, default=0.05)
    p.add_argument("--tmax", type=float, default=600.0)
    p.add_argument("--truth-every", type=int, default=2, help="write every Nth truth sample")
    p.add_argument("--out", default="sim_out")
    p.add_argument("--plot", action="store_true", help="save PNG for first 5 runs")
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    truth_rows, det_rows, run_rows = [], [], []

    for run in range(args.runs):
        af, wind, radar, kind, x0, y0, h0, alt = random_run(rng, args)
        mission = build_mission(kind, af, x0, y0, h0, rng)
        tr = simulate(mission, af, wind, x0, y0, h0, alt, args.dt, args.tmax, rng)
        dets = radar_detect(tr, radar, rng)

        keys = [k for k in tr if k != "segment"]
        for i in range(0, len(tr["t"]), args.truth_every):
            row = {"run": run, "segment": tr["segment"][i]}
            row.update({k: round(float(tr[k][i]), 4) for k in keys})
            truth_rows.append(row)
        for d in dets:
            det_rows.append({"run": run, **{k: (round(v, 4) if isinstance(v, float) else v)
                                           for k, v in d.items()}})
        run_rows.append({"run": run, "scenario": kind, "duration_s": round(float(tr["t"][-1]), 1),
                         "start_range_m": round(math.hypot(x0, y0)),
                         **{f"af_{k}": round(v, 3) for k, v in asdict(af).items()},
                         **{f"wind_{k}": round(v, 3) for k, v in asdict(wind).items()},
                         **{f"radar_{k}": v for k, v in asdict(radar).items()}})
        if args.plot and run < 5:
            plot_run(os.path.join(args.out, f"run_{run:02d}.png"), run, tr, dets)
        print(f"run {run:3d}  {kind:10s}  {tr['t'][-1]:6.0f} s  {len(dets):4d} plots  "
              f"wind {wind.speed:4.1f} m/s  gs {tr['gs'].min():4.1f}-{tr['gs'].max():4.1f} m/s  "
              f"turning {100 * tr['turning'].mean():3.0f}%")

    write_csv(os.path.join(args.out, "truth.csv"), truth_rows, list(truth_rows[0].keys()))
    det_fields = ["run", "t", "scan", "range_m", "az_deg", "x", "y", "is_false", "truth_x",
                  "truth_y", "truth_gs", "truth_course_deg", "truth_turning", "truth_bank_deg"]
    write_csv(os.path.join(args.out, "detections.csv"), det_rows, det_fields)
    write_csv(os.path.join(args.out, "runs.csv"), run_rows, list(run_rows[0].keys()))
    print(f"\nWrote {len(truth_rows)} truth rows, {len(det_rows)} detections to ./{args.out}/")


if __name__ == "__main__":
    main()
