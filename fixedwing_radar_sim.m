function [truthAll, detAll, runInfo] = fixedwing_radar_sim(varargin)
%FIXEDWING_RADAR_SIM  Slow / maneuvering fixed-wing drone + scanning radar simulator.
%
%   MATLAB port of fixedwing_radar_sim.py. Constant altitude (z fixed per run, vz = 0).
%
%   Layers:
%     1. Coordinated-turn dynamics (bank -> turn rate), roll lag + rate limit,
%        stall floor that rises with bank, speed bleed in tight turns.
%     2. Autopilot guidance (L1-style): waypoints, loiter orbits, break turns, straight legs.
%     3. Wind: steady + Gauss-Markov gusts (radar sees GROUND velocity).
%     4. Radar: rotating beam (hit when beam crosses target azimuth), range/az noise,
%        Pd < 1, optional false alarms.
%     5. Monte Carlo: random airframe, wind, scenario, start range per run.
%
%   Conventions: x = East, y = North (m), radar at (0,0).
%   Heading / course / azimuth clockwise from North. Positive bank = right turn.
%
%   Usage:
%     [truth, dets, runs] = fixedwing_radar_sim();                 % defaults
%     fixedwing_radar_sim('Runs',50,'RPM',30,'Plot',true);
%     fixedwing_radar_sim('Scenario','loiter','TruthRate',1);      % 1 Hz truth
%
%   Name-value options (defaults):
%     'Runs' 20, 'Seed' 1, 'Scenario' 'random' (lawnmower|racetrack|loiter|evasive|mixed),
%     'RPM' 0 (0 = random 20/30 per run), 'Pd' 0.9, 'FalsePerScan' 0, 'MaxWind' 10,
%     'Dt' 0.05 (sim step, s), 'TMax' 600, 'TruthRate' 10 (Hz), 'OutDir' 'sim_out',
%     'Plot' false (saves PNGs for first 5 runs)
%
%   Outputs (also written as CSV to OutDir):
%     truthAll : run, t, x, y, z, vx, vy, vz, gs, airspeed, heading_deg, course_deg,
%                bank_deg, turn_rate_dps, course_rate_dps, lat_accel_g, turning, segment
%     detAll   : radar plots with truth attached (truth_turning etc.) for scoring
%     runInfo  : parameters for each run

p = inputParser;
addParameter(p, 'Runs', 20);
addParameter(p, 'Seed', 1);
addParameter(p, 'Scenario', 'random');
addParameter(p, 'RPM', 0);
addParameter(p, 'Pd', 0.9);
addParameter(p, 'FalsePerScan', 0);
addParameter(p, 'MaxWind', 10);
addParameter(p, 'Dt', 0.05);
addParameter(p, 'TMax', 600);
addParameter(p, 'TruthRate', 10);
addParameter(p, 'OutDir', 'sim_out');
addParameter(p, 'Plot', false);
parse(p, varargin{:});
o = p.Results;

rng(o.Seed);
scenarios = {'lawnmower', 'racetrack', 'loiter', 'evasive', 'mixed'};
if ~exist(o.OutDir, 'dir'), mkdir(o.OutDir); end
truthStep = max(1, round(1 / (o.TruthRate * o.Dt)));

truthParts = cell(o.Runs, 1);
detParts   = cell(o.Runs, 1);
runRows    = cell(o.Runs, 1);

for run = 0:o.Runs-1
    [af, wind, radar, kind, x0, y0, h0, alt] = randomRun(o, scenarios);
    mission = buildMission(kind, af, x0, y0, h0);
    tr  = simulateFlight(mission, af, wind, x0, y0, h0, alt, o.Dt, o.TMax);
    det = radarDetect(tr, radar);

    idx = (1:truthStep:numel(tr.t))';
    n = numel(idx);
    truthParts{run+1} = table(repmat(run, n, 1), tr.t(idx), tr.x(idx), tr.y(idx), ...
        tr.z(idx), tr.vx(idx), tr.vy(idx), tr.vz(idx), tr.gs(idx), tr.airspeed(idx), ...
        tr.heading(idx), tr.course(idx), tr.bank(idx), tr.turnRate(idx), ...
        tr.courseRate(idx), tr.latAccelG(idx), tr.turning(idx), tr.segment(idx), ...
        'VariableNames', {'run','t','x','y','z','vx','vy','vz','gs','airspeed', ...
        'heading_deg','course_deg','bank_deg','turn_rate_dps','course_rate_dps', ...
        'lat_accel_g','turning','segment'});

    det.run = repmat(run, height(det), 1);
    detParts{run+1} = det(:, [end, 1:end-1]);

    runRows{run+1} = struct('run', run, 'scenario', kind, 'duration_s', tr.t(end), ...
        'start_range_m', hypot(x0, y0), 'z_m', alt, ...
        'v_stall', af.vStall, 'v_cruise', af.vCruise, 'bank_max_deg', af.bankMaxDeg, ...
        'roll_tau', af.rollTau, 'roll_rate_max_deg', af.rollRateMaxDeg, ...
        'wind_speed', wind.speed, 'wind_dir_from_deg', wind.dirFromDeg, ...
        'gust_sigma', wind.gustSigma, 'gust_tau', wind.gustTau, ...
        'radar_rpm', radar.rpm, 'radar_pd', radar.pd, 'false_per_scan', radar.falsePerScan);

    if o.Plot && run < 5
        plotRun(fullfile(o.OutDir, sprintf('run_%02d.png', run)), run, tr, det);
    end
    fprintf('run %3d  %-10s %6.0f s  %4d plots  wind %4.1f m/s  gs %4.1f-%4.1f m/s  turning %3.0f%%\n', ...
        run, kind, tr.t(end), height(det), wind.speed, min(tr.gs), max(tr.gs), 100*mean(tr.turning));
end

truthAll = vertcat(truthParts{:});
detAll   = vertcat(detParts{:});
runInfo  = struct2table([runRows{:}]);
writetable(truthAll, fullfile(o.OutDir, 'truth.csv'));
writetable(detAll,   fullfile(o.OutDir, 'detections.csv'));
writetable(runInfo,  fullfile(o.OutDir, 'runs.csv'));
fprintf('\nWrote %d truth rows (%.3g Hz), %d detections to %s\n', ...
    height(truthAll), 1/(truthStep*o.Dt), height(detAll), o.OutDir);
end

% =========================================================================== random run setup
function [af, wind, radar, kind, x0, y0, h0, alt] = randomRun(o, scenarios)
vStall = U(10, 15);
af = struct('vStall', vStall, ...
    'vCruise', min(U(1.4, 1.9) * vStall, 28), ...
    'vMax', 30, ...
    'bankMaxDeg', U(30, 55), ...
    'rollTau', U(0.3, 0.8), ...
    'rollRateMaxDeg', U(40, 90), ...
    'speedTau', 3, ...
    'turnBleed', 0.3, ...
    'stallMargin', 1.2);
wind = struct('speed', U(0, o.MaxWind), 'dirFromDeg', U(0, 360), ...
    'gustSigma', U(0.3, 2.0), 'gustTau', U(2, 5));
if o.RPM > 0
    rpm = o.RPM;
else
    rpms = [20 30];
    rpm = rpms(randi(2));
end
radar = struct('x', 0, 'y', 0, 'rpm', rpm, 'rangeSigma', 15, 'azSigmaDeg', 0.3, ...
    'pd', o.Pd, 'falsePerScan', o.FalsePerScan, 'falseBox', 500);
r0 = U(2000, 12000);
az0 = U(0, 2*pi);
x0 = r0 * sin(az0);
y0 = r0 * cos(az0);
h0 = U(0, 2*pi);
if strcmpi(o.Scenario, 'random')
    kind = scenarios{randi(numel(scenarios))};
else
    kind = lower(o.Scenario);
end
alt = U(100, 600);
end

% =========================================================================== missions
function mission = buildMission(kind, af, x0, y0, h0)
slow = af.stallMargin * af.vStall * U(1.05, 1.25);
switch kind
    case 'lawnmower'
        L = U(600, 1500); sp = U(100, 250); n = randi([4 6]);
        pts = zeros(2*n, 2);
        for i = 0:n-1
            if mod(i, 2) == 0, ends = [0 L]; else, ends = [L 0]; end
            pts(2*i+1, :) = [ends(1), i*sp];
            pts(2*i+2, :) = [ends(2), i*sp];
        end
        mission = {mkWp(toWorld(pts, x0, y0, h0), af.vCruise * 0.85, [])};
    case 'racetrack'
        L = U(800, 2000); W = U(150, 400);
        lap = [L 0; L W; 0 W; 0 0];
        mission = {mkWp(toWorld([lap; lap], x0, y0, h0), [], [])};
    case 'loiter'
        d = randSign();
        mission = {mkStraight(U(500, 1500), []), ...
                   mkLoiter(U(40, 150), U(2, 3.5), d, [], slow), ...
                   mkStraight(U(500, 1000), [])};
    case 'evasive'
        d = randSign();
        mission = {mkStraight(U(600, 1200), []), ...
                   mkBreak(d, U(90, 270), []), ...
                   mkStraight(U(400, 800), []), ...
                   mkBreak(-d, U(60, 180), []), ...
                   mkStraight(U(400, 800), [])};
    case 'mixed'
        mission = {mkStraight(U(300, 800), [])};
        choices = {'break', 'loiter', 'sturn', 'straight'};
        for k = 1:randi([3 5])
            c = choices{randi(4)};
            d = randSign();
            switch c
                case 'break'
                    mission{end+1} = mkBreak(d, U(45, 200), []); %#ok<AGROW>
                case 'loiter'
                    mission{end+1} = mkLoiter(U(50, 150), U(1, 2), d, [], slow); %#ok<AGROW>
                case 'sturn'
                    mission{end+1} = mkBreak(d, 90, []); %#ok<AGROW>
                    mission{end+1} = mkBreak(-d, 90, []); %#ok<AGROW>
            end
            mission{end+1} = mkStraight(U(200, 600), []); %#ok<AGROW>
        end
    otherwise
        error('Unknown scenario: %s', kind);
end
end

function w = toWorld(ptsFR, x0, y0, h0)
% (forward, right) offsets in metres -> world [x y]
f = ptsFR(:, 1); r = ptsFR(:, 2);
w = [x0 + f*sin(h0) + r*sin(h0 + pi/2), y0 + f*cos(h0) + r*cos(h0 + pi/2)];
end

function s = mkWp(pts, v, accept)
s = struct('type', 'wp', 'pts', pts, 'v', v, 'accept', accept);
end

function s = mkStraight(dist, v)
s = struct('type', 'straight', 'dist', dist, 'v', v);
end

function s = mkLoiter(radius, turns, dir, center, v)
% center = [] -> orbit starts where the aircraft is (circle beside current course)
s = struct('type', 'loiter', 'radius', radius, 'turns', turns, 'dir', dir, ...
    'center', center, 'v', v);
end

function s = mkBreak(dir, angleDeg, v)
s = struct('type', 'break', 'dir', dir, 'angle', angleDeg * pi/180, 'v', v);
end

% =========================================================================== flight sim
function tr = simulateFlight(mission, af, wind, x0, y0, h0, alt, dt, tmax)
G = 9.81; D2R = pi/180;
N = floor(tmax / dt) + 1;
t = zeros(N,1); X = t; Y = t; VX = t; VY = t; GS = t; AS = t;
HD = t; CR = t; BK = t; TRt = t; SEG = cell(N,1);

g = struct('mission', {mission}, 'i', 1, 's', struct(), 'af', af, 'L1', 60, ...
    'bankMax', af.bankMaxDeg * D2R);
x = x0; y = y0; hdg = h0; v = af.vCruise; bank = 0;
to = (wind.dirFromDeg + 180) * D2R;              % direction wind blows TO
wS = wind.speed * [sin(to), cos(to)];
gust = [0 0];
rr = af.rollRateMaxDeg * D2R;
n = 0;

for k = 0:N-1
    tk = k * dt;
    gust = gust + (-gust / wind.gustTau) * dt + ...
        wind.gustSigma * sqrt(2*dt / wind.gustTau) * randn(1, 2);
    w = wS + gust;
    vx = v*sin(hdg) + w(1);
    vy = v*cos(hdg) + w(2);
    gs = hypot(vx, vy);
    course = atan2(vx, vy);

    [g, bankCmd, vCmd, ok] = guidanceCommand(g, x, y, course, gs, hdg);
    if ~ok, break; end

    % roll dynamics: first-order lag + rate limit
    bank = bank + clamp((bankCmd - bank) / af.rollTau, -rr, rr) * dt;
    nf = 1 / cos(bank);                                    % load factor
    vMin = af.stallMargin * af.vStall * sqrt(nf);
    vT = clamp(vCmd * (1 - af.turnBleed * (nf - 1)), vMin, af.vMax);
    v = v + (vT - v) / af.speedTau * dt;
    turnRate = G * tan(bank) / v;                          % coordinated turn
    hdg = mod(hdg + turnRate * dt, 2*pi);

    n = n + 1;
    t(n) = tk; X(n) = x; Y(n) = y; VX(n) = vx; VY(n) = vy; GS(n) = gs; AS(n) = v;
    HD(n) = mod(hdg / D2R, 360); CR(n) = mod(course / D2R, 360);
    BK(n) = bank / D2R; TRt(n) = turnRate / D2R;
    SEG{n} = g.mission{min(g.i, numel(g.mission))}.type;

    x = x + vx * dt;
    y = y + vy * dt;
end

r = 1:n;
tr.t = t(r); tr.x = X(r); tr.y = Y(r);
tr.z = repmat(alt, n, 1); tr.vz = zeros(n, 1);
tr.vx = VX(r); tr.vy = VY(r); tr.gs = GS(r); tr.airspeed = AS(r);
tr.heading = HD(r); tr.course = CR(r); tr.bank = BK(r); tr.turnRate = TRt(r);
tr.segment = SEG(r);
if n > 1
    tr.courseRate = gradient(unwrap(tr.course * D2R), tr.t) / D2R;
else
    tr.courseRate = 0;
end
tr.latAccelG = tr.gs .* tr.courseRate * D2R / G;
tr.turning = double(abs(tr.bank) > 10);                    % truth label
end

% =========================================================================== guidance
function [g, bank, v, ok] = guidanceCommand(g, x, y, course, gs, hdg)
bank = 0; v = g.af.vCruise; ok = false;
while g.i <= numel(g.mission)
    seg = g.mission{g.i};
    switch seg.type
        case 'wp'
            [g.s, bank, v, active] = gWp(g.s, seg, x, y, course, gs, g.af, g.L1, g.bankMax);
        case 'straight'
            [g.s, bank, v, active] = gStraight(g.s, seg, x, y, course, gs, g.af, g.L1, g.bankMax);
        case 'loiter'
            [g.s, bank, v, active] = gLoiter(g.s, seg, x, y, course, gs, g.af, g.L1, g.bankMax);
        case 'break'
            [g.s, bank, v, active] = gBreak(g.s, seg, hdg, g.af, g.bankMax);
    end
    if active, ok = true; return; end
    g.i = g.i + 1;
    g.s = struct();
end
end

function [s, bank, v, active] = gWp(s, seg, x, y, course, gs, af, L1, bankMax)
G = 9.81;
bank = 0; v = pickV(seg, af); active = false;
if ~isfield(s, 'prev'), s.prev = [x y]; s.k = 1; end
pts = seg.pts;
while s.k <= size(pts, 1)
    a = s.prev; b = pts(s.k, :);
    L = max(hypot(b(1)-a(1), b(2)-a(2)), 1e-6);
    u = (b - a) / L;
    along = (x - a(1))*u(1) + (y - a(2))*u(2);
    rTurn = gs^2 / (G * tan(bankMax));
    if isempty(seg.accept), acc = max(30, rTurn); else, acc = seg.accept; end
    if hypot(b(1)-x, b(2)-y) < acc || along >= L
        s.prev = b; s.k = s.k + 1;
        continue
    end
    pt = a + u * min(along + L1, L);
    bank = l1Bank(course, gs, atan2(pt(1)-x, pt(2)-y), hypot(pt(1)-x, pt(2)-y), bankMax);
    active = true;
    return
end
end

function [s, bank, v, active] = gStraight(s, seg, x, y, course, gs, af, L1, bankMax)
if ~isfield(s, 'wpseg')
    s.wpseg = mkWp([x + seg.dist*sin(course), y + seg.dist*cos(course)], seg.v, []);
end
[s, bank, v, active] = gWp(s, s.wpseg, x, y, course, gs, af, L1, bankMax);
end

function [s, bank, v, active] = gLoiter(s, seg, x, y, course, gs, af, L1, bankMax)
bank = 0; v = pickV(seg, af); active = false;
R = seg.radius; d = seg.dir;
if ~isfield(s, 'c')
    if ~isempty(seg.center)
        s.c = seg.center;
    else
        side = course + d*pi/2;
        s.c = [x + R*sin(side), y + R*cos(side)];
    end
    s.acc = NaN; s.th = NaN;
end
dist = hypot(x - s.c(1), y - s.c(2));
th = atan2(x - s.c(1), y - s.c(2));
if isnan(s.acc)
    if abs(dist - R) < 0.25*R, s.acc = 0; s.th = th; end   % captured on circle
else
    s.acc = s.acc + d * wrapPi(th - s.th);
    s.th = th;
    if s.acc >= 2*pi*seg.turns, return; end
end
phi = clamp(L1 / R, 0.3, 1.5);
tt = th + d*phi;
pt = s.c + R*[sin(tt), cos(tt)];
bank = l1Bank(course, gs, atan2(pt(1)-x, pt(2)-y), hypot(pt(1)-x, pt(2)-y), bankMax);
active = true;
end

function [s, bank, v, active] = gBreak(s, seg, hdg, af, bankMax)
bank = 0; v = pickV(seg, af); active = false;
if ~isfield(s, 'acc'), s.acc = 0; s.h = hdg; end
s.acc = s.acc + seg.dir * wrapPi(hdg - s.h);
s.h = hdg;
if s.acc >= seg.angle, return; end
bank = seg.dir * bankMax;
active = true;
end

function bank = l1Bank(course, gs, brg, dist, bankMax)
% L1 guidance: lateral accel = 2 V^2 sin(eta) / L1 -> bank angle
G = 9.81;
eta = clamp(wrapPi(brg - course), -pi/2, pi/2);
aLat = 2 * gs^2 * sin(eta) / max(dist, 1);
bank = clamp(atan2(aLat, G), -bankMax, bankMax);
end

% =========================================================================== radar
function det = radarDetect(tr, radar)
D2R = pi/180;
T = 60 / radar.rpm;
ant = mod(2*pi*rand + 2*pi*tr.t / T, 2*pi);
taz = atan2(tr.x - radar.x, tr.y - radar.y);
d = mod(taz - ant + pi, 2*pi) - pi;
hits = find(d(1:end-1) > 0 & d(2:end) <= 0) + 1;          % beam sweeps past target

maxRows = numel(hits) * (2 + ceil(5*radar.falsePerScan));
R = nan(maxRows, 14); m = 0;
for scan = 0:numel(hits)-1
    i = hits(scan+1);
    truth = [tr.x(i), tr.y(i), tr.gs(i), tr.course(i), tr.turning(i), tr.bank(i)];
    if rand < radar.pd
        r  = hypot(tr.x(i) - radar.x, tr.y(i) - radar.y) + radar.rangeSigma*randn;
        az = taz(i) + radar.azSigmaDeg*D2R*randn;
        m = m + 1;
        R(m, :) = [tr.t(i), scan, r, mod(az/D2R, 360), radar.x + r*sin(az), ...
                   radar.y + r*cos(az), 0, truth];
    end
    for f = 1:poissonKnuth(radar.falsePerScan)
        fx = tr.x(i) + U(-radar.falseBox, radar.falseBox);
        fy = tr.y(i) + U(-radar.falseBox, radar.falseBox);
        r = hypot(fx - radar.x, fy - radar.y);
        az = atan2(fx - radar.x, fy - radar.y);
        m = m + 1;
        R(m, :) = [tr.t(i), scan, r, mod(az/D2R, 360), fx, fy, 1, truth];
    end
end
det = array2table(R(1:m, :), 'VariableNames', {'t','scan','range_m','az_deg','x','y', ...
    'is_false','truth_x','truth_y','truth_gs','truth_course_deg','truth_turning','truth_bank_deg'});
end

% =========================================================================== plotting
function plotRun(fname, run, tr, det)
fig = figure('Visible', 'off', 'Position', [100 100 1300 520]);
subplot(1, 2, 1); hold on; grid on; axis equal;
plot(tr.x, tr.y, 'k-', 'LineWidth', 1);
h = det.is_false == 0;
plot(det.x(h), det.y(h), 'o', 'MarkerSize', 4);
if any(~h), plot(det.x(~h), det.y(~h), 'x', 'MarkerSize', 5); end
plot(tr.x(1), tr.y(1), 'g^', 'MarkerSize', 9, 'MarkerFaceColor', 'g');
xlabel('East (m)'); ylabel('North (m)'); title(sprintf('Run %d', run));
if any(~h)
    legend('truth', 'radar hits', 'false alarms', 'start', 'Location', 'best');
else
    legend('truth', 'radar hits', 'start', 'Location', 'best');
end
subplot(1, 2, 2); hold on; grid on;
top = max(tr.gs) * 1.1;
area(tr.t, tr.turning * top, 'FaceColor', [1 0.6 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(tr.t, tr.gs); plot(tr.t, tr.airspeed);
xlabel('t (s)'); ylabel('m/s'); title('Speed vs time');
legend('turning (truth)', 'ground speed', 'airspeed', 'Location', 'best');
saveas(fig, fname);
close(fig);
end

% =========================================================================== helpers
function v = pickV(seg, af)
if isempty(seg.v), v = af.vCruise; else, v = seg.v; end
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

function a = wrapPi(a)
a = mod(a + pi, 2*pi) - pi;
end

function x = U(a, b)
x = a + (b - a) * rand;
end

function d = randSign()
d = 2*randi([0 1]) - 1;
end

function k = poissonKnuth(lambda)
% Poisson sample without Statistics Toolbox
k = 0;
if lambda <= 0, return; end
L = exp(-lambda); p = 1;
while true
    p = p * rand;
    if p <= L, return; end
    k = k + 1;
end
end
