function tables = writeRequirementReport(fileName, requirements, summaryTable)
%WRITEREQUIREMENTREPORT Export differently shaped requirement tables to Excel.
%   tables = writeRequirementReport(fileName, requirements)
%   tables = writeRequirementReport(fileName, requirements, summaryTable)
%
%   Each element of requirements is a struct with these fields:
%     ID       - Requirement ID, e.g. 'REQ-001'
%     Headers  - One column name per data vector (cell array or string array)
%     Data     - Cell array containing real numeric/logical vectors
%
%   Optional summaryTable holds your ALREADY COMPUTED summary statistics.
%   It is written first, followed by a blank row and each requirement table.
%   Statistics and pass/fail criteria are not calculated by this formatter.
%
%   Each requirement can have a different number of columns and rows.
%   Shorter vectors are padded with NaN, exported as blank Excel cells.
%   Padding is for display only: it does not align or match test/model samples.
%   Scalars occupy one row and are not repeated to fill a longer vector.
%   Empty vectors become blank columns. All-empty data gets one ID-only row.
%
%   The 'Results' worksheet is replaced on each call. Other sheets are kept.
%   Requires a MATLAB release supporting writetable WriteMode/overwritesheet.
%
%   EXAMPLE (illustrative data -- replace with your own variables):
%   req(1).ID = 'REQ-001';
%   req(1).Headers = {'UT_Doppler_Index', 'Model_Doppler_Index'};
%   req(1).Data = {[12; 18; 24], [12; 18; 24]};
%
%   req(2).ID = 'REQ-002';
%   req(2).Headers = {'UT_Doppler_Index', 'Model_Doppler_Index', ...
%                     'UT_Range_Bin_Index', 'Model_Range_Bin_Index'};
%   req(2).Data = {[12; 18], [12; 18], [100; 110], [100; 110]};
%
%   req(3).ID = 'REQ-003';
%   req(3).Headers = {'UT_SLL', 'Model_SLL', 'UT_NBar', 'Model_NBar'};
%   req(3).Data = {[-30; -31], [-30; -31], [4; 4], [4; 4]};
%
%   req(4).ID = 'REQ-004';
%   req(4).Headers = {'UT_Peaks_Survive_Match', ...
%       'Model_Peaks_Survive_Match', 'UT_Final_TOI', 'Model_Final_TOI'};
%   req(4).Data = {[1; 1; 0], [1; 1; 0], [5; 6; 7], [5; 6; 7]};
%
%   req(5).ID = 'REQ-005';
%   req(5).Headers = {'UT_Peak_To_Valley_Indices', ...
%                     'Model_Peak_To_Valley_Indices'};
%   req(5).Data = {[4; 8; 12], [4; 8]}; % Last model cell will be blank.
%
%   tables = writeRequirementReport('Requirement_Report.xlsx', req);
%   disp(tables{1});
%
%   OPTIONAL SUMMARY: use existing results, one row per requirement/metric.
%   When a requirement has several metrics (e.g. Doppler and range), include
%   a separate row per metric so statistics with different units stay separate.
%   All variables below must contain the same number of summary rows:
%
%   summaryTable = table(reqIDs, metrics, nSamples, meanError, mae, rmse, ...
%       sigma, threeSigma, maxAbsError, tolerance, pctWithin, result, ...
%       'VariableNames', {'Requirement_ID', 'Metric', 'N_Samples', ...
%       'Mean_Error', 'MAE', 'RMSE', 'Sigma', 'Three_Sigma', ...
%       'Max_Absolute_Error', 'Tolerance', 'Percent_Within_Tolerance', ...
%       'Result'});
%   tables = writeRequirementReport('Requirement_Report.xlsx', ...
%                                   req, summaryTable);

    if nargin < 3
        summaryTable = table();
    end
    if ~(ischar(fileName) && isrow(fileName)) && ...
            ~(isstring(fileName) && isscalar(fileName) && ~ismissing(fileName))
        error('fileName must be a character row or a string scalar.');
    end
    [~, ~, extension] = fileparts(fileName);
    if ~strcmpi(extension, '.xlsx')
        error('Use an .xlsx output filename.');
    end
    if ~isstruct(requirements) || isempty(requirements) || ...
            ~all(isfield(requirements, {'ID', 'Headers', 'Data'}))
        error('Supply a nonempty struct array with ID, Headers, and Data.');
    end
    if ~istable(summaryTable)
        error('summaryTable must be a MATLAB table.');
    end

    % Build and validate all tables before changing the output worksheet.
    tables = cell(numel(requirements), 1);
    for k = 1:numel(requirements)
        req = requirements(k);
        id = string(req.ID);
        if ~isscalar(id) || ismissing(id) || strlength(id) == 0
            error('Requirement %d must have one nonempty ID.', k);
        end
        headers = string(req.Headers);
        headers = headers(:).';
        columns = req.Data;
        if ~iscell(columns) || isempty(columns) || ...
                numel(headers) ~= numel(columns)
            error('Requirement %d needs one header per vector in Data.', k);
        end
        if any(ismissing(headers)) || any(strlength(headers) == 0) || ...
                numel(unique(headers)) ~= numel(headers) || ...
                any(headers == "Requirement_ID")
            error('Requirement %d has empty, duplicate, or reserved headers.', k);
        end

        nRows = max(1, max(cellfun(@numel, columns)));
        T = table(repmat(id, nRows, 1), ...
                  'VariableNames', {'Requirement_ID'});
        for j = 1:numel(columns)
            values = columns{j};
            if ~(isnumeric(values) || islogical(values)) || ...
                    ~isreal(values) || ~(isvector(values) || isempty(values))
                error('Requirement %d column %d must be a real numeric vector.', k, j);
            end
            padded = NaN(nRows, 1);
            padded(1:numel(values)) = values(:);
            T.(char(headers(j))) = padded;
        end
        tables{k} = T;
    end

    % Keep the worksheet within Excel's row/column limits.
    blocks = tables;
    if width(summaryTable) > 0
        blocks = [{summaryTable}; tables];
    end
    rowsNeeded = sum(cellfun(@height, blocks)) + 2*numel(blocks) - 1;
    if rowsNeeded > 1048576 || any(cellfun(@width, blocks) > 16384)
        error('The report exceeds the Excel worksheet size limit.');
    end

    % Clear Results only on the first write, avoiding stale rows on reruns.
    nextRow = 1;
    for k = 1:numel(blocks)
        if k == 1
            mode = 'overwritesheet';
        else
            mode = 'inplace';
        end
        writetable(blocks{k}, fileName, 'Sheet', 'Results', ...
            'Range', sprintf('A%d', nextRow), ...
            'WriteVariableNames', true, 'WriteRowNames', false, ...
            'WriteMode', mode);
        % One header row + all data rows + one blank separating row.
        nextRow = nextRow + height(blocks{k}) + 2;
    end
end
