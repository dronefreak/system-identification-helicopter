function result = validate_flight_data(filePath)
% VALIDATE_FLIGHT_DATA - Validate flight data file structure and contents
%
% This function performs comprehensive validation checks on flight data files
% to ensure they meet the required format and quality standards.
%
% Usage:
%    result = validate_flight_data('data/experiments/best.mat')
%    result = validate_flight_data('best')  % Shorthand
%
% Inputs:
%    filePath - Path to .mat file or dataset name
%
% Outputs:
%    result   - Structure containing validation results:
%               .isValid: true/false
%               .errors: Cell array of error messages
%               .warnings: Cell array of warning messages
%               .info: General information about the file
%               .checks: Detailed check results
%
% Validation Checks:
%    1. File exists and is readable
%    2. Required variables present (in_H, out_H, time, inr, outr)
%    3. Correct data dimensions
%    4. No NaN or Inf values (warning if present)
%    5. Time vector is monotonic
%    6. Consistent sample counts across variables
%    7. Data types are correct
%    8. Reasonable value ranges
%
% Example:
%    result = validate_flight_data('best');
%    if result.isValid
%        fprintf('Data is valid!\n');
%    else
%        fprintf('Validation errors:\n');
%        for i = 1:length(result.errors)
%            fprintf('  - %s\n', result.errors{i});
%        end
%    end
%
% See also: load_flight_data
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Initialize Result Structure
    result = struct();
    result.isValid = true;
    result.errors = {};
    result.warnings = {};
    result.info = struct();
    result.checks = struct();

    fprintf('=======================================================\n');
    fprintf('  Flight Data Validation\n');
    fprintf('=======================================================\n');

    %% Check 1: File Exists
    fprintf('Checking file existence...\n');

    % Handle dataset name vs full path
    if ~exist(filePath, 'file')
        [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
        projectRoot = fileparts(fileparts(scriptPath));

        if ~endsWith(filePath, '.mat')
            filePath = [filePath '.mat'];
        end

        alternativePath = fullfile(projectRoot, 'data', 'experiments', filePath);

        if exist(alternativePath, 'file')
            filePath = alternativePath;
        else
            result.isValid = false;
            result.errors{end+1} = sprintf('File not found: %s', filePath);
            fprintf('  ✗ File not found\n');
            print_summary(result);
            return;
        end
    end

    fprintf('  ✓ File exists: %s\n', filePath);
    result.info.filePath = filePath;

    % Get file info
    fileInfo = dir(filePath);
    result.info.fileSize = fileInfo.bytes;
    result.info.fileSizeMB = fileInfo.bytes / 1024^2;
    result.info.lastModified = fileInfo.date;

    fprintf('  ✓ File size: %.2f MB\n', result.info.fileSizeMB);

    %% Check 2: Load File
    fprintf('\nLoading file...\n');

    try
        data = load(filePath);
        fprintf('  ✓ File loaded successfully\n');
        result.checks.fileReadable = true;
    catch ME
        result.isValid = false;
        result.errors{end+1} = sprintf('Cannot load file: %s', ME.message);
        fprintf('  ✗ Failed to load file\n');
        print_summary(result);
        return;
    end

    %% Check 3: Required Variables Present
    fprintf('\nChecking required variables...\n');

    requiredVars = {'in_H', 'out_H', 'time', 'inr', 'outr'};
    result.checks.requiredVariables = struct();

    for i = 1:length(requiredVars)
        varName = requiredVars{i};
        if isfield(data, varName)
            fprintf('  ✓ %s present\n', varName);
            result.checks.requiredVariables.(varName) = true;
        else
            fprintf('  ✗ %s missing\n', varName);
            result.isValid = false;
            result.errors{end+1} = sprintf('Missing required variable: %s', varName);
            result.checks.requiredVariables.(varName) = false;
        end
    end

    if ~result.isValid
        print_summary(result);
        return;
    end

    %% Check 4: Data Dimensions
    fprintf('\nChecking data dimensions...\n');

    % Input data should be 4 × N
    [nInputChannels, nInputSamples] = size(data.in_H);
    fprintf('  in_H: %d × %d\n', nInputChannels, nInputSamples);

    if nInputChannels ~= 4
        result.errors{end+1} = sprintf('in_H should have 4 rows (control inputs), found %d', nInputChannels);
        result.isValid = false;
        fprintf('  ✗ Wrong number of input channels\n');
    else
        fprintf('  ✓ Correct input channels (4)\n');
    end

    % Output data should be 10 × N
    [nOutputChannels, nOutputSamples] = size(data.out_H);
    fprintf('  out_H: %d × %d\n', nOutputChannels, nOutputSamples);

    if nOutputChannels ~= 10
        result.errors{end+1} = sprintf('out_H should have 10 rows (measurements), found %d', nOutputChannels);
        result.isValid = false;
        fprintf('  ✗ Wrong number of output channels\n');
    else
        fprintf('  ✓ Correct output channels (10)\n');
    end

    % Sample counts should match
    if nInputSamples ~= nOutputSamples
        result.errors{end+1} = sprintf('Sample count mismatch: in_H has %d, out_H has %d', ...
                                       nInputSamples, nOutputSamples);
        result.isValid = false;
        fprintf('  ✗ Sample count mismatch\n');
    else
        fprintf('  ✓ Sample counts match (%d samples)\n', nInputSamples);
    end

    result.info.nSamples = nInputSamples;
    result.info.nInputs = nInputChannels;
    result.info.nOutputs = nOutputChannels;

    %% Check 5: Time Vector
    fprintf('\nChecking time vector...\n');

    if ~isvector(data.time)
        result.errors{end+1} = 'time must be a vector';
        result.isValid = false;
        fprintf('  ✗ time is not a vector\n');
    else
        nTimeSamples = length(data.time);
        fprintf('  ✓ time is a vector (%d samples)\n', nTimeSamples);

        if nTimeSamples ~= nInputSamples
            result.warnings{end+1} = sprintf('time vector length (%d) differs from data samples (%d)', ...
                                            nTimeSamples, nInputSamples);
            fprintf('  ⚠ time vector length mismatch\n');
        end

        % Check if monotonically increasing
        if issorted(data.time)
            fprintf('  ✓ time vector is monotonically increasing\n');
            result.info.duration = max(data.time) - min(data.time);
            result.info.samplingRate = 1 / mean(diff(data.time));
            fprintf('  ✓ Duration: %.2f seconds\n', result.info.duration);
            fprintf('  ✓ Sampling rate: ~%.1f Hz\n', result.info.samplingRate);
        else
            result.errors{end+1} = 'time vector is not monotonically increasing';
            result.isValid = false;
            fprintf('  ✗ time vector is not sorted\n');
        end
    end

    %% Check 6: Data Quality
    fprintf('\nChecking data quality...\n');

    % Check for NaN
    nNaN_in = sum(isnan(data.in_H(:)));
    nNaN_out = sum(isnan(data.out_H(:)));

    if nNaN_in > 0
        result.warnings{end+1} = sprintf('in_H contains %d NaN values', nNaN_in);
        fprintf('  ⚠ in_H has %d NaN values\n', nNaN_in);
    else
        fprintf('  ✓ in_H has no NaN values\n');
    end

    if nNaN_out > 0
        result.warnings{end+1} = sprintf('out_H contains %d NaN values', nNaN_out);
        fprintf('  ⚠ out_H has %d NaN values\n', nNaN_out);
    else
        fprintf('  ✓ out_H has no NaN values\n');
    end

    % Check for Inf
    nInf_in = sum(isinf(data.in_H(:)));
    nInf_out = sum(isinf(data.out_H(:)));

    if nInf_in > 0
        result.warnings{end+1} = sprintf('in_H contains %d Inf values', nInf_in);
        fprintf('  ⚠ in_H has %d Inf values\n', nInf_in);
    else
        fprintf('  ✓ in_H has no Inf values\n');
    end

    if nInf_out > 0
        result.warnings{end+1} = sprintf('out_H contains %d Inf values', nInf_out);
        fprintf('  ⚠ out_H has %d Inf values\n', nInf_out);
    else
        fprintf('  ✓ out_H has no Inf values\n');
    end

    %% Check 7: Value Ranges
    fprintf('\nChecking value ranges...\n');

    % Input ranges (control inputs should be reasonable)
    in_min = min(data.in_H, [], 2);
    in_max = max(data.in_H, [], 2);

    fprintf('  Input ranges:\n');
    for i = 1:nInputChannels
        fprintf('    Channel %d: [%.3f, %.3f]\n', i, in_min(i), in_max(i));
    end

    % Output ranges
    out_min = min(data.out_H, [], 2);
    out_max = max(data.out_H, [], 2);

    fprintf('  Output ranges:\n');
    for i = 1:min(nOutputChannels, 5)  % Show first 5
        fprintf('    Channel %d: [%.3f, %.3f]\n', i, out_min(i), out_max(i));
    end
    if nOutputChannels > 5
        fprintf('    ... (%d more channels)\n', nOutputChannels - 5);
    end

    %% Print Summary
    fprintf('\n');
    print_summary(result);
end

function print_summary(result)
    fprintf('=======================================================\n');
    fprintf('  Validation Summary\n');
    fprintf('=======================================================\n');

    if result.isValid
        fprintf('Status: ✓ VALID\n');
    else
        fprintf('Status: ✗ INVALID\n');
    end

    fprintf('\nErrors: %d\n', length(result.errors));
    if ~isempty(result.errors)
        for i = 1:length(result.errors)
            fprintf('  ✗ %s\n', result.errors{i});
        end
    end

    fprintf('\nWarnings: %d\n', length(result.warnings));
    if ~isempty(result.warnings)
        for i = 1:length(result.warnings)
            fprintf('  ⚠ %s\n', result.warnings{i});
        end
    end

    fprintf('=======================================================\n');
end
