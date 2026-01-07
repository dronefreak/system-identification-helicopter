function data = load_flight_data(datasetName, varargin)
% LOAD_FLIGHT_DATA - Load and validate flight test data
%
% This function provides a standardized interface for loading flight data
% with automatic validation and error checking.
%
% Usage:
%    data = load_flight_data('best')
%    data = load_flight_data('best2')
%    data = load_flight_data('best', 'validate', true)
%    data = load_flight_data('/full/path/to/file.mat')
%
% Inputs:
%    datasetName - Name of dataset ('best', 'best2') or full file path
%
% Optional Parameters (Name-Value pairs):
%    'validate'     - Run validation checks (default: true)
%    'verbose'      - Display loading information (default: true)
%    'dataDir'      - Base data directory (default: 'data/experiments/')
%
% Outputs:
%    data           - Structure containing flight data:
%                     .in_H: Input data (4 × N)
%                     .out_H: Output data (10 × N)
%                     .time: Time vector (N × 1)
%                     .inr: Reference input matrix
%                     .outr: Reference output matrix
%                     .metadata: File and validation metadata
%
% Examples:
%    % Load primary dataset
%    data = load_flight_data('best');
%
%    % Load without validation (faster)
%    data = load_flight_data('best2', 'validate', false);
%
%    % Load from custom path
%    data = load_flight_data('/path/to/custom_data.mat');
%
% See also: validate_flight_data, save_flight_data
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Input Arguments
    p = inputParser;
    addRequired(p, 'datasetName', @ischar);
    addParameter(p, 'validate', true, @islogical);
    addParameter(p, 'verbose', true, @islogical);
    addParameter(p, 'dataDir', 'data/experiments/', @ischar);
    parse(p, datasetName, varargin{:});

    datasetName = p.Results.datasetName;
    doValidate = p.Results.validate;
    verbose = p.Results.verbose;
    dataDir = p.Results.dataDir;

    %% Determine File Path
    if exist(datasetName, 'file')
        % Full path provided
        filePath = datasetName;
    else
        % Dataset name provided - build path
        [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
        projectRoot = fileparts(fileparts(scriptPath));

        % Add .mat extension if needed
        if ~endsWith(datasetName, '.mat')
            datasetName = [datasetName '.mat'];
        end

        filePath = fullfile(projectRoot, dataDir, datasetName);
    end

    %% Check File Exists
    if ~exist(filePath, 'file')
        error('LoadFlightData:FileNotFound', ...
              'Data file not found: %s\nAvailable datasets: %s', ...
              filePath, list_available_datasets(dataDir));
    end

    %% Load Data
    if verbose
        fprintf('Loading flight data from: %s\n', filePath);
    end

    try
        loadedData = load(filePath);
    catch ME
        error('LoadFlightData:LoadError', ...
              'Failed to load data file: %s\nError: %s', filePath, ME.message);
    end

    %% Extract Required Variables
    data = struct();
    data.metadata = struct();
    data.metadata.filePath = filePath;
    data.metadata.fileName = datasetName;
    data.metadata.loadTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    % Required variables
    requiredVars = {'in_H', 'out_H', 'time', 'inr', 'outr'};
    missingVars = {};

    for i = 1:length(requiredVars)
        varName = requiredVars{i};
        if isfield(loadedData, varName)
            data.(varName) = loadedData.(varName);
        else
            missingVars{end+1} = varName; %#ok<AGROW>
        end
    end

    if ~isempty(missingVars)
        warning('LoadFlightData:MissingVariables', ...
                'Missing required variables: %s\nFile may be incomplete.', ...
                strjoin(missingVars, ', '));
    end

    %% Add Metadata
    if isfield(data, 'in_H')
        data.metadata.nSamples = size(data.in_H, 2);
        data.metadata.nInputs = size(data.in_H, 1);
    end

    if isfield(data, 'out_H')
        data.metadata.nOutputs = size(data.out_H, 1);
    end

    if isfield(data, 'time')
        if isvector(data.time)
            data.metadata.duration = max(data.time) - min(data.time);
            data.metadata.samplingRate = 1 / mean(diff(data.time));
        end
    end

    % File size
    fileInfo = dir(filePath);
    data.metadata.fileSize = fileInfo.bytes;
    data.metadata.fileSizeMB = fileInfo.bytes / 1024^2;

    %% Validation
    if doValidate
        if verbose
            fprintf('Validating data...\n');
        end

        validationResult = validate_flight_data_struct(data);

        if ~validationResult.isValid
            warning('LoadFlightData:ValidationFailed', ...
                    'Data validation failed:\n%s', ...
                    strjoin(validationResult.errors, '\n'));
        elseif verbose
            fprintf('✓ Data validation passed\n');
        end

        data.metadata.validation = validationResult;
    end

    %% Display Summary
    if verbose
        fprintf('\n');
        fprintf('Flight Data Summary:\n');
        fprintf('-------------------\n');
        fprintf('File: %s\n', datasetName);
        fprintf('Size: %.2f MB\n', data.metadata.fileSizeMB);
        if isfield(data.metadata, 'nSamples')
            fprintf('Samples: %d\n', data.metadata.nSamples);
        end
        if isfield(data.metadata, 'duration')
            fprintf('Duration: %.2f seconds\n', data.metadata.duration);
            fprintf('Sampling Rate: %.1f Hz\n', data.metadata.samplingRate);
        end
        fprintf('\n');
    end
end

function validation = validate_flight_data_struct(data)
    % Internal validation function
    validation = struct();
    validation.isValid = true;
    validation.errors = {};
    validation.warnings = {};

    % Check required fields
    requiredFields = {'in_H', 'out_H', 'time'};
    for i = 1:length(requiredFields)
        if ~isfield(data, requiredFields{i})
            validation.isValid = false;
            validation.errors{end+1} = sprintf('Missing field: %s', requiredFields{i});
        end
    end

    if ~validation.isValid
        return;
    end

    % Check dimensions
    if size(data.in_H, 1) ~= 4
        validation.errors{end+1} = 'in_H should have 4 rows (control inputs)';
        validation.isValid = false;
    end

    if size(data.out_H, 1) ~= 10
        validation.errors{end+1} = 'out_H should have 10 rows (measurements)';
        validation.isValid = false;
    end

    if size(data.in_H, 2) ~= size(data.out_H, 2)
        validation.errors{end+1} = 'in_H and out_H must have same number of samples';
        validation.isValid = false;
    end

    % Check for NaN or Inf
    if any(isnan(data.in_H(:))) || any(isinf(data.in_H(:)))
        validation.warnings{end+1} = 'in_H contains NaN or Inf values';
    end

    if any(isnan(data.out_H(:))) || any(isinf(data.out_H(:)))
        validation.warnings{end+1} = 'out_H contains NaN or Inf values';
    end

    % Check time vector
    if ~isvector(data.time)
        validation.errors{end+1} = 'time must be a vector';
        validation.isValid = false;
    elseif ~issorted(data.time)
        validation.errors{end+1} = 'time vector is not monotonically increasing';
        validation.isValid = false;
    end
end

function datasetList = list_available_datasets(dataDir)
    % List available datasets in the data directory
    [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(scriptPath));
    fullDataDir = fullfile(projectRoot, dataDir);

    if ~exist(fullDataDir, 'dir')
        datasetList = '(Data directory not found)';
        return;
    end

    files = dir(fullfile(fullDataDir, '*.mat'));

    if isempty(files)
        datasetList = '(No .mat files found)';
    else
        datasetList = sprintf('\nAvailable datasets:\n');
        for i = 1:min(10, length(files))
            [~, name, ~] = fileparts(files(i).name);
            datasetList = sprintf('%s  - %s\n', datasetList, name);
        end
        if length(files) > 10
            datasetList = sprintf('%s  ... and %d more\n', datasetList, length(files) - 10);
        end
    end
end
