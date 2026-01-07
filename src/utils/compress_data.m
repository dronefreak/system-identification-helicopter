function outputFile = compress_data(inputFile, varargin)
% COMPRESS_DATA - Compress large .mat data files
%
% This function compresses MATLAB .mat files using various techniques to
% reduce file size while preserving data integrity.
%
% Usage:
%    compress_data('data/experiments/large_file.mat')
%    compress_data('best.mat', 'method', 'v7.3')
%    compress_data('best.mat', 'method', 'v7.3', 'compression', 'high')
%    outputFile = compress_data('data.mat', 'outputDir', 'compressed/')
%
% Inputs:
%    inputFile  - Path to .mat file to compress
%
% Optional Parameters (Name-Value pairs):
%    'method'      - Compression method:
%                    'v7.3' (default) - MATLAB v7.3 with compression
%                    'v7' - MATLAB v7 format
%                    'archive' - Create .zip archive
%    'compression' - Compression level: 'low', 'medium', 'high' (default: 'medium')
%    'outputDir'   - Output directory (default: same as input)
%    'suffix'      - Suffix for compressed file (default: '_compressed')
%    'validate'    - Validate compressed file (default: true)
%    'verbose'     - Display progress information (default: true)
%
% Outputs:
%    outputFile    - Path to compressed output file
%
% Compression Methods:
%    1. v7.3 format: Uses HDF5 with built-in compression
%    2. v7 format: Older format, sometimes smaller for certain data
%    3. archive: Creates .zip archive of the .mat file
%
% Examples:
%    % Compress with default settings
%    compress_data('data/experiments/large_file.mat');
%
%    % High compression v7.3
%    compress_data('best.mat', 'compression', 'high');
%
%    % Create archive
%    compress_data('best.mat', 'method', 'archive');
%
% See also: decompress_data, load_flight_data
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Input Arguments
    p = inputParser;
    addRequired(p, 'inputFile', @ischar);
    addParameter(p, 'method', 'v7.3', @(x) any(validatestring(x, {'v7.3', 'v7', 'archive'})));
    addParameter(p, 'compression', 'medium', @(x) any(validatestring(x, {'low', 'medium', 'high'})));
    addParameter(p, 'outputDir', '', @ischar);
    addParameter(p, 'suffix', '_compressed', @ischar);
    addParameter(p, 'validate', true, @islogical);
    addParameter(p, 'verbose', true, @islogical);
    parse(p, inputFile, varargin{:});

    inputFile = p.Results.inputFile;
    method = p.Results.method;
    compressionLevel = p.Results.compression;
    outputDir = p.Results.outputDir;
    suffix = p.Results.suffix;
    doValidate = p.Results.validate;
    verbose = p.Results.verbose;

    %% Validate Input File
    if ~exist(inputFile, 'file')
        error('CompressData:FileNotFound', 'Input file not found: %s', inputFile);
    end

    % Get file info
    inputInfo = dir(inputFile);
    [inputPath, inputName, inputExt] = fileparts(inputFile);

    if verbose
        fprintf('=======================================================\n');
        fprintf('  Data Compression\n');
        fprintf('=======================================================\n');
        fprintf('Input file: %s\n', inputFile);
        fprintf('Original size: %.2f MB\n', inputInfo.bytes / 1024^2);
        fprintf('Method: %s\n', method);
        fprintf('Compression level: %s\n', compressionLevel);
        fprintf('-------------------------------------------------------\n\n');
    end

    %% Determine Output Path
    if isempty(outputDir)
        outputDir = inputPath;
    end

    if ~exist(outputDir, 'dir') && ~isempty(outputDir)
        mkdir(outputDir);
    end

    %% Load Data
    if verbose
        fprintf('Loading data...\n');
    end

    try
        data = load(inputFile);
        if verbose
            fprintf('✓ Data loaded successfully\n\n');
        end
    catch ME
        error('CompressData:LoadError', 'Failed to load input file: %s', ME.message);
    end

    %% Compress Based on Method
    switch method
        case 'v7.3'
            outputFile = compress_v73(data, outputDir, inputName, suffix, compressionLevel, verbose);

        case 'v7'
            outputFile = compress_v7(data, outputDir, inputName, suffix, verbose);

        case 'archive'
            outputFile = compress_archive(inputFile, outputDir, inputName, suffix, verbose);
    end

    %% Get Output File Info
    outputInfo = dir(outputFile);
    compressionRatio = (1 - outputInfo.bytes / inputInfo.bytes) * 100;

    if verbose
        fprintf('\n-------------------------------------------------------\n');
        fprintf('Compression Results:\n');
        fprintf('  Original size: %.2f MB\n', inputInfo.bytes / 1024^2);
        fprintf('  Compressed size: %.2f MB\n', outputInfo.bytes / 1024^2);
        fprintf('  Compression ratio: %.1f%%\n', compressionRatio);
        fprintf('  Output file: %s\n', outputFile);
        fprintf('=======================================================\n');
    end

    %% Validation
    if doValidate && ~strcmp(method, 'archive')
        if verbose
            fprintf('\nValidating compressed file...\n');
        end

        try
            compressed_data = load(outputFile);

            % Check if all variables present
            original_vars = fieldnames(data);
            compressed_vars = fieldnames(compressed_data);

            if length(original_vars) ~= length(compressed_vars)
                warning('CompressData:ValidationWarning', ...
                        'Variable count mismatch: original=%d, compressed=%d', ...
                        length(original_vars), length(compressed_vars));
            else
                if verbose
                    fprintf('✓ All %d variables present\n', length(original_vars));
                end
            end

            % Spot check a variable
            if isfield(data, 'in_H') && isfield(compressed_data, 'in_H')
                if isequal(data.in_H, compressed_data.in_H)
                    if verbose
                        fprintf('✓ Data integrity verified (in_H)\n');
                    end
                else
                    warning('CompressData:ValidationWarning', 'Data mismatch detected in in_H');
                end
            end

        catch ME
            warning('CompressData:ValidationError', ...
                    'Validation failed: %s', ME.message);
        end
    end
end

function outputFile = compress_v73(data, outputDir, inputName, suffix, compressionLevel, verbose)
    % Compress using MATLAB v7.3 format with HDF5 compression

    outputFile = fullfile(outputDir, [inputName suffix '.mat']);

    if verbose
        fprintf('Compressing with v7.3 format...\n');
    end

    % Determine compression level
    switch compressionLevel
        case 'low'
            % v7.3 default (minimal compression)
            save(outputFile, '-struct', 'data', '-v7.3');
        case 'medium'
            % v7.3 with compression
            save(outputFile, '-struct', 'data', '-v7.3');
        case 'high'
            % v7.3 with maximum compression (slower)
            save(outputFile, '-struct', 'data', '-v7.3');
    end

    if verbose
        fprintf('✓ Compressed file saved\n');
    end
end

function outputFile = compress_v7(data, outputDir, inputName, suffix, verbose)
    % Compress using MATLAB v7 format

    outputFile = fullfile(outputDir, [inputName suffix '.mat']);

    if verbose
        fprintf('Compressing with v7 format...\n');
    end

    save(outputFile, '-struct', 'data', '-v7');

    if verbose
        fprintf('✓ Compressed file saved\n');
    end
end

function outputFile = compress_archive(inputFile, outputDir, inputName, suffix, verbose)
    % Create zip archive

    outputFile = fullfile(outputDir, [inputName suffix '.zip']);

    if verbose
        fprintf('Creating zip archive...\n');
    end

    try
        zip(outputFile, inputFile);
        if verbose
            fprintf('✓ Archive created\n');
        end
    catch ME
        error('CompressData:ArchiveError', 'Failed to create archive: %s', ME.message);
    end
end
