function outputFile = save_experiment_results(results, experimentName, outputDir)
% SAVE_EXPERIMENT_RESULTS - Save experiment results with complete metadata
%
% This function saves experiment results in a structured format with full
% metadata for reproducibility and later analysis.
%
% Usage:
%    save_experiment_results(results, 'exp001', 'results/')
%    outputFile = save_experiment_results(results, 'test_run')
%
% Inputs:
%    results        - Results structure from run_experiment
%    experimentName - Name for the experiment
%    outputDir      - Output directory (default: 'results/')
%
% Outputs:
%    outputFile     - Full path to saved results file
%
% File Structure:
%    The saved .mat file contains:
%    - bestSolution: Best parameters and cost
%    - bestCosts: Convergence history
%    - population: Final population
%    - metadata: Complete experiment metadata
%    - config: Configuration used
%    - timing: Timing information
%    - randomState: Random number generator state (if saved)
%
% Additionally, a JSON metadata file is created for easy inspection.
%
% See also: run_experiment, load_experiment_results
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Validate Inputs
    if nargin < 2
        error('SaveExperimentResults:NotEnoughInputs', ...
              'Results and experiment name are required');
    end

    if nargin < 3
        outputDir = 'results/';
    end

    %% Prepare Output Directory
    % Get project root if relative path
    if ~isempty(outputDir) && outputDir(1) ~= '/' && ~(length(outputDir) > 1 && outputDir(2) == ':')
        [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
        projectRoot = fileparts(fileparts(scriptPath));
        outputDir = fullfile(projectRoot, outputDir);
    end

    % Create output directory if needed
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
        fprintf('Created output directory: %s\n', outputDir);
    end

    %% Generate Output Filename with Timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    baseFilename = sprintf('%s_%s', experimentName, timestamp);
    outputFile = fullfile(outputDir, [baseFilename '.mat']);

    %% Extract Results Components
    if isfield(results, 'bestSolution')
        bestSolution = results.bestSolution;
    else
        error('SaveExperimentResults:MissingField', 'Results must contain bestSolution');
    end

    if isfield(results, 'bestCosts')
        bestCosts = results.bestCosts;
    else
        error('SaveExperimentResults:MissingField', 'Results must contain bestCosts');
    end

    if isfield(results, 'population')
        population = results.population;
    else
        warning('SaveExperimentResults:MissingPopulation', ...
                'Population not found in results');
        population = [];
    end

    if isfield(results, 'metadata')
        metadata = results.metadata;
    else
        warning('SaveExperimentResults:MissingMetadata', ...
                'Metadata not found in results');
        metadata = struct();
    end

    if isfield(results, 'config')
        config = results.config;
    else
        warning('SaveExperimentResults:MissingConfig', ...
                'Configuration not found in results');
        config = struct();
    end

    if isfield(results, 'timing')
        timing = results.timing;
    else
        timing = struct();
    end

    if isfield(results, 'randomState')
        randomState = results.randomState;
    end

    %% Add File Information to Metadata
    metadata.savedFile = outputFile;
    metadata.savedTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    %% Save MATLAB Results File
    fprintf('Saving results to: %s\n', outputFile);

    try
        if exist('randomState', 'var')
            save(outputFile, 'bestSolution', 'bestCosts', 'population', ...
                 'metadata', 'config', 'timing', 'randomState', '-v7.3');
        else
            save(outputFile, 'bestSolution', 'bestCosts', 'population', ...
                 'metadata', 'config', 'timing', '-v7.3');
        end
        fprintf('✓ Results saved successfully\n');
    catch ME
        error('SaveExperimentResults:SaveError', ...
              'Failed to save results: %s', ME.message);
    end

    %% Create JSON Metadata File
    jsonFile = fullfile(outputDir, [baseFilename '_metadata.json']);

    try
        % Create simplified metadata for JSON (MATLAB structures don't convert well)
        jsonMetadata = struct();
        jsonMetadata.experimentName = metadata.experimentName;
        jsonMetadata.startTime = metadata.startTime;
        jsonMetadata.endTime = metadata.endTime;
        jsonMetadata.executionTime = metadata.executionTime;
        jsonMetadata.finalCost = metadata.finalCost;
        jsonMetadata.initialCost = metadata.initialCost;
        jsonMetadata.improvement = metadata.improvement;
        jsonMetadata.platform = metadata.platform;
        jsonMetadata.matlabRelease = metadata.matlabRelease;
        jsonMetadata.numCores = metadata.numCores;

        if isfield(metadata, 'description')
            jsonMetadata.description = metadata.description;
        end

        if isfield(metadata, 'tags')
            jsonMetadata.tags = metadata.tags;
        end

        % Configuration summary
        jsonMetadata.config = struct();
        jsonMetadata.config.maxIterations = config.maxIterations;
        jsonMetadata.config.initialPopSize = config.initialPopSize;
        jsonMetadata.config.maxPopSize = config.maxPopSize;
        jsonMetadata.config.useParallel = config.useParallel;

        % Write JSON (MATLAB R2016b+)
        if exist('jsonencode', 'builtin')
            jsonText = jsonencode(jsonMetadata);
            % Pretty print
            jsonText = strrep(jsonText, ',', sprintf(',\n  '));
            jsonText = strrep(jsonText, '{', sprintf('{\n  '));
            jsonText = strrep(jsonText, '}', sprintf('\n}'));

            fid = fopen(jsonFile, 'w');
            if fid > 0
                fprintf(fid, '%s', jsonText);
                fclose(fid);
                fprintf('✓ Metadata saved to: %s\n', jsonFile);
            end
        end
    catch ME
        warning('SaveExperimentResults:JSONError', ...
                'Could not save JSON metadata: %s', ME.message);
    end

    %% Create Summary Text File
    summaryFile = fullfile(outputDir, [baseFilename '_summary.txt']);

    try
        fid = fopen(summaryFile, 'w');
        if fid > 0
            fprintf(fid, '=======================================================\n');
            fprintf(fid, '  IWO Experiment Results Summary\n');
            fprintf(fid, '=======================================================\n\n');

            fprintf(fid, 'Experiment: %s\n', metadata.experimentName);
            if isfield(metadata, 'description') && ~isempty(metadata.description)
                fprintf(fid, 'Description: %s\n', metadata.description);
            end
            fprintf(fid, '\n');

            fprintf(fid, 'Execution Details:\n');
            fprintf(fid, '  Start Time: %s\n', metadata.startTime);
            fprintf(fid, '  End Time: %s\n', metadata.endTime);
            fprintf(fid, '  Duration: %.2f minutes\n', metadata.executionTime / 60);
            fprintf(fid, '\n');

            fprintf(fid, 'Results:\n');
            fprintf(fid, '  Initial Cost: %.6f\n', metadata.initialCost);
            fprintf(fid, '  Final Cost: %.6f\n', metadata.finalCost);
            fprintf(fid, '  Improvement: %.2f%%\n', metadata.improvement);
            fprintf(fid, '\n');

            fprintf(fid, 'Configuration:\n');
            fprintf(fid, '  Iterations: %d\n', config.maxIterations);
            fprintf(fid, '  Population: %d-%d\n', config.initialPopSize, config.maxPopSize);
            fprintf(fid, '  Parallel: %s\n', mat2str(config.useParallel));
            fprintf(fid, '\n');

            fprintf(fid, 'System Information:\n');
            fprintf(fid, '  Platform: %s\n', metadata.platform);
            fprintf(fid, '  MATLAB: %s\n', metadata.matlabRelease);
            fprintf(fid, '  Cores: %d\n', metadata.numCores);
            fprintf(fid, '\n');

            fprintf(fid, 'Best Parameters (40 values):\n');
            for i = 1:length(bestSolution.Position)
                fprintf(fid, '  %2d: %12.6f\n', i, bestSolution.Position(i));
            end

            fprintf(fid, '\n');
            fprintf(fid, 'Files:\n');
            fprintf(fid, '  Results: %s\n', outputFile);
            fprintf(fid, '  Metadata: %s\n', jsonFile);
            fprintf(fid, '  Summary: %s\n', summaryFile);
            fprintf(fid, '\n');

            fprintf(fid, '=======================================================\n');

            fclose(fid);
            fprintf('✓ Summary saved to: %s\n', summaryFile);
        end
    catch ME
        warning('SaveExperimentResults:SummaryError', ...
                'Could not save summary file: %s', ME.message);
    end

    fprintf('\n');
end
