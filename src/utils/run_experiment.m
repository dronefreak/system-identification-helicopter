function results = run_experiment(experimentName, configName, varargin)
% RUN_EXPERIMENT - Run a complete IWO experiment with full tracking
%
% This function provides a complete experiment workflow with automatic tracking,
% metadata collection, and results management. It's the recommended way to run
% reproducible experiments.
%
% Usage:
%    results = run_experiment('exp001', 'default')
%    results = run_experiment('test_run', 'fast_test', 'description', 'Quick test')
%    results = run_experiment('parallel_test', 'parallel_optimized', ...
%                             'dataFile', 'best2.mat', 'saveResults', true)
%
% Inputs:
%    experimentName - Unique name for this experiment (used for output files)
%    configName     - Configuration name to load (or 'default' for config_iwo)
%
% Optional Parameters (Name-Value pairs):
%    'description'  - Text description of the experiment
%    'dataFile'     - Data file to load (default: 'data/experiments/best.mat')
%    'saveResults'  - Save results to file (default: true)
%    'outputDir'    - Directory for results (default: 'results/')
%    'tags'         - Cell array of tags for organization
%
% Outputs:
%    results        - Complete results structure with:
%                     .bestSolution: Best parameters found
%                     .bestCosts: Convergence history
%                     .population: Final population
%                     .metadata: Experiment metadata
%                     .config: Configuration used
%                     .timing: Execution time information
%
% Example:
%    % Run a reproducible experiment
%    results = run_experiment('exp_001', 'default', ...
%                             'description', 'Baseline IWO run', ...
%                             'tags', {'baseline', 'reproducible'});
%
%    % Quick test run
%    results = run_experiment('quick_test', 'fast_test', ...
%                             'saveResults', false);
%
% See also: save_experiment_results, load_config, config_iwo
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Input Arguments
    p = inputParser;
    addRequired(p, 'experimentName', @ischar);
    addRequired(p, 'configName', @ischar);
    addParameter(p, 'description', '', @ischar);
    addParameter(p, 'dataFile', 'data/experiments/best.mat', @ischar);
    addParameter(p, 'saveResults', true, @islogical);
    addParameter(p, 'outputDir', 'results/', @ischar);
    addParameter(p, 'tags', {}, @iscell);
    parse(p, experimentName, configName, varargin{:});

    experimentName = p.Results.experimentName;
    configName = p.Results.configName;
    description = p.Results.description;
    dataFile = p.Results.dataFile;
    saveResults = p.Results.saveResults;
    outputDir = p.Results.outputDir;
    tags = p.Results.tags;

    %% Display Experiment Header
    fprintf('=======================================================\n');
    fprintf('  IWO Experiment: %s\n', experimentName);
    fprintf('=======================================================\n');
    fprintf('Configuration: %s\n', configName);
    if ~isempty(description)
        fprintf('Description: %s\n', description);
    end
    fprintf('Data File: %s\n', dataFile);
    fprintf('-------------------------------------------------------\n\n');

    %% Load Configuration
    if strcmp(configName, 'default')
        config = config_iwo();
        fprintf('Using default configuration\n\n');
    else
        try
            config = load_config(configName);
            fprintf('\n');
        catch ME
            error('RunExperiment:ConfigLoadError', ...
                  'Failed to load configuration "%s": %s', configName, ME.message);
        end
    end

    %% Load Data
    fprintf('Loading flight data...\n');

    % Get project root
    [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(scriptPath));

    % Build full data path if relative
    if ~isempty(dataFile) && dataFile(1) ~= '/' && ~(length(dataFile) > 1 && dataFile(2) == ':')
        dataFile = fullfile(projectRoot, dataFile);
    end

    if ~exist(dataFile, 'file')
        error('RunExperiment:DataNotFound', ...
              'Data file not found: %s', dataFile);
    end

    try
        % Load data into base workspace (required by iwo.m)
        evalin('base', ['load(''' dataFile ''')']);
        fprintf('✓ Data loaded successfully\n\n');
    catch ME
        error('RunExperiment:DataLoadError', ...
              'Failed to load data file: %s', ME.message);
    end

    %% Collect Pre-Experiment Metadata
    metadata = struct();
    metadata.experimentName = experimentName;
    metadata.configName = configName;
    metadata.description = description;
    metadata.dataFile = dataFile;
    metadata.tags = tags;
    metadata.startTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    metadata.startTimestamp = now;

    % System information
    metadata.platform = computer;
    metadata.matlabVersion = version;
    metadata.matlabRelease = version('-release');

    % Toolbox availability
    metadata.hasParallelToolbox = license('test', 'Distrib_Computing_Toolbox');
    metadata.hasControlToolbox = license('test', 'Control_Toolbox');
    metadata.hasStatsToolbox = license('test', 'Statistics_Toolbox');

    % CPU information
    try
        metadata.numCores = feature('numcores');
    catch
        metadata.numCores = NaN;
    end

    fprintf('System Information:\n');
    fprintf('  Platform: %s\n', metadata.platform);
    fprintf('  MATLAB: %s\n', metadata.matlabRelease);
    fprintf('  CPU Cores: %d\n', metadata.numCores);
    fprintf('  Parallel Toolbox: %s\n', mat2str(metadata.hasParallelToolbox));
    fprintf('\n');

    %% Assign Configuration to Base Workspace
    % iwo.m expects config_iwo() to be called, so we override it temporarily
    assignin('base', 'EXPERIMENT_CONFIG_OVERRIDE', config);

    %% Run IWO Optimization
    fprintf('Starting optimization...\n');
    fprintf('-------------------------------------------------------\n\n');

    startTic = tic;
    try
        % Change to IWO directory and run
        currentDir = pwd;
        iwoDir = fullfile(projectRoot, 'src', 'algorithms', 'iwo', 'IWO');
        cd(iwoDir);

        % Run IWO (results stored in base workspace)
        iwo;

        % Return to original directory
        cd(currentDir);

        elapsedTime = toc(startTic);

    catch ME
        cd(currentDir);
        error('RunExperiment:OptimizationError', ...
              'Optimization failed: %s\n%s', ME.message, getReport(ME));
    end

    fprintf('\n-------------------------------------------------------\n');
    fprintf('Optimization completed in %.2f minutes\n', elapsedTime / 60);
    fprintf('-------------------------------------------------------\n\n');

    %% Collect Results
    results = struct();

    % Get results from base workspace
    try
        results.bestSolution = evalin('base', 'BestSol');
        results.bestCosts = evalin('base', 'BestCosts');
        results.population = evalin('base', 'pop');

        % Get random state if available
        if evalin('base', 'exist(''randomState'', ''var'')')
            results.randomState = evalin('base', 'randomState');
        end
    catch ME
        warning('RunExperiment:ResultsCollectionError', ...
                'Could not collect all results: %s', ME.message);
    end

    %% Complete Metadata
    metadata.endTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    metadata.endTimestamp = now;
    metadata.executionTime = elapsedTime;
    metadata.finalCost = results.bestSolution.Cost;
    metadata.initialCost = results.bestCosts(1);
    metadata.improvement = ((results.bestCosts(1) - results.bestSolution.Cost) / results.bestCosts(1)) * 100;

    % Configuration snapshot
    metadata.config = config;

    % Add to results
    results.metadata = metadata;
    results.config = config;
    results.timing = struct('total', elapsedTime, ...
                           'perIteration', elapsedTime / config.maxIterations);

    %% Display Summary
    fprintf('Experiment Results:\n');
    fprintf('  Final Cost: %.6f\n', metadata.finalCost);
    fprintf('  Initial Cost: %.6f\n', metadata.initialCost);
    fprintf('  Improvement: %.2f%%\n', metadata.improvement);
    fprintf('  Execution Time: %.2f minutes\n', elapsedTime / 60);
    fprintf('  Time/Iteration: %.3f seconds\n', results.timing.perIteration);
    fprintf('\n');

    %% Save Results
    if saveResults
        save_experiment_results(results, experimentName, outputDir);
    end

    %% Cleanup
    evalin('base', 'clear EXPERIMENT_CONFIG_OVERRIDE');

    fprintf('=======================================================\n');
    fprintf('Experiment "%s" completed successfully!\n', experimentName);
    fprintf('=======================================================\n');
end
