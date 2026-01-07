function results = load_experiment_results(resultsFile)
% LOAD_EXPERIMENT_RESULTS - Load saved experiment results
%
% This function loads experiment results that were saved using
% save_experiment_results or run_experiment.
%
% Usage:
%    results = load_experiment_results('results/exp001_20260107_143022.mat')
%    results = load_experiment_results('exp001_20260107_143022')
%
% Inputs:
%    resultsFile - Full path to results file or base filename
%
% Outputs:
%    results     - Complete results structure containing:
%                  .bestSolution: Best parameters and cost
%                  .bestCosts: Convergence history
%                  .population: Final population
%                  .metadata: Experiment metadata
%                  .config: Configuration used
%                  .timing: Timing information
%                  .randomState: RNG state (if saved)
%
% Example:
%    % Load and analyze results
%    results = load_experiment_results('results/exp001_20260107_143022.mat');
%    fprintf('Final cost: %.6f\n', results.bestSolution.Cost);
%    figure; semilogy(results.bestCosts);
%
%    % Reproduce the random state
%    if isfield(results, 'randomState')
%        rng(results.randomState.initial);
%        % Run optimization again for exact reproduction
%    end
%
% See also: run_experiment, save_experiment_results
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Validate Input
    if nargin < 1
        error('LoadExperimentResults:NoInput', 'Results file path is required');
    end

    %% Build Full Path if Needed
    if ~endsWith(resultsFile, '.mat')
        resultsFile = [resultsFile '.mat'];
    end

    % If relative path, look in results directory
    if ~exist(resultsFile, 'file')
        [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
        projectRoot = fileparts(fileparts(scriptPath));
        alternativePath = fullfile(projectRoot, 'results', resultsFile);

        if exist(alternativePath, 'file')
            resultsFile = alternativePath;
        else
            error('LoadExperimentResults:FileNotFound', ...
                  'Results file not found: %s', resultsFile);
        end
    end

    %% Load Results
    fprintf('Loading experiment results from:\n  %s\n\n', resultsFile);

    try
        loadedData = load(resultsFile);
    catch ME
        error('LoadExperimentResults:LoadError', ...
              'Failed to load results file: %s', ME.message);
    end

    %% Extract Results
    results = struct();

    % Required fields
    requiredFields = {'bestSolution', 'bestCosts', 'metadata', 'config'};
    for i = 1:length(requiredFields)
        field = requiredFields{i};
        if isfield(loadedData, field)
            results.(field) = loadedData.(field);
        else
            warning('LoadExperimentResults:MissingField', ...
                    'Required field "%s" not found in results file', field);
        end
    end

    % Optional fields
    optionalFields = {'population', 'timing', 'randomState'};
    for i = 1:length(optionalFields)
        field = optionalFields{i};
        if isfield(loadedData, field)
            results.(field) = loadedData.(field);
        end
    end

    %% Display Summary
    if isfield(results, 'metadata')
        fprintf('Experiment Information:\n');
        fprintf('  Name: %s\n', results.metadata.experimentName);

        if isfield(results.metadata, 'description') && ~isempty(results.metadata.description)
            fprintf('  Description: %s\n', results.metadata.description);
        end

        fprintf('  Date: %s\n', results.metadata.startTime);
        fprintf('  Duration: %.2f minutes\n', results.metadata.executionTime / 60);
        fprintf('\n');

        fprintf('Results Summary:\n');
        fprintf('  Final Cost: %.6f\n', results.metadata.finalCost);
        fprintf('  Initial Cost: %.6f\n', results.metadata.initialCost);
        fprintf('  Improvement: %.2f%%\n', results.metadata.improvement);
        fprintf('\n');

        fprintf('Configuration:\n');
        if isfield(results, 'config')
            fprintf('  Iterations: %d\n', results.config.maxIterations);
            fprintf('  Population: %d-%d\n', results.config.initialPopSize, results.config.maxPopSize);
            fprintf('  Parallel: %s\n', mat2str(results.config.useParallel));
        end
        fprintf('\n');
    end

    fprintf('✓ Results loaded successfully\n\n');
end
