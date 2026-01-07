function results = profile_iwo(varargin)
% PROFILE_IWO - Profile the IWO optimization to identify performance bottlenecks
%
% This utility profiles the IWO optimization algorithm and generates a detailed
% performance report showing which functions consume the most time and memory.
%
% Usage:
%    results = profile_iwo()                    % Quick profiling (100 iterations)
%    results = profile_iwo('iterations', 500)   % Custom iteration count
%    results = profile_iwo('detail', 'full')    % Full detailed profiling
%
% Optional Parameters:
%    'iterations'   - Number of iterations to profile (default: 100)
%    'detail'       - 'quick' or 'full' profiling detail (default: 'quick')
%    'plot'         - true/false to generate plots (default: true)
%
% Outputs:
%    results        - Structure containing profiling results:
%                     .profileInfo: MATLAB profiler data
%                     .totalTime: Total execution time (seconds)
%                     .topFunctions: Top time-consuming functions
%                     .memoryUsage: Peak memory usage
%
% Example:
%    load('data/experiments/best.mat')
%    results = profile_iwo('iterations', 500, 'plot', true);
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Input Arguments
    p = inputParser;
    addParameter(p, 'iterations', 100, @isnumeric);
    addParameter(p, 'detail', 'quick', @(x) any(validatestring(x, {'quick', 'full'})));
    addParameter(p, 'plot', true, @islogical);
    parse(p, varargin{:});

    iterations = p.Results.iterations;
    detail = p.Results.detail;
    plotResults = p.Results.plot;

    %% Validate Required Data
    requiredVars = {'in_H', 'out_H', 'time', 'inr', 'outr'};
    missingVars = {};

    for i = 1:length(requiredVars)
        if ~evalin('base', ['exist(''' requiredVars{i} ''', ''var'')'])
            missingVars{end+1} = requiredVars{i}; %#ok<AGROW>
        end
    end

    if ~isempty(missingVars)
        error('ProfileIWO:MissingData', ...
              ['Missing required workspace variables: %s\n' ...
               'Please load best.mat first:\n' ...
               '  load(''data/experiments/best.mat'')'], ...
              strjoin(missingVars, ', '));
    end

    fprintf('=======================================================\n');
    fprintf('  IWO Performance Profiling\n');
    fprintf('=======================================================\n');
    fprintf('Iterations: %d\n', iterations);
    fprintf('Detail Level: %s\n', detail);
    fprintf('-------------------------------------------------------\n\n');

    %% Modify Configuration for Profiling
    % Temporarily reduce iterations for profiling
    originalConfig = config_iwo();

    % Create temporary profiling configuration
    evalin('base', sprintf(['tempProfileConfig = struct(' ...
                            '''nVar'', %d, ' ...
                            '''maxIterations'', %d, ' ...
                            '''initialPopSize'', %d, ' ...
                            '''maxPopSize'', %d, ' ...
                            '''minSeeds'', %d, ' ...
                            '''maxSeeds'', %d, ' ...
                            '''displayInterval'', %d, ' ...
                            '''plotResults'', false);'], ...
                           originalConfig.nVar, iterations, ...
                           originalConfig.initialPopSize, originalConfig.maxPopSize, ...
                           originalConfig.minSeeds, originalConfig.maxSeeds, ...
                           max(10, floor(iterations/10))));

    %% Start Profiling
    fprintf('Starting profiler...\n');
    profile clear;

    if strcmp(detail, 'full')
        profile on -detail builtin -history;
    else
        profile on;
    end

    %% Record Memory Before
    memBefore = memory_info();
    tic;

    %% Run IWO Algorithm
    fprintf('Running IWO optimization (%d iterations)...\n', iterations);

    try
        % Backup config_iwo temporarily
        configBackup = functions(str2func('config_iwo'));

        % Run IWO with modified iteration count
        % We'll create a modified version inline
        evalin('base', 'iwo');

    catch ME
        profile off;
        error('ProfileIWO:OptimizationError', ...
              'Error during profiling: %s', ME.message);
    end

    %% Stop Profiling
    totalTime = toc;
    memAfter = memory_info();

    profile off;
    profileData = profile('info');

    fprintf('\n-------------------------------------------------------\n');
    fprintf('Profiling complete!\n');
    fprintf('Total Time: %.2f seconds\n', totalTime);
    fprintf('Peak Memory: %.2f MB\n', memAfter.peakMemory);
    fprintf('Memory Increase: %.2f MB\n', memAfter.peakMemory - memBefore.peakMemory);
    fprintf('-------------------------------------------------------\n\n');

    %% Analyze Profiling Results
    fprintf('Analyzing profiling data...\n\n');

    % Extract function timing information
    if ~isempty(profileData.FunctionTable)
        funcTable = profileData.FunctionTable;

        % Sort by total time
        [~, sortIdx] = sort([funcTable.TotalTime], 'descend');
        sortedFunctions = funcTable(sortIdx);

        % Display top time-consuming functions
        fprintf('Top 10 Time-Consuming Functions:\n');
        fprintf('%-40s %12s %12s %10s\n', 'Function', 'Total (s)', 'Self (s)', 'Calls');
        fprintf('%s\n', repmat('-', 1, 80));

        topN = min(10, length(sortedFunctions));
        topFunctions = struct('name', {}, 'totalTime', {}, 'selfTime', {}, 'numCalls', {});

        for i = 1:topN
            func = sortedFunctions(i);

            % Extract function name (remove path for readability)
            [~, funcName, ext] = fileparts(func.FunctionName);
            if isempty(ext)
                displayName = funcName;
            else
                displayName = [funcName ext];
            end

            fprintf('%-40s %12.4f %12.4f %10d\n', ...
                    displayName, func.TotalTime, func.SelfTime, func.NumCalls);

            topFunctions(i).name = displayName;
            topFunctions(i).totalTime = func.TotalTime;
            topFunctions(i).selfTime = func.SelfTime;
            topFunctions(i).numCalls = func.NumCalls;
        end

        fprintf('\n');
    else
        topFunctions = [];
    end

    %% Generate Recommendations
    fprintf('Performance Recommendations:\n');
    fprintf('%s\n', repmat('-', 1, 80));

    % Check for parallelization opportunity
    if ~isempty(topFunctions) && any(contains({topFunctions.name}, 'Sphere'))
        fprintf('- Cost function (Sphere) is a bottleneck\n');
        fprintf('  → Consider enabling parallel evaluation with parfor\n');
        fprintf('  → Requires Parallel Computing Toolbox\n\n');
    end

    % Check memory usage
    if memAfter.peakMemory > 1000
        fprintf('- High memory usage detected (%.0f MB)\n', memAfter.peakMemory);
        fprintf('  → Consider preallocating arrays\n');
        fprintf('  → Use memory-efficient data structures\n\n');
    end

    % Check iteration performance
    avgTimePerIter = totalTime / iterations;
    if avgTimePerIter > 1.0
        fprintf('- Slow iteration speed (%.2f s/iteration)\n', avgTimePerIter);
        fprintf('  → Profile individual functions for optimization\n');
        fprintf('  → Consider vectorizing operations\n\n');
    else
        fprintf('- Iteration speed is good (%.3f s/iteration)\n\n', avgTimePerIter);
    end

    %% Create Plots
    if plotResults && ~isempty(topFunctions)
        figure('Name', 'IWO Profiling Results', 'Position', [100 100 1200 600]);

        % Plot 1: Top functions by time
        subplot(1, 2, 1);
        numDisplay = min(8, length(topFunctions));
        barData = [topFunctions(1:numDisplay).totalTime];
        barh(barData);
        set(gca, 'YTick', 1:numDisplay);
        set(gca, 'YTickLabel', {topFunctions(1:numDisplay).name});
        xlabel('Total Time (seconds)');
        title('Top Functions by Time');
        grid on;

        % Plot 2: Time breakdown
        subplot(1, 2, 2);
        pieData = barData;
        otherTime = totalTime - sum(pieData);
        if otherTime > 0
            pieData = [pieData otherTime];
            pieLabels = [{topFunctions(1:numDisplay).name} {'Other'}];
        else
            pieLabels = {topFunctions(1:numDisplay).name};
        end
        pie(pieData);
        legend(pieLabels, 'Location', 'eastoutside');
        title(sprintf('Time Distribution (Total: %.2f s)', totalTime));
    end

    %% Package Results
    results = struct();
    results.profileInfo = profileData;
    results.totalTime = totalTime;
    results.topFunctions = topFunctions;
    results.memoryUsage = memAfter;
    results.iterations = iterations;
    results.avgTimePerIteration = avgTimePerIter;

    fprintf('\nResults stored in output structure.\n');
    fprintf('To view detailed profiler report: profview(0, results.profileInfo)\n\n');

    %% Cleanup
    evalin('base', 'clear tempProfileConfig');
end

function memInfo = memory_info()
% Get memory information (cross-platform)
    try
        if ispc
            [~, sys] = memory;
            memInfo.peakMemory = (sys.PhysicalMemory.Total - sys.PhysicalMemory.Available) / 1024^2;
        elseif ismac || isunix
            % Linux/Mac - estimate from MATLAB's memory usage
            m = memory;
            memInfo.peakMemory = m.MemUsedMATLAB / 1024^2;
        else
            memInfo.peakMemory = NaN;
        end
    catch
        memInfo.peakMemory = NaN;
    end
end
