function config = configure_performance(varargin)
% CONFIGURE_PERFORMANCE - Configure performance options for IWO optimization
%
% This utility helps you configure performance settings for the IWO algorithm,
% including parallel computing, progress tracking, and memory optimization.
%
% Usage:
%    config = configure_performance()                     % Interactive mode
%    config = configure_performance('parallel', true)     % Enable parallel
%    config = configure_performance('parallel', true, ... % Multiple options
%                                    'progress', false)
%
% Optional Parameters:
%    'parallel'    - Enable/disable parallel computing (default: auto-detect)
%    'progress'    - Enable/disable progress bar (default: true)
%    'prealloc'    - Enable/disable memory preallocation (default: true)
%
% Outputs:
%    config        - Configuration structure with performance settings
%
% Examples:
%    % Enable parallel computing
%    config = configure_performance('parallel', true);
%
%    % Disable progress bar for batch processing
%    config = configure_performance('progress', false);
%
%    % Check if parallel computing is available
%    config = configure_performance();
%    if config.useParallel
%        fprintf('Parallel computing is available!\n');
%    end
%
% Notes:
%    - Parallel computing requires Parallel Computing Toolbox
%    - Progress bar may not display correctly in some environments
%    - Memory preallocation improves performance for large populations
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Input Arguments
    p = inputParser;
    addParameter(p, 'parallel', 'auto', @(x) islogical(x) || strcmp(x, 'auto'));
    addParameter(p, 'progress', true, @islogical);
    addParameter(p, 'prealloc', true, @islogical);
    parse(p, varargin{:});

    parallelSetting = p.Results.parallel;
    progressSetting = p.Results.progress;
    preallocSetting = p.Results.prealloc;

    %% Load Base Configuration
    baseConfig = config_iwo();
    config = baseConfig;

    fprintf('=======================================================\n');
    fprintf('  IWO Performance Configuration\n');
    fprintf('=======================================================\n\n');

    %% Configure Parallel Computing
    fprintf('Checking parallel computing availability...\n');

    hasParallelToolbox = license('test', 'Distrib_Computing_Toolbox');

    if hasParallelToolbox
        fprintf('  ✓ Parallel Computing Toolbox is available\n');

        % Check for existing pool
        poolObj = gcp('nocreate');
        if ~isempty(poolObj)
            fprintf('  ✓ Parallel pool is running (%d workers)\n', poolObj.NumWorkers);
        else
            fprintf('  - No parallel pool currently running\n');
        end

        % Set parallel option
        if strcmp(parallelSetting, 'auto')
            config.useParallel = true;
            fprintf('  → Parallel computing ENABLED (auto-detected)\n');
        elseif parallelSetting
            config.useParallel = true;
            fprintf('  → Parallel computing ENABLED (user requested)\n');
        else
            config.useParallel = false;
            fprintf('  → Parallel computing DISABLED (user requested)\n');
        end
    else
        fprintf('  ✗ Parallel Computing Toolbox not available\n');
        config.useParallel = false;
        fprintf('  → Parallel computing DISABLED\n');

        if parallelSetting == true
            warning('PerformanceConfig:NoParallelToolbox', ...
                    'Parallel computing requested but toolbox not available');
        end
    end

    fprintf('\n');

    %% Configure Progress Tracking
    fprintf('Progress tracking settings:\n');
    config.showProgress = progressSetting;

    if config.showProgress
        fprintf('  → Progress bar ENABLED\n');
        fprintf('    Updates every %d iterations\n', config.progressUpdateInterval);
    else
        fprintf('  → Progress bar DISABLED\n');
    end

    fprintf('\n');

    %% Configure Memory Optimization
    fprintf('Memory optimization settings:\n');
    config.preallocateOffspring = preallocSetting;

    if config.preallocateOffspring
        fprintf('  → Memory preallocation ENABLED\n');
        fprintf('    Reduces dynamic array growth overhead\n');
    else
        fprintf('  → Memory preallocation DISABLED\n');
    end

    fprintf('\n');

    %% Performance Estimates
    fprintf('Expected performance characteristics:\n');

    % Estimate speedup from parallelization
    if config.useParallel && ~isempty(poolObj)
        estimatedSpeedup = poolObj.NumWorkers * 0.7;  % ~70% efficiency
        fprintf('  - Parallel speedup: ~%.1fx faster (%d workers)\n', ...
                estimatedSpeedup, poolObj.NumWorkers);
    elseif config.useParallel
        fprintf('  - Parallel speedup: Pool will start on first run\n');
    end

    % Estimate memory usage
    popSize = config.maxPopSize;
    nVar = config.nVar;
    if config.preallocateOffspring
        maxOffspring = popSize * config.maxSeeds;
        estimatedMemory = (popSize + maxOffspring) * nVar * 8 / 1024^2;  % MB
        fprintf('  - Estimated peak memory: ~%.1f MB\n', estimatedMemory);
    end

    fprintf('\n');

    %% Recommendations
    fprintf('Recommendations:\n');

    if hasParallelToolbox && ~config.useParallel
        fprintf('  ℹ Consider enabling parallel computing for faster optimization\n');
    end

    if ~hasParallelToolbox
        fprintf('  ℹ Install Parallel Computing Toolbox for significant speedup\n');
    end

    if config.showProgress
        fprintf('  ℹ Progress bar may slow down optimization slightly\n');
        fprintf('    Disable with: configure_performance(''progress'', false)\n');
    end

    fprintf('\n');

    %% Save Configuration (Optional)
    fprintf('Current configuration:\n');
    fprintf('  useParallel: %s\n', mat2str(config.useParallel));
    fprintf('  showProgress: %s\n', mat2str(config.showProgress));
    fprintf('  preallocateOffspring: %s\n', mat2str(config.preallocateOffspring));
    fprintf('\n');

    fprintf('To use this configuration, either:\n');
    fprintf('  1. Modify config_iwo.m directly, or\n');
    fprintf('  2. Pass config structure to a custom IWO variant\n');
    fprintf('=======================================================\n');
end
