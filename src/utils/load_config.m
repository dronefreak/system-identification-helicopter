function config = load_config(configName)
% LOAD_CONFIG - Load a saved configuration by name
%
% This function loads a previously saved configuration for IWO optimization.
% Configurations are stored in the configs/ directory.
%
% Usage:
%    config = load_config('default')
%    config = load_config('fast_test')
%    config = load_config('high_accuracy')
%
% Inputs:
%    configName - Name of the configuration to load (without .mat extension)
%
% Outputs:
%    config     - Configuration structure
%
% See also: save_config, config_iwo
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Validate Input
    if nargin < 1
        error('LoadConfig:NoInput', 'Configuration name is required');
    end

    %% Determine Configuration Path
    % Get project root directory
    [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(scriptPath));
    configDir = fullfile(projectRoot, 'configs');

    % Build full path to config file
    if ~endsWith(configName, '.mat')
        configName = [configName '.mat'];
    end
    configFile = fullfile(configDir, configName);

    %% Load Configuration
    if ~exist(configFile, 'file')
        error('LoadConfig:FileNotFound', ...
              'Configuration file not found: %s\nAvailable configs:\n%s', ...
              configFile, list_available_configs());
    end

    try
        loadedData = load(configFile);

        if isfield(loadedData, 'config')
            config = loadedData.config;
        else
            error('LoadConfig:InvalidFormat', ...
                  'Configuration file does not contain a "config" structure');
        end

        fprintf('✓ Loaded configuration: %s\n', configName);

        % Display key parameters
        if isfield(config, 'maxIterations')
            fprintf('  Iterations: %d\n', config.maxIterations);
        end
        if isfield(config, 'initialPopSize')
            fprintf('  Population: %d-%d\n', config.initialPopSize, config.maxPopSize);
        end
        if isfield(config, 'useParallel')
            fprintf('  Parallel: %s\n', mat2str(config.useParallel));
        end

    catch ME
        error('LoadConfig:LoadError', ...
              'Error loading configuration: %s', ME.message);
    end
end

function configList = list_available_configs()
    % List available configuration files
    [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(scriptPath));
    configDir = fullfile(projectRoot, 'configs');

    if ~exist(configDir, 'dir')
        configList = '(No configs directory found)';
        return;
    end

    files = dir(fullfile(configDir, '*.mat'));

    if isempty(files)
        configList = '(No configuration files found)';
    else
        configList = '';
        for i = 1:length(files)
            [~, name, ~] = fileparts(files(i).name);
            configList = sprintf('%s  - %s\n', configList, name);
        end
    end
end
