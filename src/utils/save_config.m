function save_config(config, configName, description)
% SAVE_CONFIG - Save a configuration for future use
%
% This function saves an IWO configuration to the configs/ directory for
% later reuse. Useful for maintaining different experiment configurations.
%
% Usage:
%    save_config(config, 'my_experiment')
%    save_config(config, 'fast_test', 'Quick test with 1000 iterations')
%
% Inputs:
%    config      - Configuration structure (from config_iwo or modified)
%    configName  - Name for this configuration (without .mat extension)
%    description - Optional description of this configuration
%
% The saved file includes:
%    - config: The configuration structure
%    - metadata: Creation time, MATLAB version, description
%
% Examples:
%    % Save default configuration
%    config = config_iwo();
%    save_config(config, 'default', 'Default IWO settings');
%
%    % Save fast test configuration
%    config = config_iwo();
%    config.maxIterations = 1000;
%    config.initialPopSize = 10;
%    save_config(config, 'fast_test', 'Quick test - 1000 iterations');
%
%    % Save high accuracy configuration
%    config = config_iwo();
%    config.maxIterations = 10000;
%    config.maxPopSize = 80;
%    save_config(config, 'high_accuracy', 'High accuracy - 10000 iterations');
%
% See also: load_config, config_iwo
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Validate Inputs
    if nargin < 2
        error('SaveConfig:NotEnoughInputs', ...
              'Configuration and name are required');
    end

    if nargin < 3
        description = '';
    end

    if ~isstruct(config)
        error('SaveConfig:InvalidConfig', ...
              'First argument must be a configuration structure');
    end

    %% Prepare Configuration Directory
    [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(scriptPath));
    configDir = fullfile(projectRoot, 'configs');

    % Create configs directory if it doesn't exist
    if ~exist(configDir, 'dir')
        mkdir(configDir);
        fprintf('Created configs directory: %s\n', configDir);
    end

    %% Prepare Metadata
    metadata = struct();
    metadata.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    metadata.matlabVersion = version;
    metadata.description = description;
    metadata.configName = configName;

    % Add system information
    try
        if ispc
            metadata.platform = 'Windows';
        elseif ismac
            metadata.platform = 'macOS';
        elseif isunix
            metadata.platform = 'Linux';
        else
            metadata.platform = 'Unknown';
        end
    catch
        metadata.platform = 'Unknown';
    end

    %% Save Configuration
    % Remove .mat extension if provided
    if endsWith(configName, '.mat')
        configName = configName(1:end-4);
    end

    configFile = fullfile(configDir, [configName '.mat']);

    % Check if file exists and warn
    if exist(configFile, 'file')
        warning('SaveConfig:Overwriting', ...
                'Configuration "%s" already exists and will be overwritten', configName);
    end

    try
        save(configFile, 'config', 'metadata');
        fprintf('✓ Configuration saved: %s\n', configName);
        fprintf('  Location: %s\n', configFile);

        if ~isempty(description)
            fprintf('  Description: %s\n', description);
        end

        % Display key parameters
        fprintf('  Parameters:\n');
        if isfield(config, 'maxIterations')
            fprintf('    Iterations: %d\n', config.maxIterations);
        end
        if isfield(config, 'initialPopSize')
            fprintf('    Initial population: %d\n', config.initialPopSize);
        end
        if isfield(config, 'maxPopSize')
            fprintf('    Max population: %d\n', config.maxPopSize);
        end
        if isfield(config, 'useParallel')
            fprintf('    Parallel: %s\n', mat2str(config.useParallel));
        end

    catch ME
        error('SaveConfig:SaveError', ...
              'Error saving configuration: %s', ME.message);
    end
end
