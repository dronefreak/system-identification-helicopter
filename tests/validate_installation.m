function status = validate_installation()
% VALIDATE_INSTALLATION - Validate project installation and dependencies
%
% This script checks that all required components are properly installed
% and configured for running the helicopter system identification code.
%
% Usage:
%    status = validate_installation()
%
% Outputs:
%    status - 0 if all checks pass, 1 if any checks fail
%
% Author: System Identification Project
% Date: 2026-01-06

    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('  Installation Validation\n');
    fprintf('=======================================================\n');
    fprintf('\n');

    allPassed = true;

    %% Check MATLAB Version
    fprintf('Checking MATLAB version...\n');
    matlabVersion = version('-release');
    matlabYear = str2double(matlabVersion(1:4));

    if matlabYear >= 2014
        fprintf('  ✓ MATLAB %s (OK)\n', matlabVersion);
    else
        fprintf('  ✗ MATLAB %s (too old, need R2014b or later)\n', matlabVersion);
        allPassed = false;
    end
    fprintf('\n');

    %% Check Required Toolboxes
    fprintf('Checking required toolboxes...\n');

    % Control System Toolbox
    if license('test', 'Control_Toolbox')
        try
            ver('control');
            fprintf('  ✓ Control System Toolbox (installed)\n');
        catch
            fprintf('  ✗ Control System Toolbox (license exists but not installed)\n');
            allPassed = false;
        end
    else
        fprintf('  ✗ Control System Toolbox (not available)\n');
        allPassed = false;
    end

    % Statistics and Machine Learning Toolbox (optional)
    if license('test', 'Statistics_Toolbox')
        try
            ver('stats');
            fprintf('  ✓ Statistics and Machine Learning Toolbox (installed)\n');
        catch
            fprintf('  ⚠ Statistics and Machine Learning Toolbox (recommended)\n');
        end
    else
        fprintf('  ⚠ Statistics and Machine Learning Toolbox (recommended but optional)\n');
    end
    fprintf('\n');

    %% Check Project Structure
    fprintf('Checking project structure...\n');

    requiredDirs = {
        'src',
        'src/algorithms',
        'src/algorithms/iwo',
        'src/utils',
        'data',
        'data/experiments',
        'docs',
        'tests'
    };

    for i = 1:length(requiredDirs)
        if exist(requiredDirs{i}, 'dir')
            fprintf('  ✓ %s/ (found)\n', requiredDirs{i});
        else
            fprintf('  ✗ %s/ (missing)\n', requiredDirs{i});
            allPassed = false;
        end
    end
    fprintf('\n');

    %% Check Required Files
    fprintf('Checking required files...\n');

    requiredFiles = {
        'src/algorithms/iwo/IWO/iwo.m',
        'src/algorithms/iwo/IWO/Sphere.m',
        'src/algorithms/iwo/IWO/config_iwo.m',
        'src/utils/PopulCheck.m',
        'data/experiments/best.mat',
        'README.md',
        'LICENSE'
    };

    for i = 1:length(requiredFiles)
        if exist(requiredFiles{i}, 'file')
            fprintf('  ✓ %s (found)\n', requiredFiles{i});
        else
            fprintf('  ✗ %s (missing)\n', requiredFiles{i});
            if contains(requiredFiles{i}, 'best.mat')
                fprintf('     Note: Flight data file required for running optimizations\n');
            end
            allPassed = false;
        end
    end
    fprintf('\n');

    %% Check Data File Integrity
    fprintf('Checking data file integrity...\n');

    if exist('data/experiments/best.mat', 'file')
        try
            data = load('data/experiments/best.mat');

            % Check required fields
            requiredFields = {'in_H', 'out_H', 'time', 'inr', 'outr'};
            fieldsOK = true;

            for i = 1:length(requiredFields)
                if isfield(data, requiredFields{i})
                    fprintf('  ✓ Field "%s" present\n', requiredFields{i});
                else
                    fprintf('  ✗ Field "%s" missing\n', requiredFields{i});
                    fieldsOK = false;
                end
            end

            if fieldsOK
                % Check dimensions
                if size(data.in_H, 2) == 4
                    fprintf('  ✓ Input data has correct dimensions (4 inputs)\n');
                else
                    fprintf('  ✗ Input data has wrong dimensions\n');
                    allPassed = false;
                end

                if size(data.out_H, 2) == 10
                    fprintf('  ✓ Output data has correct dimensions (10 outputs)\n');
                else
                    fprintf('  ✗ Output data has wrong dimensions\n');
                    allPassed = false;
                end
            else
                allPassed = false;
            end

        catch ME
            fprintf('  ✗ Error loading data file: %s\n', ME.message);
            allPassed = false;
        end
    else
        fprintf('  ⚠ Data file not found (skipping integrity check)\n');
    end
    fprintf('\n');

    %% Test Core Functions
    fprintf('Testing core functions...\n');

    % Test config loading
    try
        addpath(genpath('src'));
        config = config_iwo();
        if isstruct(config) && config.nVar == 40
            fprintf('  ✓ config_iwo() works\n');
        else
            fprintf('  ✗ config_iwo() returns invalid data\n');
            allPassed = false;
        end
    catch ME
        fprintf('  ✗ config_iwo() failed: %s\n', ME.message);
        allPassed = false;
    end

    % Test Sphere function with dummy data
    if exist('data/experiments/best.mat', 'file')
        try
            data = load('data/experiments/best.mat');
            testParams = zeros(1, 40);
            cost = Sphere(testParams, data.in_H, data.out_H, data.time);

            if isnumeric(cost)
                fprintf('  ✓ Sphere() cost function works\n');
            else
                fprintf('  ✗ Sphere() returned non-numeric value\n');
                allPassed = false;
            end
        catch ME
            fprintf('  ✗ Sphere() failed: %s\n', ME.message);
            allPassed = false;
        end
    else
        fprintf('  ⚠ Skipping Sphere() test (no data file)\n');
    end

    fprintf('\n');

    %% Final Summary
    fprintf('=======================================================\n');
    fprintf('  VALIDATION SUMMARY\n');
    fprintf('=======================================================\n');
    fprintf('\n');

    if allPassed
        fprintf('✓ ALL CHECKS PASSED\n');
        fprintf('\nInstallation is complete and ready to use.\n');
        fprintf('You can now run the optimization algorithms.\n');
        status = 0;
    else
        fprintf('✗ SOME CHECKS FAILED\n');
        fprintf('\nPlease fix the issues above before running the code.\n');
        fprintf('See README.md for installation instructions.\n');
        status = 1;
    end

    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('\n');

end
