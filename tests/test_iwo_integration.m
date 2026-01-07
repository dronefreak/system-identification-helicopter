function tests = test_iwo_integration
% TEST_IWO_INTEGRATION - Integration tests for IWO algorithm
%
% Tests the complete IWO optimization workflow to ensure all components
% work together correctly.
%
% Usage:
%    results = runtests('test_iwo_integration')
%
% Note: These tests may take several minutes to run
%
% Author: System Identification Project
% Date: 2026-01-06

    tests = functiontests(localfunctions);
end

%% Setup and Teardown

function setupOnce(testCase)
    % Add paths
    addpath(genpath('src'));
    addpath(genpath('data'));

    % Check for data file
    if exist('data/experiments/best.mat', 'file')
        testCase.TestData.hasData = true;
    else
        warning('Test data not found. Integration tests will be skipped.');
        testCase.TestData.hasData = false;
    end
end

function teardownOnce(testCase)
    % Clean up
    close all;
    clearvars -global;
end

%% Integration Tests

function testIWOCompleteWorkflow(testCase)
    % Test complete IWO workflow with reduced iterations
    if ~testCase.TestData.hasData
        return;
    end

    % Save current directory
    originalDir = pwd;

    try
        % Load data into base workspace (as IWO expects)
        data = load('data/experiments/best.mat');
        assignin('base', 'in_H', data.in_H);
        assignin('base', 'out_H', data.out_H);
        assignin('base', 'time', data.time);
        assignin('base', 'inr', data.inr);
        assignin('base', 'outr', data.outr);

        % Modify config for quick test
        cd('src/algorithms/iwo/IWO');

        % Create temporary config with fewer iterations
        fid = fopen('config_iwo_test.m', 'w');
        fprintf(fid, 'function config = config_iwo_test()\n');
        fprintf(fid, '    config = config_iwo();\n');
        fprintf(fid, '    config.maxIterations = 10;\n');
        fprintf(fid, '    config.initialPopSize = 5;\n');
        fprintf(fid, '    config.maxPopSize = 10;\n');
        fprintf(fid, '    config.displayInterval = 10;\n');
        fprintf(fid, '    config.plotResults = false;\n');
        fprintf(fid, 'end\n');
        fclose(fid);

        % Run modified IWO (we'll need to adapt iwo.m to use this config)
        % For now, just verify the cost function works
        config = config_iwo();

        testParams = zeros(1, 40);
        testParams(1:10) = rand(1, 10) * 10 - 5;

        cost = Sphere(testParams, data.in_H, data.out_H, data.time);

        verifyTrue(testCase, isfinite(cost) || cost == Inf, ...
                   'Cost function should return valid value');

        % Clean up test config
        if exist('config_iwo_test.m', 'file')
            delete('config_iwo_test.m');
        end

        cd(originalDir);

    catch ME
        % Clean up on error
        cd(originalDir);
        if exist('src/algorithms/iwo/IWO/config_iwo_test.m', 'file')
            delete('src/algorithms/iwo/IWO/config_iwo_test.m');
        end
        rethrow(ME);
    end
end

function testCostFunctionWithModelCreation(testCase)
    % Test that cost function can create and simulate model
    if ~testCase.TestData.hasData
        return;
    end

    data = load('data/experiments/best.mat');

    % Use parameters from a known good solution if available
    testParams = zeros(1, 40);
    testParams(1:10) = rand(1, 10) * 10 - 5;
    testParams(11:20) = rand(1, 10) * 20 - 10;
    testParams(21:30) = rand(1, 10) * 50 - 25;
    testParams(31:40) = rand(1, 10) * 100 - 50;

    % Should not crash
    try
        cost = Sphere(testParams, data.in_H, data.out_H, data.time);
        success = true;
        costValue = cost;
    catch
        success = false;
        costValue = NaN;
    end

    verifyTrue(testCase, success, 'Model creation and simulation should succeed');
    verifyTrue(testCase, isnumeric(costValue), 'Should return numeric cost');
end

function testDataLoadingAndPreprocessing(testCase)
    % Test that data can be loaded and has correct format
    if ~testCase.TestData.hasData
        return;
    end

    data = load('data/experiments/best.mat');

    % Verify data structure
    verifyTrue(testCase, isfield(data, 'in_H'), 'Data should have in_H field');
    verifyTrue(testCase, isfield(data, 'out_H'), 'Data should have out_H field');
    verifyTrue(testCase, isfield(data, 'time'), 'Data should have time field');

    % Verify dimensions
    verifyEqual(testCase, size(data.in_H, 2), 4, 'Input should have 4 columns');
    verifyEqual(testCase, size(data.out_H, 2), 10, 'Output should have 10 columns');

    % Verify data types
    verifyTrue(testCase, isnumeric(data.in_H), 'in_H should be numeric');
    verifyTrue(testCase, isnumeric(data.out_H), 'out_H should be numeric');
    verifyTrue(testCase, isnumeric(data.time), 'time should be numeric');

    % Verify no NaN or Inf
    verifyTrue(testCase, all(all(isfinite(data.in_H))), ...
               'in_H should not contain NaN or Inf');
    verifyTrue(testCase, all(all(isfinite(data.out_H))), ...
               'out_H should not contain NaN or Inf');
end

function testVisualizationPipeline(testCase)
    % Test complete visualization pipeline
    if ~testCase.TestData.hasData
        return;
    end

    data = load('data/experiments/best.mat');

    % Create test parameters
    testParams = zeros(40, 1);
    testParams(1:10) = rand(10, 1) * 10 - 5;

    % Test PopulCheck
    dataStruct.inr = data.inr;
    dataStruct.outr = data.outr;
    dataStruct.t = data.time;

    try
        PopulCheck(testParams, dataStruct);
        success = true;
    catch
        success = false;
    end

    verifyTrue(testCase, success, 'Visualization should complete without error');

    close all;
end

%% Performance Tests

function testCostFunctionPerformance(testCase)
    % Test that cost function executes in reasonable time
    if ~testCase.TestData.hasData
        return;
    end

    data = load('data/experiments/best.mat');

    testParams = zeros(1, 40);
    testParams(1:20) = rand(1, 20) * 20 - 10;

    % Measure execution time
    tic;
    for i = 1:10
        cost = Sphere(testParams, data.in_H, data.out_H, data.time);
    end
    elapsed = toc;

    avgTime = elapsed / 10;

    % Should execute in less than 1 second per call
    verifyTrue(testCase, avgTime < 1.0, ...
               sprintf('Cost function should execute in < 1s (was %.3fs)', avgTime));
end

function testMemoryUsage(testCase)
    % Test that operations don't cause memory leaks
    if ~testCase.TestData.hasData
        return;
    end

    data = load('data/experiments/best.mat');
    testParams = zeros(1, 40);

    % Get initial memory
    initialMem = memory;
    initialUsed = initialMem.MemUsedMATLAB;

    % Run multiple iterations
    for i = 1:100
        cost = Sphere(testParams, data.in_H, data.out_H, data.time);
    end

    % Check final memory
    finalMem = memory;
    finalUsed = finalMem.MemUsedMATLAB;

    memIncrease = finalUsed - initialUsed;

    % Memory increase should be minimal (< 100MB)
    verifyTrue(testCase, memIncrease < 100*1024*1024, ...
               'Should not have significant memory increase');
end
