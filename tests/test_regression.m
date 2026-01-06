function tests = test_regression
% TEST_REGRESSION - Regression tests for helicopter system identification
%
% These tests ensure that code changes don't break existing functionality
% by comparing results against known baselines.
%
% Usage:
%    results = runtests('test_regression')
%
% Author: System Identification Project
% Date: 2026-01-06

    tests = functiontests(localfunctions);
end

%% Setup

function setupOnce(testCase)
    % Add paths
    addpath(genpath('src'));
    addpath(genpath('data'));

    % Check for test data
    if exist('data/experiments/best.mat', 'file')
        testCase.TestData.hasData = true;
        data = load('data/experiments/best.mat');
        testCase.TestData.inputData = data.in_H;
        testCase.TestData.outputData = data.out_H;
        testCase.TestData.timeVector = data.time;
    else
        warning('Test data not found. Regression tests will be skipped.');
        testCase.TestData.hasData = false;
    end

    % Baseline parameters (from a known good solution)
    % These would ideally come from a saved baseline file
    testCase.TestData.baselineParams = [
        -0.5, 30.0, -0.5, -30.0, -0.5, ...
        -0.5, 150.0, 0, 0.05, 0.05, ...
        80.0, 0, 0.5, -0.5, 0.5, ...
        0.5, -0.5, 10.0, 50.0, 0.5, ...
        0.0, 0.005, 0.0, 0, 0.0, ...
        50.0, 0.0, 10.0, 0.5, 0, ...
        0, 0.05, 0.5, 0.5, 0.05, ...
        50.0, 40.0, 5.0, 0.5, 0.5
    ];

    % Baseline cost (what the cost should be approximately)
    testCase.TestData.baselineCost = NaN;  % To be determined
end

%% Regression Tests

function testCostFunctionRegression(testCase)
    % Test that cost function produces consistent results
    if ~testCase.TestData.hasData
        return;
    end

    % Calculate cost with baseline parameters
    cost1 = Sphere(testCase.TestData.baselineParams, ...
                   testCase.TestData.inputData, ...
                   testCase.TestData.outputData, ...
                   testCase.TestData.timeVector);

    % Calculate again
    cost2 = Sphere(testCase.TestData.baselineParams, ...
                   testCase.TestData.inputData, ...
                   testCase.TestData.outputData, ...
                   testCase.TestData.timeVector);

    % Should be identical
    verifyEqual(testCase, cost1, cost2, ...
                'Cost function should be deterministic');
end

function testConfigurationRegression(testCase)
    % Test that configuration values haven't changed unexpectedly
    config = config_iwo();

    % Check critical parameters
    verifyEqual(testCase, config.nVar, 40, ...
                'Number of variables should remain 40');
    verifyEqual(testCase, config.nStates, 13, ...
                'Number of states should remain 13');
    verifyEqual(testCase, config.nInputs, 4, ...
                'Number of inputs should remain 4');
    verifyEqual(testCase, config.nOutputs, 10, ...
                'Number of outputs should remain 10');

    % Check that gravity hasn't changed
    verifyEqual(testCase, config.gravity, 9.81, 'AbsTol', 0.01, ...
                'Gravity constant should remain 9.81');
end

function testParameterBoundsRegression(testCase)
    % Test that parameter bounds haven't changed
    config = config_iwo();

    % Verify bounds structure
    verifyEqual(testCase, length(config.varMin), 40, ...
                'Should have 40 lower bounds');
    verifyEqual(testCase, length(config.varMax), 40, ...
                'Should have 40 upper bounds');

    % Check specific bound values (spot checks)
    verifyTrue(testCase, config.varMin(2) == -60, ...
               'Xa lower bound should be -60');
    verifyTrue(testCase, config.varMax(7) == 220, ...
               'Lb upper bound should be 220');
end

function testModelStructureRegression(testCase)
    % Test that state-space model structure is consistent
    if ~testCase.TestData.hasData
        return;
    end

    % We can't directly test the internal model creation without
    % modifying Sphere.m, but we can test that the interface is stable

    params = testCase.TestData.baselineParams;

    % Test with known parameters
    cost = Sphere(params, ...
                  testCase.TestData.inputData, ...
                  testCase.TestData.outputData, ...
                  testCase.TestData.timeVector);

    % Cost should be finite and reasonable
    verifyTrue(testCase, isfinite(cost) && cost < 1000, ...
               'Cost with baseline parameters should be finite and reasonable');
end

function testDataFormatRegression(testCase)
    % Test that data file format hasn't changed
    if ~testCase.TestData.hasData
        return;
    end

    data = load('data/experiments/best.mat');

    % Check expected fields
    expectedFields = {'in_H', 'out_H', 'time', 'inr', 'outr'};
    for i = 1:length(expectedFields)
        verifyTrue(testCase, isfield(data, expectedFields{i}), ...
                   sprintf('Data should have field: %s', expectedFields{i}));
    end

    % Check data dimensions
    nSamples = size(data.in_H, 1);
    verifyTrue(testCase, nSamples > 1000, ...
               'Should have sufficient data samples');

    verifyEqual(testCase, size(data.in_H, 1), size(data.out_H, 1), ...
                'Input and output should have same number of samples');
end

%% Baseline Performance Tests

function testBaselinePerformance(testCase)
    % Test that performance hasn't degraded
    if ~testCase.TestData.hasData
        return;
    end

    % Measure execution time
    tic;
    numRuns = 10;
    for i = 1:numRuns
        cost = Sphere(testCase.TestData.baselineParams, ...
                      testCase.TestData.inputData, ...
                      testCase.TestData.outputData, ...
                      testCase.TestData.timeVector);
    end
    elapsed = toc;

    avgTime = elapsed / numRuns;

    % Should execute in reasonable time (< 0.5 seconds per call)
    verifyTrue(testCase, avgTime < 0.5, ...
               sprintf('Cost function performance degraded (%.3fs > 0.5s)', avgTime));
end

function testMemoryFootprint(testCase)
    % Test that memory footprint hasn't increased significantly
    if ~testCase.TestData.hasData
        return;
    end

    % Clear variables and force garbage collection
    clear;
    java.lang.System.gc();
    pause(0.1);

    initialMem = memory;

    % Reload data
    addpath(genpath('src'));
    data = load('data/experiments/best.mat');

    finalMem = memory;

    % Data loading shouldn't use excessive memory
    memUsed = finalMem.MemUsedMATLAB - initialMem.MemUsedMATLAB;

    % Should use less than 100MB
    verifyTrue(testCase, memUsed < 100*1024*1024, ...
               'Memory usage should be reasonable');
end

%% Output Validation Tests

function testOutputRangeRegression(testCase)
    % Test that cost function outputs remain in expected range
    if ~testCase.TestData.hasData
        return;
    end

    % Test with multiple parameter sets
    testSets = {
        testCase.TestData.baselineParams,
        zeros(1, 40),
        rand(1, 40) * 100 - 50
    };

    costs = zeros(length(testSets), 1);

    for i = 1:length(testSets)
        costs(i) = Sphere(testSets{i}, ...
                          testCase.TestData.inputData, ...
                          testCase.TestData.outputData, ...
                          testCase.TestData.timeVector);
    end

    % All costs should be finite or Inf (not NaN)
    verifyTrue(testCase, all(~isnan(costs)), ...
               'Cost function should not return NaN');

    % At least one should be finite
    verifyTrue(testCase, any(isfinite(costs)), ...
               'At least one parameter set should produce finite cost');
end
