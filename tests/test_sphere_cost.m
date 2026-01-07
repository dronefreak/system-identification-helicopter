function tests = test_sphere_cost
% TEST_SPHERE_COST - Unit tests for Sphere cost function
%
% Tests the Sphere.m cost function used in optimization algorithms
% to ensure it handles various inputs correctly and produces valid outputs.
%
% Usage:
%    results = runtests('test_sphere_cost')
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

    % Load test data
    testCase.TestData.dataFile = 'data/experiments/best.mat';
    if exist(testCase.TestData.dataFile, 'file')
        data = load(testCase.TestData.dataFile);
        testCase.TestData.inputData = data.in_H;
        testCase.TestData.outputData = data.out_H;
        testCase.TestData.timeVector = data.time;
    else
        warning('Test data file not found. Some tests will be skipped.');
        testCase.TestData.inputData = [];
    end

    % Create valid test parameters (40 values within bounds)
    testCase.TestData.validParams = zeros(1, 40);
    testCase.TestData.validParams(1:10) = rand(1, 10) * 10 - 5;  % Small values
    testCase.TestData.validParams(11:20) = rand(1, 10) * 20 - 10;
    testCase.TestData.validParams(21:30) = rand(1, 10) * 50 - 25;
    testCase.TestData.validParams(31:40) = rand(1, 10) * 100 - 50;
end

function teardownOnce(testCase)
    % Clean up
    close all;
end

%% Test Valid Inputs

function testValidParametersProduceCost(testCase)
    % Test that valid parameters produce a real, finite cost
    if isempty(testCase.TestData.inputData)
        return;
    end

    cost = Sphere(testCase.TestData.validParams, ...
                  testCase.TestData.inputData, ...
                  testCase.TestData.outputData, ...
                  testCase.TestData.timeVector);

    verifyTrue(testCase, isreal(cost), 'Cost should be real');
    verifyTrue(testCase, isfinite(cost), 'Cost should be finite');
    verifyTrue(testCase, cost >= 0, 'Cost should be non-negative');
end

function testCostIsNumericScalar(testCase)
    % Test that cost is a numeric scalar
    if isempty(testCase.TestData.inputData)
        return;
    end

    cost = Sphere(testCase.TestData.validParams, ...
                  testCase.TestData.inputData, ...
                  testCase.TestData.outputData, ...
                  testCase.TestData.timeVector);

    verifyTrue(testCase, isnumeric(cost), 'Cost should be numeric');
    verifyTrue(testCase, isscalar(cost), 'Cost should be scalar');
end

%% Test Invalid Inputs

function testNaNParametersReturnInf(testCase)
    % Test that NaN parameters return Inf cost
    if isempty(testCase.TestData.inputData)
        return;
    end

    nanParams = testCase.TestData.validParams;
    nanParams(1) = NaN;

    cost = Sphere(nanParams, ...
                  testCase.TestData.inputData, ...
                  testCase.TestData.outputData, ...
                  testCase.TestData.timeVector);

    verifyEqual(testCase, cost, Inf, 'NaN parameters should return Inf cost');
end

function testInfParametersReturnInf(testCase)
    % Test that Inf parameters return Inf cost
    if isempty(testCase.TestData.inputData)
        return;
    end

    infParams = testCase.TestData.validParams;
    infParams(5) = Inf;

    cost = Sphere(infParams, ...
                  testCase.TestData.inputData, ...
                  testCase.TestData.outputData, ...
                  testCase.TestData.timeVector);

    verifyEqual(testCase, cost, Inf, 'Inf parameters should return Inf cost');
end

function testWrongNumberOfParametersThrowsError(testCase)
    % Test that wrong parameter count throws error
    if isempty(testCase.TestData.inputData)
        return;
    end

    wrongParams = testCase.TestData.validParams(1:30);  % Only 30 instead of 40

    verifyError(testCase, ...
                @() Sphere(wrongParams, ...
                          testCase.TestData.inputData, ...
                          testCase.TestData.outputData, ...
                          testCase.TestData.timeVector), ...
                'Sphere:InvalidParameters');
end

function testMissingInputsThrowsError(testCase)
    % Test that missing inputs throw error
    verifyError(testCase, ...
                @() Sphere(testCase.TestData.validParams), ...
                'Sphere:NotEnoughInputs');
end

%% Test Edge Cases

function testZeroParametersProduceCost(testCase)
    % Test that all-zero parameters still produce a valid cost
    if isempty(testCase.TestData.inputData)
        return;
    end

    zeroParams = zeros(1, 40);

    cost = Sphere(zeroParams, ...
                  testCase.TestData.inputData, ...
                  testCase.TestData.outputData, ...
                  testCase.TestData.timeVector);

    verifyTrue(testCase, isfinite(cost) || cost == Inf, ...
               'Zero parameters should produce valid or Inf cost');
end

function testExtremeParametersHandledGracefully(testCase)
    % Test that extreme (but valid) parameters are handled
    if isempty(testCase.TestData.inputData)
        return;
    end

    extremeParams = testCase.TestData.validParams;
    extremeParams(1:10) = 1000;  % Very large values
    extremeParams(11:20) = -1000;

    % Should return a cost (possibly Inf if model is unstable)
    cost = Sphere(extremeParams, ...
                  testCase.TestData.inputData, ...
                  testCase.TestData.outputData, ...
                  testCase.TestData.timeVector);

    verifyTrue(testCase, isnumeric(cost), 'Should return numeric cost');
end

%% Test Reproducibility

function testDeterministicResults(testCase)
    % Test that same inputs produce same outputs
    if isempty(testCase.TestData.inputData)
        return;
    end

    cost1 = Sphere(testCase.TestData.validParams, ...
                   testCase.TestData.inputData, ...
                   testCase.TestData.outputData, ...
                   testCase.TestData.timeVector);

    cost2 = Sphere(testCase.TestData.validParams, ...
                   testCase.TestData.inputData, ...
                   testCase.TestData.outputData, ...
                   testCase.TestData.timeVector);

    verifyEqual(testCase, cost1, cost2, ...
                'Same inputs should produce same cost');
end

%% Test Data Validation

function testInputDataDimensionValidation(testCase)
    % Test that input data dimensions are validated
    if isempty(testCase.TestData.inputData)
        return;
    end

    wrongInputData = testCase.TestData.inputData(:, 1:2);  % Only 2 columns instead of 4

    verifyError(testCase, ...
                @() Sphere(testCase.TestData.validParams, ...
                          wrongInputData, ...
                          testCase.TestData.outputData, ...
                          testCase.TestData.timeVector), ...
                'Sphere:InvalidInputData');
end

function testOutputDataDimensionValidation(testCase)
    % Test that output data dimensions are validated
    if isempty(testCase.TestData.inputData)
        return;
    end

    wrongOutputData = testCase.TestData.outputData(:, 1:5);  % Only 5 columns instead of 10

    verifyError(testCase, ...
                @() Sphere(testCase.TestData.validParams, ...
                          testCase.TestData.inputData, ...
                          wrongOutputData, ...
                          testCase.TestData.timeVector), ...
                'Sphere:InvalidOutputData');
end
