function tests = test_popul_check
% TEST_POPUL_CHECK - Unit tests for PopulCheck visualization function
%
% Tests the PopulCheck.m function to ensure it handles various inputs
% correctly and generates proper visualizations.
%
% Usage:
%    results = runtests('test_popul_check')
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
    if exist('data/experiments/best.mat', 'file')
        data = load('data/experiments/best.mat');
        testCase.TestData.inputData = data.inr;
        testCase.TestData.outputData = data.outr;
        testCase.TestData.timeVector = data.t;

        % Create test parameters
        testCase.TestData.params = zeros(40, 1);
        testCase.TestData.params(1:10) = rand(10, 1) * 10 - 5;
        testCase.TestData.params(11:20) = rand(10, 1) * 20 - 10;
        testCase.TestData.params(21:30) = rand(10, 1) * 50 - 25;
        testCase.TestData.params(31:40) = rand(10, 1) * 100 - 50;

        testCase.TestData.cpop = ones(13, 1);
        testCase.TestData.hasData = true;
    else
        warning('Test data file not found. Tests will be skipped.');
        testCase.TestData.hasData = false;
    end
end

function teardownOnce(testCase)
    % Close all figures
    close all;
end

%% Test Valid Inputs

function testPopulCheckWithParameters(testCase)
    % Test PopulCheck with parameters provided
    if ~testCase.TestData.hasData
        return;
    end

    % Create data structure
    data.inr = testCase.TestData.inputData;
    data.outr = testCase.TestData.outputData;
    data.t = testCase.TestData.timeVector;
    data.cpop = testCase.TestData.cpop;

    % Should not error
    try
        PopulCheck(testCase.TestData.params, data);
        success = true;
    catch
        success = false;
    end

    verifyTrue(testCase, success, 'PopulCheck should execute without error');
end

function testPopulCheckGeneratesFigure(testCase)
    % Test that PopulCheck generates a figure
    if ~testCase.TestData.hasData
        return;
    end

    initialFigCount = length(findall(0, 'Type', 'figure'));

    data.inr = testCase.TestData.inputData;
    data.outr = testCase.TestData.outputData;
    data.t = testCase.TestData.timeVector;
    data.cpop = testCase.TestData.cpop;

    PopulCheck(testCase.TestData.params, data);

    finalFigCount = length(findall(0, 'Type', 'figure'));

    verifyTrue(testCase, finalFigCount > initialFigCount, ...
               'PopulCheck should create at least one figure');

    close all;
end

%% Test Invalid Inputs

function testPopulCheckInvalidParameterLength(testCase)
    % Test that wrong parameter length throws error
    if ~testCase.TestData.hasData
        return;
    end

    wrongParams = testCase.TestData.params(1:30);  % Only 30 instead of 40

    data.inr = testCase.TestData.inputData;
    data.outr = testCase.TestData.outputData;
    data.t = testCase.TestData.timeVector;

    verifyError(testCase, ...
                @() PopulCheck(wrongParams, data), ...
                'PopulCheck:InvalidParameters');
end

function testPopulCheckMissingDataFields(testCase)
    % Test that missing data fields throw error
    if ~testCase.TestData.hasData
        return;
    end

    incompleteData.inr = testCase.TestData.inputData;
    % Missing outr and t

    verifyError(testCase, ...
                @() PopulCheck(testCase.TestData.params, incompleteData), ...
                '');  % Should throw some error
end

function testPopulCheckInvalidDataType(testCase)
    % Test that invalid data type throws error
    verifyError(testCase, ...
                @() PopulCheck(testCase.TestData.params, 'not_a_struct'), ...
                'PopulCheck:InvalidDataType');
end

%% Test Edge Cases

function testPopulCheckWithDefaultCpop(testCase)
    % Test that PopulCheck works without cpop
    if ~testCase.TestData.hasData
        return;
    end

    data.inr = testCase.TestData.inputData;
    data.outr = testCase.TestData.outputData;
    data.t = testCase.TestData.timeVector;
    % No cpop field

    try
        PopulCheck(testCase.TestData.params, data);
        success = true;
    catch
        success = false;
    end

    verifyTrue(testCase, success, ...
               'PopulCheck should work with default cpop');

    close all;
end

function testPopulCheckDimensionMismatch(testCase)
    % Test handling of dimension mismatches
    if ~testCase.TestData.hasData
        return;
    end

    data.inr = testCase.TestData.inputData(1:100, :);  % Shorter
    data.outr = testCase.TestData.outputData;
    data.t = testCase.TestData.timeVector;

    verifyError(testCase, ...
                @() PopulCheck(testCase.TestData.params, data), ...
                'PopulCheck:DimensionMismatch');
end
