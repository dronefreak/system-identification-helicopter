function tests = test_config_iwo
% TEST_CONFIG_IWO - Unit tests for IWO configuration
%
% Tests the config_iwo.m configuration function to ensure it returns
% valid configuration parameters.
%
% Usage:
%    results = runtests('test_config_iwo')
%
% Author: System Identification Project
% Date: 2026-01-06

    tests = functiontests(localfunctions);
end

%% Setup

function setupOnce(testCase)
    % Add paths
    addpath(genpath('src'));
end

%% Test Configuration Loading

function testConfigReturnsStruct(testCase)
    % Test that config returns a structure
    config = config_iwo();

    verifyTrue(testCase, isstruct(config), 'Config should return a struct');
end

function testConfigHasRequiredFields(testCase)
    % Test that config has all required fields
    config = config_iwo();

    requiredFields = {'nVar', 'maxIterations', 'initialPopSize', ...
                      'maxPopSize', 'varMin', 'varMax', 'gravity'};

    for i = 1:length(requiredFields)
        verifyTrue(testCase, isfield(config, requiredFields{i}), ...
                   sprintf('Config missing field: %s', requiredFields{i}));
    end
end

%% Test Parameter Values

function testNVarIsPositive(testCase)
    % Test that nVar is positive
    config = config_iwo();

    verifyTrue(testCase, config.nVar > 0, 'nVar should be positive');
    verifyEqual(testCase, config.nVar, 40, 'nVar should be 40 for helicopter model');
end

function testIterationsArePositive(testCase)
    % Test that iterations are positive
    config = config_iwo();

    verifyTrue(testCase, config.maxIterations > 0, 'maxIterations should be positive');
end

function testPopulationSizesValid(testCase)
    % Test that population sizes are valid
    config = config_iwo();

    verifyTrue(testCase, config.initialPopSize > 0, ...
               'initialPopSize should be positive');
    verifyTrue(testCase, config.maxPopSize > 0, ...
               'maxPopSize should be positive');
    verifyTrue(testCase, config.maxPopSize >= config.initialPopSize, ...
               'maxPopSize should be >= initialPopSize');
end

function testSeedsRangeValid(testCase)
    % Test that seeds range is valid
    config = config_iwo();

    verifyTrue(testCase, config.minSeeds > 0, 'minSeeds should be positive');
    verifyTrue(testCase, config.maxSeeds > config.minSeeds, ...
               'maxSeeds should be > minSeeds');
end

function testSigmaRangeValid(testCase)
    % Test that sigma range is valid
    config = config_iwo();

    verifyTrue(testCase, config.sigmaInitial > config.sigmaFinal, ...
               'sigmaInitial should be > sigmaFinal');
    verifyTrue(testCase, config.sigmaFinal > 0, ...
               'sigmaFinal should be positive');
end

%% Test Parameter Bounds

function testVarMinMaxLengthsMatch(testCase)
    % Test that varMin and varMax have correct length
    config = config_iwo();

    verifyEqual(testCase, length(config.varMin), config.nVar, ...
                'varMin length should match nVar');
    verifyEqual(testCase, length(config.varMax), config.nVar, ...
                'varMax length should match nVar');
end

function testVarMinLessThanVarMax(testCase)
    % Test that varMin < varMax for all parameters
    config = config_iwo();

    for i = 1:config.nVar
        verifyTrue(testCase, config.varMin(i) < config.varMax(i), ...
                   sprintf('varMin(%d) should be < varMax(%d)', i, i));
    end
end

function testBoundsAreFinite(testCase)
    % Test that all bounds are finite
    config = config_iwo();

    verifyTrue(testCase, all(isfinite(config.varMin)), ...
               'All varMin values should be finite');
    verifyTrue(testCase, all(isfinite(config.varMax)), ...
               'All varMax values should be finite');
end

%% Test Physical Constants

function testGravityIsReasonable(testCase)
    % Test that gravity is approximately 9.81 m/s^2
    config = config_iwo();

    verifyTrue(testCase, abs(config.gravity - 9.81) < 0.1, ...
               'Gravity should be approximately 9.81 m/s^2');
end

%% Test Cost Function Parameters

function testCostParametersValid(testCase)
    % Test that cost function parameters are valid
    config = config_iwo();

    verifyTrue(testCase, isinf(config.maxCostPenalty), ...
               'maxCostPenalty should be Inf');
    verifyTrue(testCase, config.targetOutputs > 0, ...
               'targetOutputs should be positive');
    verifyTrue(testCase, config.correlationTarget > 0, ...
               'correlationTarget should be positive');
end
