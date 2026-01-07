function results = run_all_tests(varargin)
% RUN_ALL_TESTS - Execute all test suites for helicopter system identification
%
% This script runs all unit tests, integration tests, and regression tests
% for the project and generates a comprehensive report.
%
% Usage:
%    results = run_all_tests()              - Run all tests
%    results = run_all_tests('quick')       - Run only fast unit tests
%    results = run_all_tests('full')        - Run all tests including slow ones
%    results = run_all_tests('verbose')     - Run with detailed output
%
% Outputs:
%    results - Test results structure with summary statistics
%
% Examples:
%    % Run all tests
%    results = run_all_tests();
%
%    % Run quick tests only
%    results = run_all_tests('quick');
%
%    % Run with verbose output
%    results = run_all_tests('verbose');
%
% Author: System Identification Project
% Date: 2026-01-06

    %% Parse Input Arguments
    p = inputParser;
    addOptional(p, 'mode', 'full', @(x) ismember(x, {'quick', 'full', 'verbose'}));
    parse(p, varargin{:});

    mode = p.Results.mode;
    verbose = strcmp(mode, 'verbose');

    %% Setup
    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('  Helicopter System Identification - Test Suite\n');
    fprintf('=======================================================\n');
    fprintf('\n');

    % Add project to path
    projectRoot = fileparts(mfilename('fullpath'));
    cd(projectRoot);
    cd('..');  % Go to project root
    addpath(genpath('src'));
    addpath(genpath('data'));
    addpath(genpath('tests'));

    % Record start time
    startTime = tic;

    %% Define Test Suites
    allTests = {
        % Unit Tests
        'test_config_iwo.m',
        'test_sphere_cost.m',
        'test_popul_check.m',

        % Integration Tests (slower)
        'test_iwo_integration.m',

        % Regression Tests
        'test_regression.m'
    };

    quickTests = {
        'test_config_iwo.m',
        'test_sphere_cost.m',
        'test_popul_check.m'
    };

    % Select tests based on mode
    if strcmp(mode, 'quick')
        testsToRun = quickTests;
        fprintf('Running QUICK test suite (unit tests only)...\n\n');
    else
        testsToRun = allTests;
        fprintf('Running FULL test suite (all tests)...\n\n');
    end

    %% Run Tests
    allResults = [];
    testSummary = struct();
    testSummary.passed = 0;
    testSummary.failed = 0;
    testSummary.incomplete = 0;
    testSummary.total = 0;

    for i = 1:length(testsToRun)
        testFile = testsToRun{i};
        testPath = fullfile(projectRoot, testFile);

        if ~exist(testPath, 'file')
            warning('Test file not found: %s', testFile);
            continue;
        end

        fprintf('-------------------------------------------------------\n');
        fprintf('Running: %s\n', testFile);
        fprintf('-------------------------------------------------------\n');

        try
            % Run the test
            if verbose
                result = runtests(testPath, 'OutputDetail', 'Detailed');
            else
                result = runtests(testPath, 'OutputDetail', 'Concise');
            end

            % Aggregate results
            allResults = [allResults; result]; %#ok<AGROW>

            % Update summary
            testSummary.passed = testSummary.passed + sum([result.Passed]);
            testSummary.failed = testSummary.failed + sum([result.Failed]);
            testSummary.incomplete = testSummary.incomplete + sum([result.Incomplete]);
            testSummary.total = testSummary.total + length(result);

            % Display summary for this test file
            fprintf('\n');
            fprintf('  Passed:     %d/%d\n', sum([result.Passed]), length(result));
            if sum([result.Failed]) > 0
                fprintf('  FAILED:     %d\n', sum([result.Failed]));
            end
            if sum([result.Incomplete]) > 0
                fprintf('  Incomplete: %d\n', sum([result.Incomplete]));
            end
            fprintf('\n');

        catch ME
            fprintf('ERROR running %s:\n', testFile);
            fprintf('  %s\n', ME.message);
            fprintf('\n');
            testSummary.failed = testSummary.failed + 1;
        end
    end

    %% Generate Report
    totalTime = toc(startTime);

    fprintf('=======================================================\n');
    fprintf('  TEST SUMMARY\n');
    fprintf('=======================================================\n');
    fprintf('\n');
    fprintf('Total Tests:    %d\n', testSummary.total);
    fprintf('Passed:         %d (%.1f%%)\n', ...
            testSummary.passed, ...
            100 * testSummary.passed / max(testSummary.total, 1));
    fprintf('Failed:         %d\n', testSummary.failed);
    fprintf('Incomplete:     %d\n', testSummary.incomplete);
    fprintf('\n');
    fprintf('Execution Time: %.2f seconds\n', totalTime);
    fprintf('\n');

    if testSummary.failed == 0
        fprintf('✓ ALL TESTS PASSED\n');
        exitCode = 0;
    else
        fprintf('✗ SOME TESTS FAILED\n');
        exitCode = 1;
    end

    fprintf('=======================================================\n');
    fprintf('\n');

    %% Return Results
    results.summary = testSummary;
    results.details = allResults;
    results.executionTime = totalTime;
    results.exitCode = exitCode;

    % Save detailed results
    resultFile = fullfile(projectRoot, sprintf('test_results_%s.mat', ...
                                               datestr(now, 'yyyy-mm-dd_HH-MM-SS')));
    save(resultFile, 'results');
    fprintf('Detailed results saved to: %s\n\n', resultFile);

end
