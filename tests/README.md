# Tests Directory

This directory contains test files for validating the helicopter system identification code.

## Test Coverage

Implemented tests:

### Unit Tests

1. **`test_sphere_cost.m`** - Cost Function Tests ✅
   - Valid input handling
   - NaN/Inf parameter handling
   - Input validation
   - Data dimension checking
   - Deterministic results
   - Edge cases (zero parameters, extreme values)

2. **`test_config_iwo.m`** - Configuration Tests ✅
   - Configuration structure validation
   - Parameter value checks
   - Population sizes validation
   - Seed range validation
   - Sigma range validation
   - Bounds validation
   - Physical constants

3. **`test_popul_check.m`** - Visualization Tests ✅
   - Parameter input handling
   - Figure generation
   - Invalid input handling
   - Missing data fields
   - Dimension mismatch handling
   - Default cpop handling

### Integration Tests

4. **`test_iwo_integration.m`** - IWO Algorithm Tests ✅
   - Complete workflow testing
   - Model creation and simulation
   - Data loading and preprocessing
   - Visualization pipeline
   - Performance benchmarks
   - Memory usage tests

### Regression Tests

5. **`test_regression.m`** - Baseline Comparison ✅
   - Cost function consistency
   - Configuration stability
   - Parameter bounds stability
   - Model structure consistency
   - Data format validation
   - Performance regression
   - Memory footprint
   - Output range validation

### Validation Scripts

6. **`validate_installation.m`** - Installation Validator ✅
   - MATLAB version check
   - Toolbox availability
   - Project structure verification
   - Required files check
   - Data file integrity
   - Core function testing

7. **`run_all_tests.m`** - Test Runner ✅
   - Execute all test suites
   - Generate comprehensive reports
   - Quick mode (unit tests only)
   - Full mode (all tests)
   - Verbose output option

## Running Tests

### Quick Start

```matlab
% Navigate to project root
cd /path/to/system-identification-helicopter

% Run all tests
cd tests
run_all_tests

% Run only fast unit tests
run_all_tests('quick')

% Run with verbose output
run_all_tests('verbose')

% Validate installation
validate_installation
```

### Running Individual Tests

```matlab
% Run specific test file
results = runtests('test_sphere_cost')

% Run specific test file with detailed output
results = runtests('test_sphere_cost', 'OutputDetail', 'Detailed')

% Run all tests in directory
results = runtests('tests/')
```

### Understanding Results

Test results are saved automatically to:
```
tests/test_results_YYYY-MM-DD_HH-MM-SS.mat
```

Load and analyze:
```matlab
load('test_results_2026-01-06_15-30-00.mat')
disp(results.summary)
```

## Test Data

Small test datasets will be stored in `tests/data/` for:
- Quick validation
- Edge case testing
- Regression testing

## Writing Tests

When writing tests:

1. Use MATLAB's testing framework
2. Include both positive and negative test cases
3. Test edge cases (NaN, Inf, empty inputs)
4. Mock external dependencies when possible
5. Keep tests fast (<1 second per test when possible)
6. Document test purpose

## Test Template

```matlab
function tests = test_example
    tests = functiontests(localfunctions);
end

function testBasicFunctionality(testCase)
    % Test basic functionality
    expected = 10;
    actual = myFunction(5);
    verifyEqual(testCase, actual, expected);
end

function testInvalidInput(testCase)
    % Test error handling
    verifyError(testCase, @()myFunction([]), 'MATLAB:expectedNonempty');
end

function testNaNHandling(testCase)
    % Test NaN handling
    result = myFunction(NaN);
    verifyEqual(testCase, result, Inf);
end
```

## Continuous Integration

When CI is set up:
- Tests run automatically on every commit
- Pull requests require passing tests
- Coverage reports generated
- Performance benchmarks tracked

## Contributing Tests

To add tests:

1. Create test file with `test_` prefix
2. Follow testing best practices
3. Ensure tests are deterministic
4. Update this README
5. Submit via pull request
