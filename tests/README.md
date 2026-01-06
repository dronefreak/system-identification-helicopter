# Tests Directory

This directory will contain test files for validating the helicopter system identification code.

## Planned Test Coverage

Future tests to be added:

### Unit Tests

1. **Cost Function Tests**
   - `test_sphere_cost.m` - Validate Sphere cost function
   - Input validation
   - NaN/Inf handling
   - Correlation calculations

2. **State-Space Model Tests**
   - `test_model_creation.m` - Validate model construction
   - Matrix dimensions
   - Parameter bounds
   - Simulation stability

3. **Configuration Tests**
   - `test_config.m` - Validate configuration loading
   - Parameter defaults
   - Bound checking

### Integration Tests

4. **Algorithm Tests**
   - `test_iwo.m` - Test IWO algorithm end-to-end
   - `test_ga.m` - Test GA algorithm
   - `test_abc.m` - Test ABC algorithm
   - Convergence verification
   - Result consistency

5. **Visualization Tests**
   - `test_popul_check.m` - Test PopulCheck function
   - Plot generation
   - Data compatibility

### Regression Tests

6. **Baseline Comparison**
   - `test_regression.m` - Compare against known good results
   - Ensure code changes don't break functionality

## Running Tests

Future test framework setup:

```matlab
% Run all tests
runtests('tests/')

% Run specific test
runtests('tests/test_sphere_cost.m')

% Generate coverage report
coverage = testCoverage('tests/')
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
