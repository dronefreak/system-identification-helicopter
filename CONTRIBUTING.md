# Contributing to System Identification for Unmanned Helicopter

First off, thank you for considering contributing to this project! It's people like you that make this research more valuable to the community.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Enhancements](#suggesting-enhancements)
  - [Adding New Algorithms](#adding-new-algorithms)
  - [Improving Documentation](#improving-documentation)
  - [Code Contributions](#code-contributions)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Commit Messages](#commit-messages)
- [Pull Request Process](#pull-request-process)
- [Testing Guidelines](#testing-guidelines)

---

## Code of Conduct

This project and everyone participating in it is governed by basic principles of respect and professionalism. By participating, you are expected to uphold these standards:

- Use welcoming and inclusive language
- Be respectful of differing viewpoints and experiences
- Gracefully accept constructive criticism
- Focus on what is best for the community
- Show empathy towards other community members

---

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates. When you create a bug report, include as many details as possible:

**Bug Report Template:**

```markdown
**Description**
A clear description of the bug.

**To Reproduce**
Steps to reproduce the behavior:
1. Navigate to '...'
2. Run script '....'
3. See error

**Expected Behavior**
What you expected to happen.

**Actual Behavior**
What actually happened.

**Environment:**
- MATLAB Version: [e.g., R2020b]
- Operating System: [e.g., Windows 10, macOS 12.0]
- Toolboxes Installed: [e.g., Control System Toolbox, Statistics Toolbox]

**Error Messages**
```matlab
% Paste full error message here
```

**Additional Context**
Any other information that might be helpful.
```

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion:

1. **Use a clear and descriptive title**
2. **Provide a detailed description** of the suggested enhancement
3. **Explain why this enhancement would be useful** to most users
4. **List any alternatives** you've considered
5. **Include code examples** if applicable

### Adding New Algorithms

We welcome implementations of new optimization algorithms! To add a new algorithm:

1. **Create a new directory** under the root with a descriptive name
2. **Follow the existing structure**:
   - Main algorithm file (e.g., `algorithm_name.m`)
   - Cost function adapter (use `Sphere.m` pattern)
   - README in the directory explaining the algorithm
   - License file if using third-party code

3. **Standardize the interface**:
   ```matlab
   % Standard parameters all algorithms should accept
   CostFunction = @(x) Sphere(x, in, ou, t);
   nVar = 40;
   VarMin = [...];  % Same bounds as other algorithms
   VarMax = [...];
   MaxIt = 5000;    % Configurable iterations
   ```

4. **Document the algorithm**:
   - Add header comments explaining the algorithm
   - Reference the original paper
   - Explain parameter choices
   - Add to main README.md

5. **Test thoroughly**:
   - Run on the standard helicopter problem
   - Compare results with existing algorithms
   - Document performance characteristics

### Improving Documentation

Documentation improvements are always welcome! This includes:

- Fixing typos or grammatical errors
- Clarifying confusing sections
- Adding examples
- Improving code comments
- Translating documentation
- Adding diagrams or visualizations

**Small fixes**: Feel free to submit a PR directly.

**Large changes**: Open an issue first to discuss the approach.

### Code Contributions

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/your-feature-name`
3. **Make your changes** following our coding standards
4. **Test your changes** thoroughly
5. **Commit with clear messages** (see [Commit Messages](#commit-messages))
6. **Push to your fork**: `git push origin feature/your-feature-name`
7. **Open a Pull Request**

---

## Development Setup

### Prerequisites

- MATLAB R2014b or later
- Control System Toolbox
- Statistics and Machine Learning Toolbox (optional)
- Git

### Setting Up Your Development Environment

```bash
# 1. Fork and clone
git clone https://github.com/YOUR-USERNAME/system-identification-helicopter.git
cd system-identification-helicopter

# 2. Add upstream remote
git remote add upstream https://github.com/dronefreak/system-identification-helicopter.git

# 3. Create a branch
git checkout -b feature/your-feature-name

# 4. Keep your fork synced
git fetch upstream
git merge upstream/main
```

### Running Tests

```matlab
% Navigate to test directory (when available)
cd tests

% Run all tests
run_all_tests

% Run specific test
test_iwo_optimization
```

---

## Coding Standards

### MATLAB Style Guide

Follow these conventions for consistency:

#### 1. Naming Conventions

```matlab
% Variables: camelCase
numberOfIterations = 100;
bestSolution = [];
inputData = load('data.mat');

% Functions: camelCase
function result = calculateCost(parameters)
    % ...
end

% Constants: UPPER_CASE
MAX_ITERATIONS = 5000;
DEFAULT_POPULATION = 40;

% File names: match function name
% calculateCost.m contains function calculateCost()
```

#### 2. Code Formatting

```matlab
% Indentation: 4 spaces (no tabs)
for i = 1:nPop
    if condition
        doSomething();
    end
end

% Spaces around operators
x = a + b;
result = (x * 2) / 3;

% Function spacing
function output = myFunction(input1, input2)
    % Function description

    % Initialize
    output = [];

    % Process
    for i = 1:length(input1)
        output(i) = input1(i) + input2(i);
    end
end
```

#### 3. Comments

```matlab
% File header comment
% FILENAME - Brief description of what the file does
%
% Syntax: output = functionName(input1, input2)
%
% Inputs:
%    input1 - Description (type, dimensions)
%    input2 - Description (type, dimensions)
%
% Outputs:
%    output - Description (type, dimensions)
%
% Example:
%    result = functionName([1 2 3], 5);
%
% Author: Your Name
% Date: YYYY-MM-DD

% Inline comments
% Explain WHY, not WHAT
% Good:
% Use correlation coefficient because it's scale-invariant
correlation = corrcoef(y1, y2);

% Bad:
% Calculate correlation
correlation = corrcoef(y1, y2);
```

#### 4. Code Organization

```matlab
%% Main Section Title
% Description of this section

%% Initialize Variables
nVar = 40;
MaxIt = 5000;

%% Main Algorithm Loop
for iter = 1:MaxIt
    % Algorithm implementation
end

%% Results Visualization
figure;
plot(results);
```

#### 5. Error Handling

```matlab
% Validate inputs
function output = processData(inputData)
    % Check input arguments
    if nargin < 1
        error('processData:missingInput', 'Input data is required');
    end

    if ~isnumeric(inputData)
        error('processData:invalidType', 'Input must be numeric');
    end

    if any(isnan(inputData))
        warning('processData:nanValues', 'Input contains NaN values');
    end

    % Process data
    output = someOperation(inputData);
end
```

#### 6. Avoid Code Smells

```matlab
% DON'T: Use global variables
global myData;  % Bad

% DO: Pass as parameters
function result = myFunction(myData)  % Good

% DON'T: Use magic numbers
if x > 40  % Bad

% DO: Use named constants
MAX_PARAMETERS = 40;
if x > MAX_PARAMETERS  % Good

% DON'T: Grow arrays in loops
for i = 1:1000
    array(i) = i;  % Bad (if array not pre-allocated)
end

% DO: Pre-allocate
array = zeros(1000, 1);  % Good
for i = 1:1000
    array(i) = i;
end
```

---

## Commit Messages

### Format

```
<type>(<scope>): <subject>

<body>

<footer>
```

### Types

- **feat**: New feature
- **fix**: Bug fix
- **docs**: Documentation changes
- **style**: Formatting, missing semicolons, etc.
- **refactor**: Code restructuring without changing behavior
- **perf**: Performance improvements
- **test**: Adding or modifying tests
- **chore**: Maintenance tasks

### Examples

```
feat(iwo): add adaptive sigma parameter

Implement adaptive standard deviation that adjusts based on
convergence rate. This improves exploration in early iterations
and exploitation in later iterations.

Closes #42
```

```
fix(sphere): handle NaN in cost function

Add validation to check for NaN values in state-space matrices
before simulation. Returns Inf cost for invalid parameters.

Fixes #38
```

```
docs(readme): add troubleshooting section

Include common issues and solutions for:
- Missing toolbox errors
- NaN initialization problems
- Memory errors

Related to #45
```

---

## Pull Request Process

### Before Submitting

1. **Update documentation** if you've changed APIs or added features
2. **Add tests** for new functionality
3. **Ensure all tests pass**
4. **Update CHANGELOG.md** with a note describing your changes
5. **Ensure code follows style guidelines**
6. **Remove debugging code** and commented-out sections

### PR Description Template

```markdown
## Description
Brief description of what this PR does.

## Motivation
Why is this change needed? What problem does it solve?

## Changes Made
- Change 1
- Change 2
- Change 3

## Testing
How was this tested?
- [ ] Ran on standard helicopter problem
- [ ] Compared results with existing algorithms
- [ ] Tested edge cases

## Checklist
- [ ] Code follows style guidelines
- [ ] Documentation updated
- [ ] Tests added/updated
- [ ] CHANGELOG.md updated
- [ ] No debugging code left in

## Related Issues
Closes #issue_number
```

### Review Process

1. **Maintainer review**: A maintainer will review your PR
2. **Feedback**: Address any requested changes
3. **Approval**: Once approved, a maintainer will merge
4. **Recognition**: Your contribution will be acknowledged

---

## Testing Guidelines

### Manual Testing

Before submitting:

```matlab
% 1. Test on standard problem
cd 'YPEA119 Invasive Weed Optimization/IWO'
load('best.mat')
iwo  % Should complete without errors

% 2. Verify results are reasonable
% BestSol.Cost should be < 10 for good convergence

% 3. Test visualization
cd ../..
PopulCheck  % Should generate plots without errors
```

### Performance Testing

```matlab
% Profile your code
profile on
yourFunction();
profile viewer

% Check for bottlenecks
% Optimize critical sections
```

### Regression Testing

When modifying existing code:
1. **Run before changes**: Record outputs
2. **Run after changes**: Compare outputs
3. **Ensure behavior unchanged** (unless that's the goal)

---

## Questions?

Don't hesitate to ask questions:

- **Open an issue** with the "question" label
- **Start a discussion** on GitHub Discussions
- **Contact maintainers** via GitHub

---

## Recognition

Contributors will be:
- Listed in CHANGELOG.md for their contributions
- Acknowledged in the README.md
- Credited in any resulting publications (for significant contributions)

---

Thank you for contributing to advancing helicopter system identification research! 🚁
