# Examples Directory

This directory will contain example scripts and tutorials for using the helicopter system identification tools.

## Planned Examples

Future examples to be added:

1. **Quick Start Example**
   - `quickstart.m` - Basic usage of IWO algorithm
   - Load data, run optimization, visualize results

2. **Algorithm Comparison**
   - `compare_algorithms.m` - Compare all available algorithms
   - Statistical analysis of convergence rates

3. **Custom Cost Function**
   - `custom_cost_example.m` - Implement your own cost function
   - Tutorial on modifying optimization objectives

4. **Parameter Sensitivity**
   - `sensitivity_analysis.m` - Analyze parameter sensitivity
   - Identify which parameters most affect model fit

5. **Multi-Dataset Testing**
   - `multi_dataset_validation.m` - Test on multiple flight datasets
   - Cross-validation techniques

## Contributing Examples

To add an example:

1. Create a well-commented MATLAB script
2. Include a header with:
   - Purpose
   - Required data
   - Expected runtime
   - Expected output
3. Test thoroughly
4. Add to this README
5. Submit via pull request

## Usage

To run an example:

```matlab
% Add project to path
addpath(genpath('/path/to/system-identification-helicopter'));

% Navigate to examples directory
cd examples

% Run example
exampleScript
```

## Example Template

```matlab
%% EXAMPLE_NAME - Brief description
%
% Purpose:
%    Describe what this example demonstrates
%
% Required Data:
%    - data/experiments/best.mat
%
% Runtime: ~5 minutes
%
% Expected Output:
%    - Convergence plot
%    - Final cost value
%    - Optimized parameters
%
% Author: Your Name
% Date: YYYY-MM-DD

%% Setup
clear; clc; close all;

% Add paths
addpath('../src/algorithms/iwo/IWO');

% Load data
load('../data/experiments/best.mat');

%% Run Example
% Your code here

%% Display Results
fprintf('Example complete!\n');
```
