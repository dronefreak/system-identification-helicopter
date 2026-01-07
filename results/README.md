# Results Directory

This directory stores optimization results with complete experiment tracking and metadata.

## Purpose

Automatically stores:
- Optimization results (parameters, costs, population)
- Complete experiment metadata (timestamps, configuration, system info)
- Reproducibility information (random seeds, RNG state)
- Performance metrics (execution time, convergence)
- Summary reports (text and JSON formats)

## Quick Start

### Run a Tracked Experiment

```matlab
% Run with automatic tracking
results = run_experiment('my_exp', 'default', ...
                         'description', 'Testing new parameters');

% Results automatically saved to: results/my_exp_YYYYMMDD_HHMMSS.mat
```

### Load Previous Results

```matlab
% Load experiment results
results = load_experiment_results('my_exp_20260107_143022');

% Access data
bestParams = results.bestSolution.Position;
convergence = results.bestCosts;
metadata = results.metadata;
```

## File Structure

The experiment tracking system automatically creates:

```
results/
├── experimentName_20260107_143022.mat           # Complete results (MATLAB)
├── experimentName_20260107_143022_metadata.json # Metadata (JSON)
└── experimentName_20260107_143022_summary.txt   # Human-readable summary
```

### What's Saved

**MATLAB Results File (.mat)**:
- `bestSolution` - Best parameters and cost
- `bestCosts` - Convergence history (all iterations)
- `population` - Final population
- `metadata` - Complete experiment metadata
- `config` - Configuration structure used
- `timing` - Execution time details
- `randomState` - RNG state (for reproducibility)

**Metadata JSON File (.json)**:
- Experiment name, timestamps, duration
- Cost values (initial, final, improvement)
- System information (platform, MATLAB version, cores)
- Configuration summary

**Summary Text File (.txt)**:
- Human-readable experiment summary
- All parameters, results, and system info
- Best parameters (40 values)
- File locations

## Usage Examples

### Basic Experiment

```matlab
% Run experiment with default configuration
results = run_experiment('baseline_run', 'default');
```

### Custom Configuration

```matlab
% Run with custom configuration
results = run_experiment('fast_test', 'fast_test', ...
                         'description', 'Quick validation run', ...
                         'tags', {'test', 'validation'});
```

### Using Alternative Data

```matlab
% Run with different dataset
results = run_experiment('alt_data', 'default', ...
                         'dataFile', 'data/experiments/best2.mat');
```

### Reproducible Experiment

```matlab
% Configure for reproducibility
config = config_iwo();
config.useRandomSeed = true;
config.randomSeed = 12345;
save_config(config, 'reproducible');

% Run reproducible experiment
results = run_experiment('reproducible_run', 'reproducible', ...
                         'description', 'Fixed seed for reproducibility');

% Later, reproduce exact results
results2 = run_experiment('reproduce', 'reproducible');
% Should get identical results
```

### Batch Experiments

```matlab
% Run multiple experiments with different configurations
configs = {'fast_test', 'default', 'high_accuracy'};
for i = 1:length(configs)
    experimentName = sprintf('batch_%s_%02d', configs{i}, i);
    results = run_experiment(experimentName, configs{i}, ...
                             'description', sprintf('Batch run %d', i));
end
```

### Analyzing Results

```matlab
% Load and compare multiple experiments
exp1 = load_experiment_results('baseline_run_20260107_120000');
exp2 = load_experiment_results('optimized_run_20260107_130000');

% Compare convergence
figure;
semilogy(exp1.bestCosts, 'b-', 'LineWidth', 2);
hold on;
semilogy(exp2.bestCosts, 'r-', 'LineWidth', 2);
legend('Baseline', 'Optimized');
xlabel('Iteration');
ylabel('Cost');
title('Convergence Comparison');

% Compare final results
fprintf('Baseline final cost: %.6f\n', exp1.metadata.finalCost);
fprintf('Optimized final cost: %.6f\n', exp2.metadata.finalCost);
fprintf('Improvement: %.2f%%\n', ...
        (exp1.metadata.finalCost - exp2.metadata.finalCost) / exp1.metadata.finalCost * 100);
```

## Experiment Naming

Use descriptive experiment names:
- `baseline_run` - Standard reference runs
- `param_study_01` - Parameter sensitivity studies
- `comparison_ga_iwo` - Algorithm comparisons
- `validation_dataset2` - Validation runs
- `reproduce_paper_fig3` - Paper figure reproduction

Timestamps are added automatically: `experimentName_YYYYMMDD_HHMMSS`

## .gitignore

This directory is typically excluded from git (especially for large result files).
Consider using Git LFS for version controlling important results.

## Data Management

- Keep results organized by date and algorithm
- Document experimental conditions
- Archive old results periodically
- Back up important results externally
- Delete temporary/test results regularly

## Reproducibility

### Enabling Reproducibility

```matlab
% Create reproducible configuration
config = config_iwo();
config.useRandomSeed = true;          % Enable fixed seed
config.randomSeed = 42;               % Specific seed value
config.saveRandomState = true;        % Save RNG state
save_config(config, 'reproducible');

% Run experiment
results = run_experiment('paper_fig1', 'reproducible');
```

### Reproducing Results

```matlab
% Load previous results
results = load_experiment_results('paper_fig1_20260107_143022');

% Restore random state and rerun
if isfield(results, 'randomState')
    rng(results.randomState.initial);
    % Now run optimization again - should get identical results
end
```

## Best Practices

### For Research Papers

1. **Always use fixed random seeds** for published results
2. **Save complete metadata** with all experiments
3. **Use descriptive experiment names** (e.g., `paper_fig3_baseline`)
4. **Tag experiments** for easy organization
5. **Include configuration files** in supplementary materials

### For Development

1. **Use `fast_test` configuration** for quick iterations
2. **Tag experiments** with version numbers or features tested
3. **Clean up test results** periodically
4. **Keep only successful runs** for comparison

### For Production

1. **Use tracked experiments** (`run_experiment`) not raw `iwo` calls
2. **Save all results** with metadata
3. **Document experiment purpose** in description field
4. **Review metadata** before long runs

## Advanced Features

### Custom Output Directory

```matlab
results = run_experiment('exp001', 'default', ...
                         'outputDir', 'results/paper_revision/');
```

### Without Saving

```matlab
% Run without saving (testing only)
results = run_experiment('test', 'fast_test', ...
                         'saveResults', false);
```

### Accessing Metadata

```matlab
results = load_experiment_results('exp001_20260107_143022');

% System information
fprintf('Platform: %s\n', results.metadata.platform);
fprintf('MATLAB: %s\n', results.metadata.matlabRelease);
fprintf('Cores: %d\n', results.metadata.numCores);

% Timing
fprintf('Total time: %.2f min\n', results.metadata.executionTime / 60);
fprintf('Time per iteration: %.3f sec\n', results.timing.perIteration);

% Configuration
fprintf('Population: %d-%d\n', ...
        results.config.initialPopSize, results.config.maxPopSize);
fprintf('Parallel: %s\n', mat2str(results.config.useParallel));
```

## See Also

- `src/utils/run_experiment.m` - Main experiment runner
- `src/utils/save_experiment_results.m` - Results saver
- `src/utils/load_experiment_results.m` - Results loader
- `configs/README.md` - Configuration management
- `docs/PERFORMANCE.md` - Performance optimization
