# Configuration Files

This directory contains saved configurations for IWO optimization experiments.

## Usage

### Saving Configurations

```matlab
% Create and save a configuration
config = config_iwo();
config.maxIterations = 1000;
save_config(config, 'my_experiment', 'Quick test configuration');
```

### Loading Configurations

```matlab
% Load a saved configuration
config = load_config('my_experiment');

% Use in IWO
load('data/experiments/best.mat')
% Modify iwo.m to use loaded config or pass it manually
```

## Pre-defined Configurations

### default.mat
- Standard IWO parameters
- 5000 iterations
- Population: 20-40
- Use this for production runs

### fast_test.mat
- Quick testing configuration
- 1000 iterations
- Population: 10-20
- Use for debugging and testing

### high_accuracy.mat
- Maximum accuracy configuration
- 10000 iterations
- Population: 40-80
- Use when computational time is not a concern

### parallel_optimized.mat
- Optimized for parallel computing
- 5000 iterations
- Population: 40-80
- Requires Parallel Computing Toolbox

## Configuration Structure

Each .mat file contains:
- `config` - Configuration structure with all IWO parameters
- `metadata` - Creation time, MATLAB version, platform, description

## Creating Custom Configurations

```matlab
% Start with default
config = config_iwo();

% Modify parameters
config.maxIterations = 3000;
config.initialPopSize = 30;
config.useParallel = true;
config.showProgress = true;

% Save with description
save_config(config, 'custom', 'My custom configuration');
```

## Configuration Parameters

Key parameters you can modify:

### Algorithm Parameters
- `nVar` - Number of variables (40 for helicopter model)
- `maxIterations` - Maximum iterations
- `initialPopSize` - Initial population size
- `maxPopSize` - Maximum population size
- `minSeeds` - Minimum seeds per plant
- `maxSeeds` - Maximum seeds per plant

### Performance Options
- `useParallel` - Enable parallel computing
- `showProgress` - Show progress bar
- `progressUpdateInterval` - Progress update frequency

### Memory Options
- `preallocateOffspring` - Preallocate arrays
- `clearIntermediateData` - Clear intermediate data

### Display Options
- `displayInterval` - Console output frequency
- `plotResults` - Generate plots
- `useSemilogy` - Use log scale for plots

### Reproducibility
- `randomSeed` - Random seed for reproducibility

## Tips

1. **Always save configurations** before long runs
2. **Use descriptive names** for configurations
3. **Add descriptions** to remember the purpose
4. **Version control** configurations for experiments
5. **Test with fast_test** before production runs

## See Also

- `src/utils/save_config.m` - Save configuration function
- `src/utils/load_config.m` - Load configuration function
- `src/algorithms/iwo/IWO/config_iwo.m` - Default configuration generator
