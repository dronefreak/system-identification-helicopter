# Results Directory

This directory is reserved for storing optimization results and analysis outputs from experimental runs.

## Purpose

Use this directory to store:
- New optimization run results
- Comparative analysis outputs
- Generated plots and figures
- Performance benchmarks
- Exported data for publications

## Structure

Suggested organization:

```
results/
├── YYYY-MM-DD_algorithm_name/  - Results from specific runs
│   ├── parameters.mat          - Optimized parameters
│   ├── convergence.png         - Convergence plot
│   ├── validation.png          - Model vs. actual comparison
│   └── summary.txt             - Run summary and statistics
│
├── comparisons/                - Algorithm comparison studies
└── publications/               - Figures and data for papers
```

## Usage

### Saving Results

```matlab
% After running optimization
resultsDir = 'results/2026-01-06_iwo_run1/';
mkdir(resultsDir);

% Save parameters
save([resultsDir 'parameters.mat'], 'BestSol', 'BestCosts');

% Save convergence plot
saveas(gcf, [resultsDir 'convergence.png']);

% Save summary
fid = fopen([resultsDir 'summary.txt'], 'w');
fprintf(fid, 'Algorithm: IWO\n');
fprintf(fid, 'Final Cost: %.6f\n', BestSol.Cost);
fprintf(fid, 'Iterations: %d\n', length(BestCosts));
fclose(fid);
```

### Loading Results

```matlab
% Load previous results
load('results/2026-01-06_iwo_run1/parameters.mat');

% Display
fprintf('Best cost: %.6f\n', BestSol.Cost);
```

## Naming Conventions

Use descriptive directory names:
- `YYYY-MM-DD_algorithm_description/`
- `comparison_ga_vs_iwo_2026-01/`
- `sensitivity_analysis_tf_parameter/`

## .gitignore

This directory is typically excluded from git (especially for large result files).
Consider using Git LFS for version controlling important results.

## Data Management

- Keep results organized by date and algorithm
- Document experimental conditions
- Archive old results periodically
- Back up important results externally
- Delete temporary/test results regularly

## Publishing Results

Before publishing:
1. Verify results are reproducible
2. Document random seeds used
3. Include algorithm parameters
4. Save high-resolution figures
5. Export data in standard formats
