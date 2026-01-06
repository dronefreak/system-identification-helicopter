# Data Directory

This directory contains all data files for the helicopter system identification project.

## Structure

```
data/
├── flight_data/    - Raw flight test data from TREX 550 helicopter
│   └── raw/       - Original .bin and .mat files from flight tests
│
├── experiments/    - Experimental results and analysis
│   ├── *.mat      - Optimization results (best.mat, best2.mat, etc.)
│   └── *.pdf      - Analysis documents and plots
│
└── logs/           - Flight and simulation logs
    ├── flight_logs/ - Organized flight sessions (sess001-010, dated logs)
    └── reference/   - Reference logs for comparison
```

## Flight Data

Raw flight test data includes:
- Control inputs (lateral, longitudinal, pedal, collective)
- Measured outputs (velocities, rates, angles)
- Time vectors
- Test conditions and metadata

### Key Files

- `experiments/best.mat` - Primary dataset with pre-defined initial population
  - Contains: `in_H`, `out_H`, `time`, `inr`, `outr`, `population`
  - Used by IWO and other algorithms

- `experiments/best2.mat` - Alternative dataset for testing

## Experiment Results

Results from optimization runs:
- Converged parameter sets
- Cost history
- Statistical analysis
- Correlation data

## Data Format

All `.mat` files are MATLAB binary format. Load with:

```matlab
load('data/experiments/best.mat')
```

## Data Usage

**Important**: Due to file sizes (some >20MB), these files may be tracked with Git LFS in future updates.

Current total data size: ~60MB

## Adding New Data

When adding new flight data:
1. Place raw data in `flight_data/raw/`
2. Process and save results in `experiments/`
3. Use descriptive filenames with dates
4. Document data collection conditions
5. Update this README

## Data Privacy

All data in this directory is from research UAV flights and contains no personally identifiable information.
