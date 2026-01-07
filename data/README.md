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

All `.mat` files are MATLAB binary format.

### Loading Data

**Recommended** - Use data loading utility:
```matlab
% Load with automatic validation
data = load_flight_data('best');

% Access data
in_H = data.in_H;
out_H = data.out_H;
time = data.time;
```

**Alternative** - Direct load:
```matlab
load('data/experiments/best.mat')
```

### Validation

Validate data files before use:
```matlab
% Comprehensive validation
result = validate_flight_data('best');

% Check validation result
if result.isValid
    fprintf('Data is valid!\n');
else
    fprintf('Validation errors found\n');
end
```

### Data Compression

For large files, use compression utilities:
```matlab
% Compress large files
compress_data('data/experiments/large_file.mat');

% High compression
compress_data('large_file.mat', 'compression', 'high');
```

## Data Management

### Git LFS

Large `.mat` files (>5MB) are configured for Git LFS tracking:
- See `.gitattributes` for LFS configuration
- Refer to `docs/GIT_LFS_SETUP.md` for setup instructions
- Current total data size: ~65MB

**Files recommended for Git LFS**:
- `11_05_09.bin-1077528.mat` (30 MB)
- `SysID405705308.mat` (22 MB)
- `best2.mat` (5.0 MB)
- `best.mat` (4.2 MB)

### Data Catalog

Complete documentation of all data files:
- See `docs/DATA_CATALOG.md` for comprehensive catalog
- Includes file sizes, descriptions, and usage examples

## Data Utilities

### Available Functions

Located in `src/utils/`:

- **load_flight_data.m** - Load and validate flight data
- **validate_flight_data.m** - Comprehensive data validation
- **compress_data.m** - Compress large .mat files

### Usage Examples

```matlab
% Load with automatic validation
data = load_flight_data('best');

% Validate existing file
result = validate_flight_data('data/experiments/best.mat');

% Compress large file
compress_data('data/experiments/large_file.mat', 'method', 'v7.3');
```

## Adding New Data

When adding new flight data:
1. Place raw data in `flight_data/raw/`
2. Process and save results in `experiments/`
3. Use descriptive filenames with dates
4. **Validate data**: Run `validate_flight_data(filename)`
5. **Document**: Update `docs/DATA_CATALOG.md`
6. **Use Git LFS** for files >5MB
7. Update this README if structure changes

## Data Privacy

All data in this directory is from research UAV flights and contains no personally identifiable information.

## See Also

- `docs/DATA_CATALOG.md` - Complete data file documentation
- `docs/GIT_LFS_SETUP.md` - Git LFS setup and usage
- `src/utils/load_flight_data.m` - Data loading utility
- `src/utils/validate_flight_data.m` - Data validation utility
- `src/utils/compress_data.m` - Data compression utility
