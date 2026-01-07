# Data Catalog

Comprehensive catalog of all data files in the system identification helicopter repository.

## Overview

This repository contains flight test data, experimental results, and optimization outputs for the TREX 550 flybarless helicopter system identification project.

**Total Data Size**: ~65 MB (across all `.mat` files)
**Data Types**: Flight logs, experimental results, optimization outputs
**Format**: MATLAB `.mat` files (v7.3 compatible)

---

## Primary Datasets

### 1. best.mat
**Location**: `data/experiments/best.mat`
**Size**: 4.2 MB
**Status**: ✅ Primary dataset (recommended)

**Description**: Main dataset for IWO optimization containing preprocessed flight data.

**Contents**:
- `in_H` - Helicopter input data (4 control inputs × time steps)
- `out_H` - Helicopter output data (10 measurements × time steps)
- `time` - Time vector for simulation
- `inr` - Reference input matrix
- `outr` - Reference output matrix
- Initial population data for optimization

**Usage**:
```matlab
load('data/experiments/best.mat')
% Access data: in_H, out_H, time, inr, outr
```

**Recommended For**:
- Standard IWO optimization runs
- Algorithm testing and validation
- Reproducible research

---

### 2. best2.mat
**Location**: `data/experiments/best2.mat`
**Size**: 5.0 MB
**Status**: ✅ Alternative dataset

**Description**: Alternative flight dataset with different flight maneuvers or conditions.

**Contents**: Similar structure to `best.mat`
- Input/output data from different flight session
- May have different time duration or sampling rate

**Usage**:
```matlab
load('data/experiments/best2.mat')
% Validate with different data
```

**Recommended For**:
- Cross-validation
- Testing robustness
- Comparing algorithm performance

---

## Experimental Results

### 3. 11_05_09.bin-1077528.mat
**Location**: `data/experiments/11_05_09.bin-1077528.mat`
**Size**: 30 MB ⚠️ (Largest file - use Git LFS)
**Status**: ⚠️ Large file

**Description**: Raw flight log converted from binary format.

**Date**: May 9, 2011 (based on filename)
**Origin**: Binary flight log (`.bin` file)

**Recommendation**: Use Git LFS for this file

---

### 4. SysID405705308.mat
**Location**: `data/experiments/SysID405705308.mat`
**Size**: 22 MB ⚠️ (Second largest - use Git LFS)
**Status**: ⚠️ Large file

**Description**: System identification experimental results.

**ID**: 405705308 (likely timestamp or experiment ID)

**Recommendation**: Use Git LFS for this file

---

### 5. ParamCABAc.mat
**Location**: `data/experiments/ParamCABAc.mat`
**Size**: 2.0 MB
**Status**: ✅ Medium size

**Description**: Parameter set from CABC (Chaotic Artificial Bee Colony) algorithm.

**Contents**: Optimized 40-parameter set

**Usage**:
```matlab
load('data/experiments/ParamCABAc.mat')
% Examine CABC optimization results
```

---

### 6. Record.mat
**Location**: `data/experiments/Record.mat`
**Size**: 1.3 KB (tiny)
**Status**: ✅ Small

**Description**: Record or log file, possibly containing metadata or summary statistics.

---

## Experiment Data Subdirectory

### Location: `data/experiments/Experiment_data/`

This subdirectory contains results from algorithm comparison experiments.

#### Combined.mat
**Description**: Combined results from multiple algorithm runs
**Purpose**: Comparative analysis

#### Data_graphready.mat
**Description**: Preprocessed data ready for plotting/visualization
**Purpose**: Figure generation for papers

#### GA_1000.mat
**Description**: Genetic Algorithm results (1000 iterations)
**Contents**: GA optimization output
**Usage**: Compare GA vs IWO performance

#### IWO_1000.mat
**Description**: IWO results (1000 iterations)
**Contents**: IWO optimization output
**Usage**: Algorithm comparison baseline

#### gaout.mat
**Description**: GA output data
**Purpose**: Genetic Algorithm experiment results

#### iwout.mat
**Description**: IWO output data
**Purpose**: IWO experiment results

---

## Flight Data

### Location: `data/flight_data/raw/`

Raw flight test data from the TREX 550 helicopter.

#### Binary Logs (Converted to .mat)

1. **12_10_16.bin-497185.mat** - Flight log from Dec 16, 2010
2. **12_11_59.bin-289768.mat** - Flight log from Nov 59, 2012
3. **26_11_2016.bin-349438.mat** - Flight log from Nov 26, 2016
4. **9-1-17.mat** - Flight log from Jan 9, 2017

**Format**: Converted from binary `.bin` logs to MATLAB `.mat` format

**Contents**: Time-series flight data including:
- Control inputs (4 channels)
- Sensor measurements (10 outputs)
- Timestamps

---

## Data Variable Naming Conventions

### Input Data
- `in` or `in_H` - Helicopter control inputs (4 channels)
  - Channel 1: Lateral cyclic
  - Channel 2: Longitudinal cyclic
  - Channel 3: Pedal (yaw)
  - Channel 4: Collective (thrust)

### Output Data
- `out` or `out_H` - Helicopter state measurements (10 channels)
  - Velocities: u, v, w
  - Angular rates: p, q, r
  - Euler angles: φ, θ, ψ
  - Flapping angle: a

### Time Data
- `time` or `t` - Time vector (seconds)

### Reference Data
- `inr` - Reference input matrix
- `outr` - Reference output matrix

### Parameters
- `popul` or `params` - 40-parameter vector
- `BestSol` - Best solution structure
- `BestCosts` - Convergence history

---

## Data Loading

### Quick Start

```matlab
% Load primary dataset
load('data/experiments/best.mat')

% Or use data loading utility
data = load_flight_data('best');
```

### Validation

```matlab
% Validate data structure
validate_flight_data('data/experiments/best.mat');
```

See `src/utils/load_data.m` and `src/utils/validate_data.m` for data utilities.

---

## Data File Sizes and Git LFS

| File | Size | Git LFS? | Priority |
|------|------|----------|----------|
| 11_05_09.bin-1077528.mat | 30 MB | ✅ Required | Medium |
| SysID405705308.mat | 22 MB | ✅ Required | Medium |
| best2.mat | 5.0 MB | ⚠️ Recommended | High |
| best.mat | 4.2 MB | ⚠️ Recommended | Critical |
| ParamCABAc.mat | 2.0 MB | Optional | Low |
| Record.mat | 1.3 KB | ❌ Not needed | Low |

**Legend**:
- ✅ Required: File should use Git LFS
- ⚠️ Recommended: Consider Git LFS for better performance
- ❌ Not needed: File is small enough for regular Git

---

## Data Provenance

### Flight Test Hardware
- **Aircraft**: TREX 550 flybarless helicopter
- **IMU**: Onboard inertial measurement unit
- **Control**: RC transmitter (4-channel)
- **Data Logger**: Custom flight data recorder

### Data Collection
- **Location**: Flight test facility
- **Conditions**: Various flight maneuvers
- **Sampling Rate**: Variable (check time vector)
- **Duration**: Variable per flight

### Data Processing
1. Raw binary logs (`.bin` files)
2. Conversion to MATLAB (`.mat` format)
3. Preprocessing and filtering
4. Packaging for optimization (`best.mat`)

---

## Data Quality

### Validation Checks

Run validation before using data:
```matlab
% Check data integrity
validate_flight_data('data/experiments/best.mat');

% Expected output:
% ✓ File exists and is readable
% ✓ Required variables present
% ✓ Data dimensions consistent
% ✓ No NaN or Inf values
% ✓ Time vector monotonic
```

### Known Issues

- Some flight logs may have sensor dropouts
- Time synchronization varies between logs
- Control input saturation in some maneuvers

---

## Data Usage Guidelines

### For Research Papers
1. Always cite data source
2. Document which dataset used (`best.mat` or `best2.mat`)
3. Include preprocessing details
4. Report data quality metrics
5. Make results reproducible with specific data version

### For Development
1. Use `best.mat` for standard testing
2. Validate with `best2.mat` for robustness
3. Keep test data separate from production data
4. Document any data modifications

### For Experiments
1. Use `run_experiment()` with proper data file specification
2. Save metadata linking results to data files
3. Track data provenance in experiment notes

---

## Data Backup and Archival

### Backup Strategy
- Keep original binary logs archived separately
- Version control `.mat` files with Git LFS
- Back up large files to external storage
- Maintain data checksums for integrity

### Archival
- Important datasets: Store in institutional repository
- Published data: Upload to Zenodo/Figshare with DOI
- Working data: Git LFS with regular backups

---

## Adding New Data

### Checklist
1. Convert to MATLAB `.mat` format (v7.3)
2. Document in this catalog
3. Add to appropriate directory
4. Use Git LFS if > 5 MB
5. Update data loading utilities if needed
6. Run validation checks
7. Commit with descriptive message

### Example
```matlab
% Save new dataset
save('data/experiments/new_flight_20260107.mat', ...
     'in_H', 'out_H', 'time', 'inr', 'outr', '-v7.3');

% Validate
validate_flight_data('data/experiments/new_flight_20260107.mat');

% Commit
git add data/experiments/new_flight_20260107.mat
git commit -m "data: add flight test from 2026-01-07"
```

---

## See Also

- `docs/GIT_LFS_SETUP.md` - Git LFS configuration guide
- `src/utils/load_data.m` - Data loading utilities
- `src/utils/validate_data.m` - Data validation functions
- `src/utils/compress_data.m` - Data compression utilities
- `data/README.md` - Data directory overview
- `data/experiments/README.md` - Experiments data guide
- `data/flight_data/README.md` - Flight data guide

---

## Contact

For questions about specific datasets or data collection procedures, please open an issue or contact the repository maintainers.

**Last Updated**: 2026-01-07
