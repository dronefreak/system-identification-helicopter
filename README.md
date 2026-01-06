# System Identification for Unmanned Helicopter

[![MATLAB](https://img.shields.io/badge/MATLAB-R2014b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-success.svg)]()

> Parameter optimization for a TREX 550 flybarless helicopter using metaheuristic algorithms and state-space modeling

## Table of Contents

- [About](#about)
- [Features](#features)
- [Installation](#installation)
- [Prerequisites](#prerequisites)
- [Quick Start](#quick-start)
- [Algorithms Implemented](#algorithms-implemented)
- [Usage Guide](#usage-guide)
  - [Running IWO](#running-iwo)
  - [Running Genetic Algorithm](#running-genetic-algorithm)
  - [Running Other Algorithms](#running-other-algorithms)
- [Parameter Description](#parameter-description)
- [Project Structure](#project-structure)
- [Results Visualization](#results-visualization)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Authors](#authors)

---

## About

This project implements **system identification** for an unmanned TREX 550 flybarless helicopter using the state-space structure based on Bernard Mettler's Yamaha RMAX helicopter model. The system optimizes 40 aerodynamic and control parameters to maximize the correlation between actual flight test data and simulated model outputs.

### Objective

Maximize the **Pearson correlation coefficient** between:
- Actual flight test data from the TREX 550 helicopter
- Simulated outputs from the state-space model

The optimization problem involves a **40-dimensional parameter space** representing the A, B, C, and D matrices of the state-space model with 13 states, 4 inputs, and 10 measurable outputs.

---

## Features

- **Multiple Optimization Algorithms**: Six different metaheuristic algorithms implemented
- **Real Flight Data**: Validated using actual flight test data from TREX 550 helicopter
- **State-Space Modeling**: Based on Bernard Mettler's proven helicopter dynamics model
- **3D Visualization**: Advanced animation and visualization capabilities
- **Comprehensive Results**: Detailed comparison and statistical analysis tools
- **Modular Design**: Easy to extend with new algorithms or cost functions

---

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/dronefreak/system-identification-helicopter.git
cd system-identification-helicopter
```

### 2. Verify MATLAB Installation

Ensure you have MATLAB R2014b or later installed with the required toolboxes (see [Prerequisites](#prerequisites)).

### 3. Add to MATLAB Path

Open MATLAB and add the project directory to your path:

```matlab
addpath(genpath('/path/to/system-identification-helicopter'));
savepath;
```

---

## Prerequisites

### Software Requirements

- **MATLAB**: R2014b or later (tested on R2014b, R2016a, R2018b)
- **Operating System**: Windows, macOS, or Linux

### Required MATLAB Toolboxes

- **Control System Toolbox** - For state-space models (`ss()`, `lsim()`)
- **Statistics and Machine Learning Toolbox** - For correlation functions (recommended)

### Verify Toolbox Installation

```matlab
% Check if toolboxes are installed
ver('control')
ver('stats')
```

### Third-Party Dependencies

The following third-party libraries are included in the repository:

- **Yarpiz YPEA Framework** - Optimization algorithms (BSD License)
- **STLRead** - 3D model import (File Exchange)
- **MatlabPlaneGraphics** - 3D visualization framework

---

## Quick Start

### Run IWO Optimization (Recommended)

```matlab
% 1. Navigate to IWO directory
cd src/algorithms/iwo/IWO

% 2. Load the workspace with flight data
load('../../../../data/experiments/best.mat')

% 3. Run the optimization (5000 iterations)
iwo

% 4. Visualize results
cd ../../../..
addpath('src/utils')
PopulCheck
```

Expected runtime: **30-60 minutes** (depends on your CPU)

---

## Algorithms Implemented

| Algorithm | Status | Directory | Iterations | Population Size |
|-----------|--------|-----------|------------|-----------------|
| **Invasive Weed Optimization (IWO)** | ✅ Tested | `src/algorithms/iwo/` | 5000 | 40 |
| **Genetic Algorithm (GA)** | ✅ Tested | `src/algorithms/ga/` | 1000 | Variable |
| **Artificial Bee Colony (ABC)** | ⚠️ Experimental | `src/algorithms/abc/` | Configurable | Variable |
| **Bees Behavior Algorithm (BBA)** | ⚠️ Experimental | `src/algorithms/bba/` | Configurable | Variable |
| **Particle Swarm Optimization (PSO)** | 📝 Example | `src/algorithms/pso_example.m` | - | - |
| **Simulated Annealing (SA)** | 🔨 In Development | `src/algorithms/sa/` | Configurable | Variable |

**Legend**: ✅ Fully tested | ⚠️ Experimental | 📝 Example code | 🔨 Under development

---

## Usage Guide

### Running IWO

**Invasive Weed Optimization** is the recommended algorithm for this problem.

```matlab
% Step 1: Load the workspace
cd src/algorithms/iwo/IWO
load('../../../../data/experiments/best.mat')  % Contains: inr, outr, t, population

% Step 2: Run optimization
iwo  % Runs for 5000 iterations

% Step 3: Extract best parameters
bestParameters = BestSol.Position;

% Step 4: Visualize results
cd ../../../..
addpath('src/utils')
popul = bestParameters';  % Copy to popul variable
PopulCheck  % Generates comparison plots
```

**What happens during optimization:**
- The `Sphere.m` cost function creates a state-space model using the 40 parameters
- Simulates helicopter response using `lsim()` with actual input data
- Compares simulated output with actual flight test data
- Returns a cost value (to be minimized) based on correlation and least squares error

### Running Genetic Algorithm

```matlab
cd src/algorithms/ga
MainCode  % Run genetic algorithm optimization
```

### Running Other Algorithms

**Artificial Bee Colony:**
```matlab
cd src/algorithms/abc
main  % Run ABC optimization
```

**Bees Behavior Algorithm:**
```matlab
cd src/algorithms/bba
Main  % Run BBA optimization
```

**Simulated Annealing:**
```matlab
cd 'src/algorithms/sa/01 TSP using SA (Standard)'
% Configure and run
```

---

## Parameter Description

The optimization identifies **40 parameters** in the state-space model:

### State-Space Model Structure

```
ẋ = Ax + Bu
y = Cx + Du
```

Where:
- **States (13)**: u, v, w, p, q, r, φ, θ, ψ, a, b, rfb, c
  - u, v, w: Velocity components
  - p, q, r: Angular rates
  - φ, θ, ψ: Euler angles
  - a, b, c: Rotor flapping angles
  - rfb: Rotor feedback

- **Inputs (4)**: δlat, δlon, δped, δcol
  - Lateral cyclic, Longitudinal cyclic, Pedal, Collective

- **Outputs (10)**: Measurable states (excluding rfb, c, d)

### Parameter Ranges

| Parameter Type | Min Value | Max Value | Description |
|---------------|-----------|-----------|-------------|
| Velocity derivatives | -60 | +60 | Xu, Xa, Yv, Yb, etc. |
| Rotational dynamics | -220 | +220 | Angular rate coefficients |
| Time constants | 0 | 220 | Tf, Ts |
| Control gains | -100 | +100 | Control effectiveness |

**Total**: 40 parameters optimized simultaneously

---

## Project Structure

```
system-identification-helicopter/
├── src/                                  # Source code
│   ├── algorithms/                       # Optimization algorithms
│   │   ├── iwo/                         # Invasive Weed Optimization
│   │   │   └── IWO/
│   │   │       ├── iwo.m                # Main IWO algorithm
│   │   │       ├── Sphere.m             # Cost function
│   │   │       └── config_iwo.m         # Configuration
│   │   ├── ga/                          # Genetic Algorithm
│   │   ├── abc/                         # Artificial Bee Colony
│   │   ├── bba/                         # Bees Behavior Algorithm
│   │   ├── sa/                          # Simulated Annealing
│   │   └── pso_example.m                # PSO example
│   │
│   ├── models/                          # Helicopter dynamics models
│   ├── utils/                           # Utility functions
│   │   └── PopulCheck.m                 # Results visualization
│   └── visualization/                   # 3D visualization & animation
│       └── Animations/                  # Animation framework
│
├── data/                                # All data files
│   ├── flight_data/                     # Raw flight test data
│   │   └── raw/                         # .bin and .mat files
│   ├── experiments/                     # Experimental results
│   │   ├── best.mat                     # Primary dataset
│   │   ├── best2.mat                    # Alternative dataset
│   │   └── *.mat                        # Other results
│   └── logs/                            # Flight and simulation logs
│       ├── flight_logs/                 # 17 flight sessions
│       └── reference/                   # Reference logs
│
├── docs/                                # Documentation & papers
│   ├── IFAC_heli_weed.doc              # Academic paper
│   └── figures/                         # Plots and diagrams
│
├── examples/                            # Example scripts (planned)
├── tests/                               # Test files (planned)
├── results/                             # Output directory for new results
│
├── README.md                            # This file
├── CONTRIBUTING.md                      # Contribution guidelines
├── LICENSE                              # Project license
├── CHANGELOG.md                         # Version history
└── CITATION.cff                         # Citation information
```

**Note**: Each directory contains its own README.md with detailed information.

---

## Results Visualization

### Using PopulCheck.m

After optimization completes, visualize the identified parameters:

```matlab
% Copy best parameters to popul variable
popul = BestSol.Position';

% Run visualization
PopulCheck
```

This generates:
- Time-domain comparison plots (actual vs. simulated)
- State trajectories for all 10 outputs
- Correlation coefficients for each output
- Model fit quality metrics

### Accessing Raw Results

```matlab
% Best solution structure
BestSol.Position  % 40 optimized parameters
BestSol.Cost      % Final cost value

% Convergence history
figure;
semilogy(BestCosts);
xlabel('Iteration');
ylabel('Best Cost');
title('Optimization Convergence');
```

---

## Troubleshooting

### Common Issues

#### 1. NaN Values During Initialization

**Problem**: Random initialization sometimes produces NaN cost values.

**Solution**: The code uses pre-defined initial populations in `best.mat` to avoid this issue. If you modify initialization:

```matlab
% Check for NaN before evaluation
if any(isnan(pop(i).Position))
    pop(i).Position = unifrnd(VarMin, VarMax);
end
```

#### 2. Missing Toolbox Error

**Problem**: `Undefined function 'ss' for input arguments of type 'double'`

**Solution**: Install Control System Toolbox:
```matlab
% Check installation
ver('control')

% Or use MATLAB Add-On Explorer to install
```

#### 3. File Not Found Errors

**Problem**: `Unable to read file 'best.mat'`

**Solution**: Ensure you're in the correct directory:
```matlab
cd 'YPEA119 Invasive Weed Optimization/IWO'
ls  % Verify best.mat exists
```

#### 4. Out of Memory Errors

**Problem**: MATLAB runs out of memory during optimization

**Solution**: Reduce population size or iterations:
```matlab
% In iwo.m, modify:
MaxIt = 1000;   % Instead of 5000
nPop = 20;      % Instead of 40
```

#### 5. Slow Performance

**Problem**: Optimization takes too long

**Optimization tips**:
- Close unnecessary MATLAB figures
- Use MATLAB's parallel computing toolbox
- Reduce iteration count for testing
- Profile the cost function: `profile on; iwo; profile viewer`

---

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Areas for Contribution

- Additional optimization algorithms
- Performance improvements
- Better visualization tools
- Documentation improvements
- Test coverage
- Bug fixes

---

## Citation

If you use this code in your research, please cite:

```bibtex
@software{krishnan2019sysid,
  author = {Krishnan, Navaneeth and Kumaar, Saumya},
  title = {System Identification for Unmanned Helicopter using Metaheuristic Algorithms},
  year = {2019},
  publisher = {GitHub},
  url = {https://github.com/dronefreak/system-identification-helicopter}
}
```

**Related Publications:**
- Bernard Mettler's helicopter dynamics papers
- IFAC paper (see `IFAC_heli_weed.doc`)

Also see [CITATION.cff](CITATION.cff) for structured citation data.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Third-Party Licenses

- **Yarpiz YPEA Framework**: BSD License (see individual algorithm directories)
- **STLRead**: BSD License
- **MatlabPlaneGraphics**: See respective directories for licenses

---

## Authors

**Creator**: [Navaneeth Krishnan](https://github.com/Navaneeth-krishnan)
**Maintainer**: [Saumya Kumaar](https://dronefreak.bitbucket.io/)

### Acknowledgments

- Bernard Mettler for the helicopter state-space model framework
- Yarpiz team for optimization algorithm implementations
- Contributors to STLRead and MatlabPlaneGraphics

---

## Support

For questions, issues, or feature requests:
- **Issues**: [GitHub Issues](https://github.com/dronefreak/system-identification-helicopter/issues)
- **Email**: Contact the maintainer via GitHub profile

---

**Last Updated**: 2026-01-06
**Version**: 2.0.0
