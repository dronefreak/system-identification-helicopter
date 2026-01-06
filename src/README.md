# Source Code Directory

This directory contains all source code for the helicopter system identification project.

## Structure

```
src/
├── algorithms/      - Optimization algorithm implementations
│   ├── iwo/        - Invasive Weed Optimization
│   ├── ga/         - Genetic Algorithm
│   ├── abc/        - Artificial Bee Colony
│   ├── bba/        - Bees Behavior Algorithm
│   ├── sa/         - Simulated Annealing
│   └── pso_example.m - Particle Swarm Optimization example
│
├── models/         - Helicopter dynamics models and cost functions
├── utils/          - Utility functions and helper scripts
└── visualization/  - Visualization and animation tools
```

## Algorithms

Each algorithm directory contains:
- Main algorithm implementation
- Cost/fitness functions adapted for helicopter system identification
- Configuration files
- Algorithm-specific utilities
- License and documentation

## Usage

To use an algorithm:

```matlab
% Navigate to the algorithm directory
cd src/algorithms/iwo/IWO

% Load flight data
load('../../../../data/experiments/best.mat')

% Run optimization
iwo
```

## Adding New Algorithms

When adding a new optimization algorithm:

1. Create a new directory under `src/algorithms/`
2. Follow the naming convention (lowercase, descriptive)
3. Include a README explaining the algorithm
4. Use the standard cost function interface
5. Add configuration file for parameters
6. Update this README with the new algorithm

## License

See individual algorithm directories for third-party licenses.
Main project code is under MIT License (see root LICENSE file).
