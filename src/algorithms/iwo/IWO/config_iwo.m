function config = config_iwo()
% CONFIG_IWO - Configuration parameters for Invasive Weed Optimization
%
% Returns a structure containing all configuration parameters for the IWO
% algorithm applied to helicopter system identification.
%
% Outputs:
%    config - Structure containing all configuration parameters
%
% Author: System Identification Project
% Date: 2026-01-06

    %% Problem Parameters
    config.nVar = 40;                    % Number of decision variables (parameters to optimize)
    config.nOutputs = 10;                % Number of measurable outputs
    config.nStates = 13;                 % Number of states in state-space model
    config.nInputs = 4;                  % Number of control inputs

    %% Physical Constants
    config.gravity = 9.81;               % Gravitational acceleration (m/s^2)

    %% IWO Algorithm Parameters
    config.maxIterations = 5000;         % Maximum number of iterations
    config.initialPopSize = 20;          % Initial population size
    config.maxPopSize = 40;              % Maximum population size
    config.minSeeds = 2;                 % Minimum number of seeds per plant
    config.maxSeeds = 7;                 % Maximum number of seeds per plant
    config.varianceExponent = 0.5;       % Variance reduction exponent
    config.sigmaInitial = 0.9;           % Initial standard deviation
    config.sigmaFinal = 0.001;           % Final standard deviation

    %% Parameter Bounds (Mettler's range)
    % Lower bounds for the 40 parameters
    config.varMin = [-1   -60   -1    -60   -1 ...
                     -1    120   0    -0.1  -0.1 ...
                      40   0     0.01 -1     0 ...
                     -1   -1    -20   -160  -1 ...
                     -2    0    -10    0    -10 ...
                    -100 -5    -20    0.01  0 ...
                      0   -0.1  -1    -1    -0.1 ...
                    -100 -80   -10   -1    -1];

    % Upper bounds for the 40 parameters
    config.varMax = [ 1     60    1    60    0 ...
                      1     220   0    0.1   0.1 ...
                      120   0     1    0     1 ...
                      1     1     20   100   1 ...
                      2     0.01  10   0     10 ...
                      100   5     20   1     0 ...
                      0     0.1   1    1     0.1 ...
                      100   80    10   1     1];

    %% Cost Function Parameters
    config.maxCostPenalty = Inf;         % Penalty for invalid parameters (NaN)
    config.targetOutputs = 8;            % Target number of outputs for correlation sum
    config.correlationTarget = 8;        % Perfect correlation sum would be 8 (for 8 outputs)

    %% Display Options
    config.displayInterval = 1;          % Display every N iterations (1 = every iteration)
    config.plotResults = true;           % Generate convergence plot
    config.useSemilogy = true;           % Use logarithmic scale for cost plot

    %% Performance Options
    config.useParallel = false;          % Enable parallel cost function evaluation
                                          % Requires Parallel Computing Toolbox
    config.showProgress = true;          % Show progress bar during optimization
    config.progressUpdateInterval = 10;  % Update progress every N iterations

    %% Memory Optimization Options
    config.preallocateOffspring = true;  % Preallocate offspring arrays
    config.clearIntermediateData = true; % Clear intermediate data to reduce memory

    %% Reproducibility Options
    config.useRandomSeed = false;        % Enable fixed random seed for reproducibility
    config.randomSeed = 42;              % Random seed value (used if useRandomSeed = true)
                                          % Set to specific value for reproducible results
                                          % Set to 'shuffle' to use current time
    config.saveRandomState = true;       % Save random number generator state with results

end
