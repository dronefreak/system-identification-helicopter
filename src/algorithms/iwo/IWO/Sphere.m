function cost = Sphere(parameters, inputData, outputData, timeVector)
% SPHERE - Cost function for helicopter system identification
%
% This function evaluates the cost of a given parameter set by:
% 1. Creating a state-space model from the 40 parameters
% 2. Simulating the model with actual input data
% 3. Comparing simulated output with actual flight test data
% 4. Computing cost based on correlation coefficients
%
% Inputs:
%    parameters  - 1×40 array of state-space model parameters
%    inputData   - N×4 matrix of control inputs (lateral, longitudinal, pedal, collective)
%    outputData  - N×10 matrix of measured outputs
%    timeVector  - N×1 time vector for simulation
%
% Outputs:
%    cost        - Scalar cost value (lower is better)
%                  Returns Inf for invalid parameters
%
% Author: System Identification Project
% Date: 2026-01-06

    %% Input Validation
    if nargin < 4
        error('Sphere:NotEnoughInputs', ...
              'This function requires 4 inputs: parameters, inputData, outputData, timeVector');
    end

    % Validate parameter dimensions
    if ~isvector(parameters) || length(parameters) ~= 40
        error('Sphere:InvalidParameters', ...
              'Parameters must be a vector of length 40');
    end

    % Validate data dimensions
    if size(inputData, 2) ~= 4
        error('Sphere:InvalidInputData', ...
              'inputData must have 4 columns (control inputs)');
    end

    if size(outputData, 2) ~= 10
        error('Sphere:InvalidOutputData', ...
              'outputData must have 10 columns (measurements)');
    end

    % Check for NaN or Inf in parameters
    if any(isnan(parameters)) || any(isinf(parameters))
        cost = Inf;
        return;
    end

    %% Load Configuration
    config = config_iwo();
    gravity = config.gravity;

    %% Extract Parameters
    % Ensure parameters is a column vector
    params = parameters(:);

    % Velocity derivatives
    Xu = params(1);
    Xa = params(2);
    Yv = params(3);
    Yb = params(4);

    % Rotational derivatives
    Lu = params(5);
    Lv = params(6);
    Lb = params(7);
    Lw = params(8);
    Mu = params(9);
    Mv = params(10);
    Ma = params(11);
    Mw = params(12);

    % Time constants
    Tf = params(13);
    Ts = params(29);

    % Flapping dynamics
    Ab = params(14);
    Ac = params(15);
    Ba = params(16);
    Bd = params(17);

    % Vertical dynamics
    Za = params(18);
    Zb = params(19);
    Zw = params(20);
    Zr = params(21);

    % Yaw dynamics
    Nv = params(22);
    Np = params(23);
    Nw = params(24);
    Nr = params(25);
    Nrfb = params(26);
    Kr = params(27);
    Krfb = params(28);

    % Control effectiveness parameters
    Yped = params(30);
    Mcol = params(31);
    Alat = params(32);
    Alon = params(33);
    Blat = params(34);
    Blon = params(35);
    Zcol = params(36);
    Nped = params(37);
    Ncol = params(38);
    Clon = params(39);
    Dlat = params(40);

    %% Construct State-Space Matrices
    % A matrix (13×13) - System dynamics
    A = [Xu    0      0      0      0      -gravity  Xa        0        0      0      0        0        0
         0     Yv     0      0      gravity 0        0         Yb       0      0      0        0        0
         Lu    Lv     0      0      0       0        0         Lb       Lw     0      0        0        0
         Mu    Mv     0      0      0       0        Ma        0        Mw     0      0        0        0
         0     0      1      0      0       0        0         0        0      0      0        0        0
         0     0      0      1      0       0        0         0        0      0      0        0        0
         0     0      0     -1      0       0       -1/Tf      Ab/Tf    0      0      0        Ac/Tf    0
         0     0     -1      0      0       0        Ba/Tf    -1/Tf     0      0      0        0        Bd/Tf
         0     0      0      0      0       0        Za        Zb       Zw     Zr     0        0        0
         0     Nv     Np     0      0       0        0         0        Nw     Nr     Nrfb     0        0
         0     0      0      0      0       0        0         0        0      Kr     Krfb     0        0
         0     0      0     -1      0       0        0         0        0      0      0       -1/Ts     0
         0     0     -1      0      0       0        0         0        0      0      0        0       -1/Ts];

    % B matrix (13×4) - Control inputs
    B = [0         0         0        0
         0         0         0        0
         0         0         Yped     0
         0         0         0        Mcol
         0         0         0        0
         0         0         0        0
         Alat/Tf   Alon/Tf   0        0
         Blat/Tf   Blon/Tf   0        0
         0         0         0        Zcol
         0         0         Nped     Ncol
         0         0         0        0
         0         Clon/Ts   0        0
         Dlat/Ts   0         0        0];

    % C matrix (13×13) - Output mapping (identity for full state observation)
    C = eye(13);

    % D matrix - No direct feedthrough
    D = 0;

    %% Create and Simulate State-Space Model
    try
        % Create state-space model
        model = ss(A, B, C, D);

        % Simulate model response
        simulatedOutput = lsim(model, inputData, timeVector);

    catch ME
        % If model creation or simulation fails, return high cost
        warning('Sphere:SimulationError', ...
                'Error in model simulation: %s', ME.message);
        cost = Inf;
        return;
    end

    %% Calculate Correlation Coefficients
    % We use only the first 10 outputs (exclude states 11-13 which cannot be measured)
    correlationSum = 0;
    numOutputs = min(10, size(outputData, 2));

    for i = 1:numOutputs
        try
            % Calculate correlation coefficient matrix
            corrMatrix = corrcoef(simulatedOutput(:, i), outputData(:, i));

            % Extract correlation coefficient (off-diagonal element)
            correlationCoeff = abs(corrMatrix(1, 2));

            % Check for NaN (can occur if std dev is zero)
            if isnan(correlationCoeff)
                correlationCoeff = 0;
            end

            % Sum correlations (skip outputs 7 and 8 as per original logic)
            if i ~= 7 && i ~= 8
                correlationSum = correlationSum + correlationCoeff;
            end

        catch ME
            warning('Sphere:CorrelationError', ...
                    'Error calculating correlation for output %d: %s', i, ME.message);
            % If correlation calculation fails, assume no correlation
            correlationSum = correlationSum + 0;
        end
    end

    %% Compute Cost
    % Target correlation sum is 8 (perfect correlation for 8 outputs)
    % Cost = difference from perfect correlation
    if isreal(correlationSum) && ~isnan(correlationSum)
        cost = config.correlationTarget - correlationSum;
    else
        cost = Inf;
    end

    % Ensure cost is non-negative
    if cost < 0
        cost = 0;
    end

end
