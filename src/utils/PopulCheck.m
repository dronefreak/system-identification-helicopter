function PopulCheck(varargin)
% POPULCHECK - Visualize optimized helicopter model vs. actual flight data
%
% This function creates comparison plots between simulated helicopter
% dynamics (using optimized parameters) and actual flight test data.
%
% Usage:
%    PopulCheck                 - Uses workspace variables
%    PopulCheck(params)         - Uses provided parameter vector
%    PopulCheck(params, data)   - Uses provided parameters and data structure
%
% Required workspace variables (if no arguments provided):
%    popul  - 40×1 vector of optimized parameters
%    inr    - N×4 matrix of control inputs
%    outr   - N×10 matrix of measured outputs
%    t      - N×1 time vector
%    cpop   - 13×1 vector of output scaling factors (optional, defaults to ones)
%
% Inputs (optional):
%    params - 40×1 vector of state-space parameters
%    data   - Structure with fields: inr, outr, t, cpop (optional)
%
% Outputs:
%    Creates a 4×3 subplot figure showing:
%    - Control inputs (4 plots)
%    - State responses (6 plots)
%    - Comparison between simulated (dashed) and actual (solid) data
%
% Author: System Identification Project
% Date: 2026-01-06

    %% Constants
    GRAVITY = 9.81;  % m/s^2
    NUM_PARAMETERS = 40;
    NUM_STATES = 13;
    NUM_INPUTS = 4;
    NUM_OUTPUTS = 10;

    %% Parse Inputs
    if nargin == 0
        % Load from workspace
        try
            params = evalin('base', 'popul');
            inputData = evalin('base', 'inr');
            outputData = evalin('base', 'outr');
            timeVector = evalin('base', 't');

            % Check if cpop exists, otherwise use ones
            if evalin('base', 'exist(''cpop'', ''var'')')
                outputScaling = evalin('base', 'cpop');
            else
                outputScaling = ones(NUM_STATES, 1);
            end
        catch ME
            error('PopulCheck:MissingVariables', ...
                  ['Required workspace variables not found.\n' ...
                   'Please ensure popul, inr, outr, and t are loaded.\n' ...
                   'Original error: %s'], ME.message);
        end

    elseif nargin == 1
        params = varargin{1};
        % Load data from workspace
        try
            inputData = evalin('base', 'inr');
            outputData = evalin('base', 'outr');
            timeVector = evalin('base', 't');
            if evalin('base', 'exist(''cpop'', ''var'')')
                outputScaling = evalin('base', 'cpop');
            else
                outputScaling = ones(NUM_STATES, 1);
            end
        catch ME
            error('PopulCheck:MissingData', ...
                  'Could not load data from workspace: %s', ME.message);
        end

    elseif nargin == 2
        params = varargin{1};
        data = varargin{2};

        % Extract from data structure
        if ~isstruct(data)
            error('PopulCheck:InvalidDataType', ...
                  'Second argument must be a structure with fields: inr, outr, t');
        end

        inputData = data.inr;
        outputData = data.outr;
        timeVector = data.t;

        if isfield(data, 'cpop')
            outputScaling = data.cpop;
        else
            outputScaling = ones(NUM_STATES, 1);
        end

    else
        error('PopulCheck:TooManyInputs', ...
              'Too many input arguments. Use 0, 1, or 2 arguments.');
    end

    %% Validate Inputs
    if ~isvector(params) || length(params) ~= NUM_PARAMETERS
        error('PopulCheck:InvalidParameters', ...
              'Parameters must be a %d-element vector', NUM_PARAMETERS);
    end

    if size(inputData, 2) ~= NUM_INPUTS
        error('PopulCheck:InvalidInputData', ...
              'Input data must have %d columns', NUM_INPUTS);
    end

    if size(outputData, 2) ~= NUM_OUTPUTS
        error('PopulCheck:InvalidOutputData', ...
              'Output data must have %d columns', NUM_OUTPUTS);
    end

    if length(timeVector) ~= size(inputData, 1)
        error('PopulCheck:DimensionMismatch', ...
              'Time vector length must match input data rows');
    end

    %% Extract Parameters
    % Ensure column vector
    p = params(:);
    j = 1;  % Column index for compatibility with original code

    % State-space parameters
    Xu = p(1);
    Xa = p(2);
    Yv = p(3);
    Yb = p(4);
    Lu = p(5);
    Lv = p(6);
    Lb = p(7);
    Lw = p(8);
    Mu = p(9);
    Mv = p(10);
    Ma = p(11);
    Mw = p(12);
    Tf = p(13);
    Ab = p(14);
    Ac = p(15);
    Ba = p(16);
    Bd = p(17);
    Za = p(18);
    Zb = p(19);
    Zw = p(20);
    Zr = p(21);
    Zrfb = 0;  % Fixed parameter
    Nv = p(22);
    Np = p(23);
    Nw = p(24);
    Nr = p(25);
    Nrfb = p(26);
    Kr = p(27);
    Krfb = p(28);
    Ts = p(29);
    Yped = p(30);
    Mcol = p(31);
    Alat = p(32);
    Alon = p(33);
    Blat = p(34);
    Blon = p(35);
    Zcol = p(36);
    Nped = p(37);
    Ncol = p(38);
    Clon = p(39);
    Dlat = p(40);

    %% Construct State-Space Model
    A = [Xu, 0, 0, 0, 0, -GRAVITY, Xa, 0, 0, 0, 0, 0, 0;
         0, Yv, 0, 0, GRAVITY, 0, 0, Yb, 0, 0, 0, 0, 0;
         Lu, Lv, 0, 0, 0, 0, 0, Lb, Lw, 0, 0, 0, 0;
         Mu, Mv, 0, 0, 0, 0, Ma, 0, Mw, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, -Tf, 0, 0, -1, Ab, 0, 0, 0, Ac, 0;
         0, 0, -Tf, 0, 0, 0, Ba, -1, 0, 0, 0, 0, Bd;
         0, 0, 0, 0, 0, 0, Za, Zb, Zw, Zr, Zrfb, 0, 0;
         0, Nv, Np, 0, 0, 0, 0, 0, Nw, Nr, Nrfb, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0, 0, Kr, Krfb, 0, 0;
         0, 0, 0, -Ts, 0, 0, 0, 0, 0, 0, 0, -1, 0;
         0, 0, -Ts, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1];

    B = [0, 0, 0, 0;
         0, 0, Yped, 0;
         0, 0, 0, 0;
         0, 0, 0, Mcol;
         0, 0, 0, 0;
         0, 0, 0, 0;
         Alat, Alon, 0, 0;
         Blat, Blon, 0, 0;
         0, 0, 0, Zcol;
         0, 0, Nped, Ncol;
         0, 0, 0, 0;
         0, Clon, 0, 0;
         Dlat, 0, 0, 0];

    % Output scaling matrix
    C = diag(outputScaling(:));

    D = 0;

    %% Simulate Model
    try
        model = ss(A, B, C, D);
        simulatedOutput = lsim(model, inputData, timeVector);
    catch ME
        error('PopulCheck:SimulationFailed', ...
              'Model simulation failed: %s', ME.message);
    end

    %% Create Visualization
    figure('Name', 'Parameter Validation: Model vs. Flight Data', ...
           'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    % Subplot 1: Pitch input
    subplot(4, 3, 1);
    plot(timeVector, inputData(:, 2), 'LineWidth', 1.5);
    title('Pitch Input');
    xlabel('Time (s)');
    ylabel('Input');
    grid on;

    % Subplot 2: Roll input
    subplot(4, 3, 2);
    plot(timeVector, inputData(:, 1), 'LineWidth', 1.5);
    title('Roll Input');
    xlabel('Time (s)');
    ylabel('Input');
    grid on;

    % Subplot 3: Yaw input
    subplot(4, 3, 3);
    plot(timeVector, inputData(:, 3), 'LineWidth', 1.5);
    title('Yaw Input');
    xlabel('Time (s)');
    ylabel('Input');
    grid on;

    % Subplot 4: Roll rate
    subplot(4, 3, 4);
    plot(timeVector, simulatedOutput(:, 3), '-.', ...
         timeVector, outputData(:, 3), '-', 'LineWidth', 1.5);
    ylim([-2, 2]);
    title('Roll Rate');
    xlabel('Time (s)');
    ylabel('Rate (rad/s)');
    legend('Simulated', 'Actual', 'Location', 'best');
    grid on;

    % Subplot 5: Pitch rate
    subplot(4, 3, 5);
    plot(timeVector, simulatedOutput(:, 4), '-.', ...
         timeVector, outputData(:, 4), '-', 'LineWidth', 1.5);
    title('Pitch Rate');
    xlabel('Time (s)');
    ylabel('Rate (rad/s)');
    legend('Simulated', 'Actual', 'Location', 'best');
    grid on;

    % Subplot 6: Yaw rate
    subplot(4, 3, 6);
    plot(timeVector, simulatedOutput(:, 10), '-.', ...
         timeVector, outputData(:, 10), '-', 'LineWidth', 1.5);
    ylim([-5, 5]);
    title('Yaw Rate');
    xlabel('Time (s)');
    ylabel('Rate (rad/s)');
    legend('Simulated', 'Actual', 'Location', 'best');
    grid on;

    % Subplot 7: Pitch angle
    subplot(4, 3, 7);
    plot(timeVector, simulatedOutput(:, 6), '-.', ...
         timeVector, outputData(:, 6), '-', 'LineWidth', 1.5);
    ylim([-1, 1]);
    title('Pitch Angle');
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('Simulated', 'Actual', 'Location', 'best');
    grid on;

    % Subplot 8: Roll angle
    subplot(4, 3, 8);
    plot(timeVector, simulatedOutput(:, 5), '-.', ...
         timeVector, outputData(:, 5), '-', 'LineWidth', 1.5);
    ylim([-2, 2]);
    title('Roll Angle');
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('Simulated', 'Actual', 'Location', 'best');
    grid on;

    % Subplot 9: Collective input
    subplot(4, 3, 9);
    plot(timeVector, inputData(:, 4), 'LineWidth', 1.5);
    title('Collective Input');
    xlabel('Time (s)');
    ylabel('Input');
    grid on;

    % Subplot 10: u velocity
    subplot(4, 3, 10);
    plot(timeVector, simulatedOutput(:, 1), '-.', ...
         timeVector, outputData(:, 1), '-', 'LineWidth', 1.5);
    ylim([-20, 20]);
    title('u Velocity');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    legend('Simulated', 'Actual', 'Location', 'best');
    grid on;

    % Subplot 11: v velocity
    subplot(4, 3, 11);
    plot(timeVector, simulatedOutput(:, 2), '-.', ...
         timeVector, outputData(:, 2), '-', 'LineWidth', 1.5);
    ylim([-10, 10]);
    title('v Velocity');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    legend('Simulated', 'Actual', 'Location', 'best');
    grid on;

    % Subplot 12: w velocity
    subplot(4, 3, 12);
    plot(timeVector, simulatedOutput(:, 9), '-.', ...
         timeVector, outputData(:, 9), '-', 'LineWidth', 1.5);
    ylim([-10, 10]);
    title('w Velocity');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    legend('Simulated', 'Actual', 'Location', 'best');
    grid on;

    %% Print Summary Statistics
    fprintf('\n=== Model Validation Summary ===\n');
    fprintf('Parameters: %d optimized values\n', NUM_PARAMETERS);
    fprintf('Time span: %.2f to %.2f seconds (%d data points)\n', ...
            timeVector(1), timeVector(end), length(timeVector));

    % Calculate and display correlation coefficients
    fprintf('\nCorrelation Coefficients:\n');
    stateNames = {'u', 'v', 'p', 'q', 'φ', 'θ', 'a', 'b', 'w', 'r'};

    for i = 1:min(NUM_OUTPUTS, size(simulatedOutput, 2))
        corrMatrix = corrcoef(simulatedOutput(:, i), outputData(:, i));
        corrCoeff = corrMatrix(1, 2);
        fprintf('  %6s: %.4f\n', stateNames{i}, corrCoeff);
    end

    fprintf('\n');

end
