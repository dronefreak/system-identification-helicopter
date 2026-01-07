function fig = plot_comparison(actualData, simulatedData, timeVector, varargin)
% PLOT_COMPARISON - Plot actual vs simulated helicopter response
%
% This function creates standardized comparison plots between actual flight
% data and simulated model outputs for system identification validation.
%
% Usage:
%    plot_comparison(actual, simulated, time)
%    plot_comparison(actual, simulated, time, 'outputs', [1 2 3])
%    fig = plot_comparison(actual, simulated, time, 'saveas', 'comparison.png')
%
% Inputs:
%    actualData    - Actual flight data (10 × N matrix)
%    simulatedData - Simulated model output (10 × N or 13 × N matrix)
%    timeVector    - Time vector (N × 1 or 1 × N)
%
% Optional Parameters (Name-Value pairs):
%    'outputs'      - Output indices to plot (default: 1:10)
%    'outputNames'  - Cell array of output names
%    'title'        - Main title (default: 'Actual vs Simulated Response')
%    'subplot'      - Layout [rows, cols] (default: auto)
%    'linewidth'    - Line width (default: 1.5)
%    'saveas'       - Save to file
%    'showCorr'     - Show correlation coefficients (default: true)
%
% Outputs:
%    fig            - Figure handle
%
% Example:
%    % Plot first 4 outputs
%    plot_comparison(out_H, simulated, time, ...
%                    'outputs', [1 2 3 4], ...
%                    'outputNames', {'u', 'v', 'w', 'p'});
%
% See also: plot_convergence, PopulCheck, generate_report
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Inputs
    p = inputParser;
    addRequired(p, 'actualData', @isnumeric);
    addRequired(p, 'simulatedData', @isnumeric);
    addRequired(p, 'timeVector', @isvector);
    addParameter(p, 'outputs', 1:10, @isnumeric);
    addParameter(p, 'outputNames', {}, @iscell);
    addParameter(p, 'title', 'Actual vs Simulated Response', @ischar);
    addParameter(p, 'subplot', [], @isnumeric);
    addParameter(p, 'linewidth', 1.5, @isnumeric);
    addParameter(p, 'saveas', '', @ischar);
    addParameter(p, 'showCorr', true, @islogical);
    parse(p, actualData, simulatedData, timeVector, varargin{:});

    actualData = p.Results.actualData;
    simulatedData = p.Results.simulatedData;
    timeVector = p.Results.timeVector;
    outputIndices = p.Results.outputs;
    outputNames = p.Results.outputNames;
    plotTitle = p.Results.title;
    subplotLayout = p.Results.subplot;
    lineWidth = p.Results.linewidth;
    saveFile = p.Results.saveas;
    showCorr = p.Results.showCorr;

    %% Validate Data
    timeVector = timeVector(:);
    nSamples = length(timeVector);

    if size(actualData, 2) ~= nSamples
        actualData = actualData';
    end

    if size(simulatedData, 2) ~= nSamples
        simulatedData = simulatedData';
    end

    % Extract relevant outputs from simulated data
    if size(simulatedData, 1) == 13
        simulatedData = simulatedData(1:10, :);  % Use first 10 states
    end

    nOutputs = length(outputIndices);

    % Default output names
    if isempty(outputNames)
        defaultNames = {'u (m/s)', 'v (m/s)', 'w (m/s)', ...
                       'p (rad/s)', 'q (rad/s)', 'r (rad/s)', ...
                       'φ (rad)', 'θ (rad)', 'ψ (rad)', 'a (rad)'};
        outputNames = defaultNames(outputIndices);
    end

    % Auto subplot layout
    if isempty(subplotLayout)
        if nOutputs <= 4
            subplotLayout = [2, 2];
        elseif nOutputs <= 6
            subplotLayout = [2, 3];
        elseif nOutputs <= 9
            subplotLayout = [3, 3];
        else
            subplotLayout = [4, 3];
        end
    end

    %% Create Figure
    fig = figure('Name', plotTitle, 'NumberTitle', 'off', ...
                 'Position', [50 50 1200 800]);

    %% Plot Each Output
    for i = 1:nOutputs
        idx = outputIndices(i);

        subplot(subplotLayout(1), subplotLayout(2), i);
        hold on;

        % Plot actual data
        plot(timeVector, actualData(idx, :), 'b-', ...
             'LineWidth', lineWidth, 'DisplayName', 'Actual');

        % Plot simulated data
        plot(timeVector, simulatedData(idx, :), 'r--', ...
             'LineWidth', lineWidth, 'DisplayName', 'Simulated');

        hold off;

        % Labels and title
        xlabel('Time (s)', 'FontSize', 9);
        ylabel(outputNames{i}, 'FontSize', 9);

        % Add correlation coefficient if requested
        if showCorr
            try
                corrMatrix = corrcoef(actualData(idx, :), simulatedData(idx, :));
                corrCoef = abs(corrMatrix(1, 2));
                titleStr = sprintf('%s (ρ=%.3f)', outputNames{i}, corrCoef);
            catch
                titleStr = outputNames{i};
            end
        else
            titleStr = outputNames{i};
        end

        title(titleStr, 'FontSize', 10, 'FontWeight', 'bold');

        % Format
        grid on;
        legend('Location', 'best', 'FontSize', 8);
        ax = gca;
        ax.FontSize = 8;
    end

    % Super title
    sgtitle(plotTitle, 'FontSize', 14, 'FontWeight', 'bold');

    %% Save if requested
    if ~isempty(saveFile)
        [~, ~, ext] = fileparts(saveFile);
        if isempty(ext)
            saveFile = [saveFile '.png'];
        end

        try
            saveas(fig, saveFile);
            fprintf('✓ Comparison plot saved to: %s\n', saveFile);
        catch ME
            warning('PlotComparison:SaveError', ...
                    'Failed to save plot: %s', ME.message);
        end
    end
end
