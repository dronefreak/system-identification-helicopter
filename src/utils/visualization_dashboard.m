function fig = visualization_dashboard(results, varargin)
% VISUALIZATION_DASHBOARD - Create comprehensive optimization results dashboard
%
% This function creates a multi-panel dashboard displaying key metrics and
% visualizations of IWO optimization results in a single figure.
%
% Usage:
%    visualization_dashboard(results)
%    visualization_dashboard(results, 'layout', '2x3')
%    fig = visualization_dashboard(results, 'saveas', 'dashboard.png')
%
% Inputs:
%    results    - Results structure from IWO or run_experiment
%
% Optional Parameters (Name-Value pairs):
%    'layout'     - Dashboard layout: '2x2', '2x3', '3x3' (default: '2x3')
%    'panels'     - Cell array of panels to show (default: all available)
%                   Options: 'convergence', 'cost_stats', 'improvement',
%                           'parameter_dist', 'timing', 'metadata'
%    'title'      - Dashboard title (default: 'IWO Optimization Dashboard')
%    'saveas'     - Save dashboard to file (.png, .pdf, .fig) (default: '')
%    'fontSize'   - Base font size (default: 10)
%    'colorScheme'- Color scheme: 'default', 'colorful', 'grayscale' (default: 'default')
%
% Available Panels:
%    - convergence: Convergence plot over iterations
%    - cost_stats: Statistical summary of costs
%    - improvement: Improvement percentage over iterations
%    - parameter_dist: Parameter value distribution
%    - timing: Execution time breakdown
%    - metadata: Experiment metadata and configuration
%
% Outputs:
%    fig          - Figure handle
%
% Examples:
%    % Default dashboard
%    visualization_dashboard(results);
%
%    % Custom layout and save
%    visualization_dashboard(results, 'layout', '3x3', 'saveas', 'dashboard.png');
%
%    % Specific panels only
%    visualization_dashboard(results, 'panels', {'convergence', 'cost_stats'});
%
% See also: plot_convergence, generate_report, animate_convergence
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Inputs
    p = inputParser;
    addRequired(p, 'results', @isstruct);
    addParameter(p, 'layout', '2x3', @(x) any(validatestring(x, {'2x2', '2x3', '3x3'})));
    addParameter(p, 'panels', {}, @iscell);
    addParameter(p, 'title', 'IWO Optimization Dashboard', @ischar);
    addParameter(p, 'saveas', '', @ischar);
    addParameter(p, 'fontSize', 10, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'colorScheme', 'default', @(x) any(validatestring(x, {'default', 'colorful', 'grayscale'})));
    parse(p, results, varargin{:});

    layout = p.Results.layout;
    requestedPanels = p.Results.panels;
    dashboardTitle = p.Results.title;
    saveFile = p.Results.saveas;
    fontSize = p.Results.fontSize;
    colorScheme = p.Results.colorScheme;

    %% Extract Data
    if isfield(results, 'bestCosts')
        bestCosts = results.bestCosts;
    elseif isfield(results, 'BestCosts')
        bestCosts = results.BestCosts;
    else
        error('VisualizationDashboard:MissingData', 'Results must contain bestCosts field');
    end

    if isfield(results, 'bestSolution')
        bestSolution = results.bestSolution;
    elseif isfield(results, 'BestSol')
        bestSolution = results.BestSol;
    else
        bestSolution = struct();
    end

    %% Determine Available Panels
    availablePanels = {'convergence', 'cost_stats', 'improvement'};

    if isfield(bestSolution, 'Position') && ~isempty(bestSolution.Position)
        availablePanels{end+1} = 'parameter_dist';
    end

    if isfield(results, 'metadata') && isfield(results.metadata, 'executionTime')
        availablePanels{end+1} = 'timing';
    end

    if isfield(results, 'metadata') || isfield(results, 'config')
        availablePanels{end+1} = 'metadata';
    end

    % Use requested panels or all available
    if isempty(requestedPanels)
        panels = availablePanels;
    else
        panels = requestedPanels;
    end

    %% Setup Layout
    switch layout
        case '2x2'
            nRows = 2; nCols = 2;
        case '2x3'
            nRows = 2; nCols = 3;
        case '3x3'
            nRows = 3; nCols = 3;
    end

    maxPanels = nRows * nCols;
    if length(panels) > maxPanels
        warning('VisualizationDashboard:TooManyPanels', ...
                'Requested %d panels but layout only supports %d. Using first %d panels.', ...
                length(panels), maxPanels, maxPanels);
        panels = panels(1:maxPanels);
    end

    %% Color Scheme
    switch colorScheme
        case 'default'
            colorPrimary = [0, 0.4470, 0.7410];
            colorSecondary = [0.8500, 0.3250, 0.0980];
            colorTertiary = [0.9290, 0.6940, 0.1250];
        case 'colorful'
            colorPrimary = [0.2, 0.6, 0.8];
            colorSecondary = [0.9, 0.2, 0.4];
            colorTertiary = [0.4, 0.8, 0.3];
        case 'grayscale'
            colorPrimary = [0.2, 0.2, 0.2];
            colorSecondary = [0.5, 0.5, 0.5];
            colorTertiary = [0.7, 0.7, 0.7];
    end

    %% Create Figure
    fig = figure('Position', [50, 50, 1200, 800], 'Color', 'white');

    % Main title
    annotation('textbox', [0, 0.96, 1, 0.04], 'String', dashboardTitle, ...
               'FontSize', fontSize + 4, 'FontWeight', 'bold', ...
               'HorizontalAlignment', 'center', 'EdgeColor', 'none');

    %% Create Each Panel
    for i = 1:length(panels)
        subplot(nRows, nCols, i);

        switch panels{i}
            case 'convergence'
                plot_panel_convergence(bestCosts, colorPrimary, fontSize);

            case 'cost_stats'
                plot_panel_cost_stats(bestCosts, colorPrimary, colorSecondary, fontSize);

            case 'improvement'
                plot_panel_improvement(bestCosts, colorPrimary, fontSize);

            case 'parameter_dist'
                plot_panel_parameter_dist(bestSolution.Position, colorPrimary, fontSize);

            case 'timing'
                plot_panel_timing(results, colorPrimary, colorSecondary, colorTertiary, fontSize);

            case 'metadata'
                plot_panel_metadata(results, fontSize);
        end
    end

    %% Save Figure
    if ~isempty(saveFile)
        fprintf('Saving dashboard to: %s\n', saveFile);
        [~, ~, ext] = fileparts(saveFile);

        switch ext
            case '.png'
                print(fig, saveFile, '-dpng', '-r300');
            case '.pdf'
                print(fig, saveFile, '-dpdf', '-r300');
            case '.fig'
                savefig(fig, saveFile);
            otherwise
                warning('VisualizationDashboard:UnknownFormat', ...
                        'Unknown file format: %s. Supported: .png, .pdf, .fig', ext);
        end
    end

    fprintf('✓ Dashboard created successfully!\n');
end

%% Panel Functions

function plot_panel_convergence(bestCosts, color, fontSize)
    semilogy(1:length(bestCosts), bestCosts, 'Color', color, 'LineWidth', 2);
    xlabel('Iteration', 'FontSize', fontSize);
    ylabel('Best Cost (log)', 'FontSize', fontSize);
    title('Convergence History', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
    grid on;
    xlim([0, length(bestCosts)]);
end

function plot_panel_cost_stats(bestCosts, color1, color2, fontSize)
    stats = [bestCosts(1), bestCosts(end), min(bestCosts), mean(bestCosts)];
    labels = {'Initial', 'Final', 'Best', 'Mean'};

    bar(stats, 'FaceColor', color1, 'EdgeColor', color2, 'LineWidth', 1.5);
    set(gca, 'XTickLabel', labels, 'FontSize', fontSize);
    ylabel('Cost', 'FontSize', fontSize);
    title('Cost Statistics', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
    grid on;

    % Add value labels on bars
    for i = 1:length(stats)
        text(i, stats(i), sprintf('%.2f', stats(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', fontSize - 2);
    end
end

function plot_panel_improvement(bestCosts, color, fontSize)
    improvement = (bestCosts(1) - bestCosts) ./ bestCosts(1) * 100;
    plot(1:length(improvement), improvement, 'Color', color, 'LineWidth', 2);
    xlabel('Iteration', 'FontSize', fontSize);
    ylabel('Improvement (%)', 'FontSize', fontSize);
    title('Cost Improvement Over Time', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
    grid on;
    xlim([0, length(improvement)]);
    ylim([0, max(improvement) * 1.1]);
end

function plot_panel_parameter_dist(parameters, color, fontSize)
    histogram(parameters, 20, 'FaceColor', color, 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    xlabel('Parameter Value', 'FontSize', fontSize);
    ylabel('Frequency', 'FontSize', fontSize);
    title(sprintf('Parameter Distribution (n=%d)', length(parameters)), ...
          'FontSize', fontSize + 2, 'FontWeight', 'bold');
    grid on;
end

function plot_panel_timing(results, color1, color2, color3, fontSize)
    if isfield(results, 'timing') && isstruct(results.timing)
        timing = results.timing;
        fields = fieldnames(timing);
        values = zeros(length(fields), 1);

        for i = 1:length(fields)
            values(i) = timing.(fields{i});
        end

        colors = [color1; color2; color3];
        if length(values) > 3
            colors = repmat(colors, ceil(length(values)/3), 1);
        end
        colors = colors(1:length(values), :);

        pie(values, fields);
        colormap(colors);
        title('Execution Time Breakdown', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
    elseif isfield(results, 'metadata') && isfield(results.metadata, 'executionTime')
        text(0.5, 0.5, sprintf('Total Time: %.2f sec', results.metadata.executionTime), ...
             'HorizontalAlignment', 'center', 'FontSize', fontSize);
        title('Execution Time', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
        axis off;
    end
end

function plot_panel_metadata(results, fontSize)
    axis off;
    yPos = 0.9;
    lineHeight = 0.12;

    text(0.05, yPos, '\bfMetadata', 'FontSize', fontSize + 2);
    yPos = yPos - lineHeight;

    if isfield(results, 'metadata')
        meta = results.metadata;

        if isfield(meta, 'experimentName')
            text(0.05, yPos, sprintf('Experiment: %s', meta.experimentName), 'FontSize', fontSize);
            yPos = yPos - lineHeight;
        end

        if isfield(meta, 'finalCost')
            text(0.05, yPos, sprintf('Final Cost: %.6f', meta.finalCost), 'FontSize', fontSize);
            yPos = yPos - lineHeight;
        end

        if isfield(meta, 'improvement')
            text(0.05, yPos, sprintf('Improvement: %.2f%%', meta.improvement), 'FontSize', fontSize);
            yPos = yPos - lineHeight;
        end

        if isfield(meta, 'executionTime')
            text(0.05, yPos, sprintf('Time: %.2f min', meta.executionTime / 60), 'FontSize', fontSize);
            yPos = yPos - lineHeight;
        end

        if isfield(meta, 'platform')
            text(0.05, yPos, sprintf('Platform: %s', meta.platform), 'FontSize', fontSize - 1);
            yPos = yPos - lineHeight;
        end
    end

    if isfield(results, 'config')
        yPos = yPos - lineHeight/2;
        text(0.05, yPos, '\bfConfiguration', 'FontSize', fontSize + 2);
        yPos = yPos - lineHeight;

        config = results.config;
        if isfield(config, 'nPop')
            text(0.05, yPos, sprintf('Population: %d', config.nPop), 'FontSize', fontSize);
            yPos = yPos - lineHeight;
        end

        if isfield(config, 'MaxIt')
            text(0.05, yPos, sprintf('Iterations: %d', config.MaxIt), 'FontSize', fontSize);
        end
    end
end
