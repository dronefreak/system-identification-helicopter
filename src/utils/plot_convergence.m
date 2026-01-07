function fig = plot_convergence(varargin)
% PLOT_CONVERGENCE - Plot optimization convergence history
%
% This function creates standardized convergence plots for optimization results
% with customizable styling and comparison support.
%
% Usage:
%    plot_convergence(bestCosts)
%    plot_convergence(bestCosts, 'title', 'IWO Convergence')
%    plot_convergence(bestCosts1, bestCosts2, 'labels', {'Run 1', 'Run 2'})
%    fig = plot_convergence(results.bestCosts, 'scale', 'semilogy')
%
% Inputs:
%    bestCosts  - Cost history vector or multiple vectors for comparison
%
% Optional Parameters (Name-Value pairs):
%    'title'      - Plot title (default: 'Optimization Convergence')
%    'xlabel'     - X-axis label (default: 'Iteration')
%    'ylabel'     - Y-axis label (default: 'Best Cost')
%    'scale'      - 'linear' or 'semilogy' (default: 'semilogy')
%    'labels'     - Cell array of labels for multiple plots
%    'colors'     - Cell array of colors
%    'linewidth'  - Line width (default: 2)
%    'grid'       - Show grid (default: true)
%    'legend'     - Show legend (default: true if multiple plots)
%    'markers'    - Show start/end markers (default: true)
%    'saveas'     - Save to file (e.g., 'convergence.png')
%
% Outputs:
%    fig          - Figure handle
%
% Examples:
%    % Single convergence plot
%    plot_convergence(BestCosts);
%
%    % Compare multiple runs
%    plot_convergence(run1_costs, run2_costs, run3_costs, ...
%                     'labels', {'IWO', 'GA', 'ABC'});
%
%    % Linear scale with custom title
%    plot_convergence(BestCosts, 'scale', 'linear', ...
%                     'title', 'My Experiment');
%
% See also: plot_comparison, plot_correlation, generate_report
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Inputs
    % Find where name-value pairs start
    nvPairStart = 1;
    costData = {};

    for i = 1:nargin
        if ischar(varargin{i}) || isstring(varargin{i})
            nvPairStart = i;
            break;
        end
        costData{end+1} = varargin{i}; %#ok<AGROW>
    end

    % Parse name-value pairs
    p = inputParser;
    addParameter(p, 'title', 'Optimization Convergence', @ischar);
    addParameter(p, 'xlabel', 'Iteration', @ischar);
    addParameter(p, 'ylabel', 'Best Cost', @ischar);
    addParameter(p, 'scale', 'semilogy', @(x) any(validatestring(x, {'linear', 'semilogy'})));
    addParameter(p, 'labels', {}, @iscell);
    addParameter(p, 'colors', {}, @iscell);
    addParameter(p, 'linewidth', 2, @isnumeric);
    addParameter(p, 'grid', true, @islogical);
    addParameter(p, 'legend', [], @(x) islogical(x) || isempty(x));
    addParameter(p, 'markers', true, @islogical);
    addParameter(p, 'saveas', '', @ischar);

    if nvPairStart <= nargin
        parse(p, varargin{nvPairStart:end});
    else
        parse(p);
    end

    plotTitle = p.Results.title;
    xLabel = p.Results.xlabel;
    yLabel = p.Results.ylabel;
    plotScale = p.Results.scale;
    labels = p.Results.labels;
    colors = p.Results.colors;
    lineWidth = p.Results.linewidth;
    showGrid = p.Results.grid;
    showLegend = p.Results.legend;
    showMarkers = p.Results.markers;
    saveFile = p.Results.saveas;

    %% Validate Input
    if isempty(costData)
        error('PlotConvergence:NoData', 'At least one cost vector is required');
    end

    nPlots = length(costData);

    % Auto-generate labels if not provided
    if isempty(labels)
        if nPlots == 1
            labels = {'Best Cost'};
        else
            labels = cell(1, nPlots);
            for i = 1:nPlots
                labels{i} = sprintf('Run %d', i);
            end
        end
    elseif length(labels) ~= nPlots
        warning('PlotConvergence:LabelMismatch', ...
                'Number of labels (%d) does not match number of plots (%d)', ...
                length(labels), nPlots);
        % Pad or truncate
        if length(labels) < nPlots
            for i = (length(labels)+1):nPlots
                labels{i} = sprintf('Run %d', i);
            end
        end
    end

    % Default colors
    if isempty(colors)
        defaultColors = {'b', 'r', 'g', 'm', 'c', 'k', ...
                        [0.8 0.4 0], [0.4 0 0.8], [0 0.6 0.4]};
        colors = defaultColors(mod(0:nPlots-1, length(defaultColors)) + 1);
    end

    % Auto legend decision
    if isempty(showLegend)
        showLegend = (nPlots > 1);
    end

    %% Create Figure
    fig = figure('Name', plotTitle, 'NumberTitle', 'off', ...
                 'Position', [100 100 800 600]);

    hold on;

    %% Plot Data
    plotHandles = zeros(1, nPlots);

    for i = 1:nPlots
        costs = costData{i};

        % Validate data
        if ~isvector(costs)
            warning('PlotConvergence:InvalidData', ...
                    'Cost data %d is not a vector, skipping', i);
            continue;
        end

        costs = costs(:);  % Ensure column vector
        iterations = 1:length(costs);

        % Plot based on scale
        if strcmp(plotScale, 'semilogy')
            plotHandles(i) = semilogy(iterations, costs, ...
                                     'Color', colors{i}, ...
                                     'LineWidth', lineWidth, ...
                                     'DisplayName', labels{i});
        else
            plotHandles(i) = plot(iterations, costs, ...
                                 'Color', colors{i}, ...
                                 'LineWidth', lineWidth, ...
                                 'DisplayName', labels{i});
        end

        % Add markers for start and end
        if showMarkers
            % Start marker
            plot(1, costs(1), 'o', ...
                 'Color', colors{i}, ...
                 'MarkerFaceColor', colors{i}, ...
                 'MarkerSize', 8, ...
                 'HandleVisibility', 'off');

            % End marker
            plot(length(costs), costs(end), '*', ...
                 'Color', colors{i}, ...
                 'MarkerSize', 10, ...
                 'LineWidth', 2, ...
                 'HandleVisibility', 'off');
        end
    end

    hold off;

    %% Format Plot
    xlabel(xLabel, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(yLabel, 'FontSize', 12, 'FontWeight', 'bold');
    title(plotTitle, 'FontSize', 14, 'FontWeight', 'bold');

    if showGrid
        grid on;
    end

    if showLegend && nPlots > 0
        leg = legend(plotHandles(plotHandles ~= 0), 'Location', 'best');
        leg.FontSize = 10;
    end

    % Improve axes
    ax = gca;
    ax.FontSize = 10;
    ax.Box = 'on';

    %% Save if requested
    if ~isempty(saveFile)
        [~, ~, ext] = fileparts(saveFile);
        if isempty(ext)
            saveFile = [saveFile '.png'];
        end

        try
            saveas(fig, saveFile);
            fprintf('✓ Plot saved to: %s\n', saveFile);
        catch ME
            warning('PlotConvergence:SaveError', ...
                    'Failed to save plot: %s', ME.message);
        end
    end
end
