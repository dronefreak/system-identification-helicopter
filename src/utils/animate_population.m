function animate_population(results, varargin)
% ANIMATE_POPULATION - Animate population evolution in parameter space
%
% This function creates an animated visualization showing how the IWO
% population evolves in parameter space over iterations.
%
% Usage:
%    animate_population(results)
%    animate_population(results, 'params', [1, 2])
%    animate_population(results, 'saveas', 'population.gif')
%
% Inputs:
%    results    - Results structure from IWO with population history
%
% Optional Parameters (Name-Value pairs):
%    'params'      - Parameter indices to plot [param1, param2] (default: [1, 2])
%    'speed'       - Animation speed in seconds per frame (default: 0.1)
%    'title'       - Plot title (default: 'Population Evolution')
%    'saveas'      - Save animation to file (.gif or .avi) (default: '')
%    'quality'     - Video quality: 'low', 'medium', 'high' (default: 'medium')
%    'fps'         - Frames per second for video (default: 10)
%    'markerSize'  - Size of population markers (default: 50)
%    'showBest'    - Highlight best solution (default: true)
%    'showTrails'  - Show movement trails (default: false)
%    'trailLength' - Number of iterations for trails (default: 5)
%
% Note: This function requires the results structure to contain population
%       history, which can be saved by setting config.savePopulation = true
%       in the IWO configuration.
%
% Examples:
%    % Animate first two parameters
%    animate_population(results);
%
%    % Animate specific parameters
%    animate_population(results, 'params', [3, 5]);
%
%    % Save with trails
%    animate_population(results, 'showTrails', true, 'saveas', 'pop.gif');
%
% See also: animate_convergence, plot_convergence
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Inputs
    p = inputParser;
    addRequired(p, 'results', @isstruct);
    addParameter(p, 'params', [1, 2], @(x) isnumeric(x) && length(x) == 2);
    addParameter(p, 'speed', 0.1, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'title', 'Population Evolution', @ischar);
    addParameter(p, 'saveas', '', @ischar);
    addParameter(p, 'quality', 'medium', @(x) any(validatestring(x, {'low', 'medium', 'high'})));
    addParameter(p, 'fps', 10, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'markerSize', 50, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'showBest', true, @islogical);
    addParameter(p, 'showTrails', false, @islogical);
    addParameter(p, 'trailLength', 5, @(x) isnumeric(x) && x > 0);
    parse(p, results, varargin{:});

    paramIndices = p.Results.params;
    speed = p.Results.speed;
    titleStr = p.Results.title;
    saveFile = p.Results.saveas;
    quality = p.Results.quality;
    fps = p.Results.fps;
    markerSize = p.Results.markerSize;
    showBest = p.Results.showBest;
    showTrails = p.Results.showTrails;
    trailLength = p.Results.trailLength;

    %% Check for Population History
    if ~isfield(results, 'populationHistory')
        error('AnimatePopulation:MissingData', ...
              ['Results must contain populationHistory field.\n' ...
               'Set config.savePopulation = true before running IWO.']);
    end

    popHistory = results.populationHistory;
    nIterations = length(popHistory);

    fprintf('Creating population evolution animation...\n');
    fprintf('  Iterations: %d\n', nIterations);
    fprintf('  Parameters: [%d, %d]\n', paramIndices(1), paramIndices(2));

    %% Setup Figure
    fig = figure('Position', [100, 100, 800, 700], 'Color', 'white');

    %% Determine Frame Interval
    if nIterations > 200
        frameInterval = ceil(nIterations / 200);
    else
        frameInterval = 1;
    end

    frames = 1:frameInterval:nIterations;
    if frames(end) ~= nIterations
        frames = [frames, nIterations];
    end

    fprintf('  Animation frames: %d\n', length(frames));

    %% Setup Video Writer (if saving)
    if ~isempty(saveFile)
        [~, ~, ext] = fileparts(saveFile);

        if strcmp(ext, '.gif')
            saveAsGif = true;
            fprintf('  Saving as: %s (GIF)\n', saveFile);
        elseif strcmp(ext, '.avi')
            saveAsGif = false;

            switch quality
                case 'low'
                    videoQuality = 50;
                case 'medium'
                    videoQuality = 75;
                case 'high'
                    videoQuality = 95;
            end

            videoWriter = VideoWriter(saveFile);
            videoWriter.FrameRate = fps;
            videoWriter.Quality = videoQuality;
            open(videoWriter);
            fprintf('  Saving as: %s (AVI, quality=%d)\n', saveFile, videoQuality);
        else
            error('AnimatePopulation:InvalidFormat', 'Save format must be .gif or .avi');
        end
    end

    %% Determine Axis Limits
    allPositions1 = [];
    allPositions2 = [];
    for i = 1:length(popHistory)
        if ~isempty(popHistory{i})
            positions = popHistory{i};
            allPositions1 = [allPositions1; positions(:, paramIndices(1))]; %#ok<AGROW>
            allPositions2 = [allPositions2; positions(:, paramIndices(2))]; %#ok<AGROW>
        end
    end

    range1 = max(allPositions1) - min(allPositions1);
    range2 = max(allPositions2) - min(allPositions2);
    xlims = [min(allPositions1) - range1*0.1, max(allPositions1) + range1*0.1];
    ylims = [min(allPositions2) - range2*0.1, max(allPositions2) + range2*0.1];

    %% Create Animation
    for i = 1:length(frames)
        idx = frames(i);
        clf;

        % Show trails if enabled
        if showTrails && i > 1
            trailStart = max(1, i - trailLength);
            for t = trailStart:i-1
                trailIdx = frames(t);
                if ~isempty(popHistory{trailIdx})
                    positions = popHistory{trailIdx};
                    alpha = (t - trailStart + 1) / trailLength;
                    scatter(positions(:, paramIndices(1)), ...
                           positions(:, paramIndices(2)), ...
                           markerSize/2, 'MarkerEdgeColor', 'none', ...
                           'MarkerFaceColor', [0.7, 0.7, 0.7], ...
                           'MarkerFaceAlpha', alpha * 0.3);
                    hold on;
                end
            end
        end

        % Plot current population
        if ~isempty(popHistory{idx})
            positions = popHistory{idx};
            popSize = size(positions, 1);

            % Color by fitness if available
            if isfield(results, 'costHistory') && idx <= length(results.costHistory)
                costs = results.costHistory{idx};
                scatter(positions(:, paramIndices(1)), ...
                       positions(:, paramIndices(2)), ...
                       markerSize, costs, 'filled', 'MarkerEdgeColor', 'k');
                colorbar;
                colormap('jet');
            else
                scatter(positions(:, paramIndices(1)), ...
                       positions(:, paramIndices(2)), ...
                       markerSize, 'b', 'filled', 'MarkerEdgeColor', 'k');
            end

            hold on;

            % Highlight best solution
            if showBest && isfield(results, 'bestCosts')
                if isfield(results, 'bestSolution')
                    bestPos = results.bestSolution.Position;
                elseif isfield(results, 'BestSol')
                    bestPos = results.BestSol.Position;
                end

                if exist('bestPos', 'var')
                    scatter(bestPos(paramIndices(1)), bestPos(paramIndices(2)), ...
                           markerSize*2, 'r', 'filled', 'MarkerEdgeColor', 'k', ...
                           'LineWidth', 2);
                end
            end

            hold off;

            % Labels
            xlabel(sprintf('Parameter %d', paramIndices(1)), 'FontSize', 12);
            ylabel(sprintf('Parameter %d', paramIndices(2)), 'FontSize', 12);
            title(sprintf('%s\nIteration: %d/%d | Population Size: %d', ...
                         titleStr, idx, nIterations, popSize), 'FontSize', 14);
        end

        grid on;
        xlim(xlims);
        ylim(ylims);
        axis square;

        drawnow;

        % Save frame
        if ~isempty(saveFile)
            frame = getframe(fig);

            if saveAsGif
                im = frame2im(frame);
                [imind, cm] = rgb2ind(im, 256);

                if i == 1
                    imwrite(imind, cm, saveFile, 'gif', 'Loopcount', inf, 'DelayTime', speed);
                else
                    imwrite(imind, cm, saveFile, 'gif', 'WriteMode', 'append', 'DelayTime', speed);
                end
            else
                writeVideo(videoWriter, frame);
            end
        end

        % Pause for animation
        if i < length(frames)
            pause(speed);
        end

        % Progress indicator
        if mod(i, max(1, floor(length(frames)/10))) == 0 || i == length(frames)
            fprintf('  Progress: %.0f%%\n', (i/length(frames))*100);
        end
    end

    %% Cleanup
    if ~isempty(saveFile) && ~saveAsGif
        close(videoWriter);
    end

    fprintf('✓ Animation complete!\n');
    if ~isempty(saveFile)
        fprintf('  Saved to: %s\n', saveFile);
    end
end
