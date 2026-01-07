function animate_convergence(results, varargin)
% ANIMATE_CONVERGENCE - Create animated visualization of optimization convergence
%
% This function creates an animated plot showing the convergence of the IWO
% algorithm over iterations. Useful for presentations and analysis.
%
% Usage:
%    animate_convergence(results)
%    animate_convergence(results, 'speed', 0.05)
%    animate_convergence(results, 'saveas', 'convergence.gif')
%
% Inputs:
%    results    - Results structure from IWO or run_experiment
%
% Optional Parameters (Name-Value pairs):
%    'speed'    - Animation speed in seconds per frame (default: 0.1)
%    'scale'    - Y-axis scale: 'linear' or 'log' (default: 'log')
%    'title'    - Plot title (default: 'IWO Convergence Animation')
%    'saveas'   - Save animation to file (.gif or .avi) (default: '')
%    'quality'  - Video quality: 'low', 'medium', 'high' (default: 'medium')
%    'fps'      - Frames per second for video (default: 10)
%    'showBest' - Show best cost value in legend (default: true)
%
% Examples:
%    % Basic animation
%    animate_convergence(results);
%
%    % Save as GIF
%    animate_convergence(results, 'saveas', 'convergence.gif');
%
%    % Save as video with high quality
%    animate_convergence(results, 'saveas', 'convergence.avi', 'quality', 'high');
%
% See also: plot_convergence, animate_population
%
% Author: System Identification Project
% Date: 2026-01-07

    %% Parse Inputs
    p = inputParser;
    addRequired(p, 'results', @isstruct);
    addParameter(p, 'speed', 0.1, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'scale', 'log', @(x) any(validatestring(x, {'linear', 'log'})));
    addParameter(p, 'title', 'IWO Convergence Animation', @ischar);
    addParameter(p, 'saveas', '', @ischar);
    addParameter(p, 'quality', 'medium', @(x) any(validatestring(x, {'low', 'medium', 'high'})));
    addParameter(p, 'fps', 10, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'showBest', true, @islogical);
    parse(p, results, varargin{:});

    speed = p.Results.speed;
    scale = p.Results.scale;
    titleStr = p.Results.title;
    saveFile = p.Results.saveas;
    quality = p.Results.quality;
    fps = p.Results.fps;
    showBest = p.Results.showBest;

    %% Extract Best Costs
    if isfield(results, 'bestCosts')
        bestCosts = results.bestCosts;
    elseif isfield(results, 'BestCosts')
        bestCosts = results.BestCosts;
    else
        error('AnimateConvergence:MissingData', 'Results must contain bestCosts field');
    end

    nIterations = length(bestCosts);

    fprintf('Creating convergence animation...\n');
    fprintf('  Iterations: %d\n', nIterations);
    fprintf('  Speed: %.2f sec/frame\n', speed);

    %% Setup Figure
    fig = figure('Position', [100, 100, 800, 600], 'Color', 'white');

    %% Determine Frame Interval
    % For large datasets, skip frames to keep animation manageable
    if nIterations > 500
        frameInterval = ceil(nIterations / 500);
    else
        frameInterval = 1;
    end

    frames = 1:frameInterval:nIterations;
    if frames(end) ~= nIterations
        frames = [frames, nIterations];  % Always include final frame
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

            % Set video quality
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
            error('AnimateConvergence:InvalidFormat', 'Save format must be .gif or .avi');
        end
    end

    %% Create Animation
    for i = 1:length(frames)
        idx = frames(i);

        % Plot data up to current iteration
        if strcmp(scale, 'log')
            semilogy(1:idx, bestCosts(1:idx), 'b-', 'LineWidth', 2);
        else
            plot(1:idx, bestCosts(1:idx), 'b-', 'LineWidth', 2);
        end

        hold on;

        % Mark current point
        if strcmp(scale, 'log')
            semilogy(idx, bestCosts(idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        else
            plot(idx, bestCosts(idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        end

        hold off;

        % Labels and title
        xlabel('Iteration', 'FontSize', 12);
        if strcmp(scale, 'log')
            ylabel('Best Cost (log scale)', 'FontSize', 12);
        else
            ylabel('Best Cost', 'FontSize', 12);
        end

        if showBest
            title(sprintf('%s\nIteration: %d/%d | Best Cost: %.6f', ...
                         titleStr, idx, nIterations, bestCosts(idx)), 'FontSize', 14);
        else
            title(sprintf('%s\nIteration: %d/%d', titleStr, idx, nIterations), 'FontSize', 14);
        end

        grid on;
        xlim([0, nIterations * 1.05]);

        % Set consistent y-limits
        if strcmp(scale, 'log')
            ylim([min(bestCosts) * 0.5, max(bestCosts) * 2]);
        else
            yRange = max(bestCosts) - min(bestCosts);
            ylim([min(bestCosts) - yRange*0.1, max(bestCosts) + yRange*0.1]);
        end

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
