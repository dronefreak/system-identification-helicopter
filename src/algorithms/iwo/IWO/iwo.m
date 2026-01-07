%
% IWO - Invasive Weed Optimization for Helicopter System Identification
%
% This script implements the Invasive Weed Optimization algorithm to identify
% 40 aerodynamic and control parameters of an unmanned helicopter state-space
% model. The optimization maximizes the correlation between simulated and
% actual flight test data.
%
% Required workspace variables (load from best.mat):
%    in_H     - Input data for helicopter (4 control inputs × time steps)
%    out_H    - Output data for helicopter (10 measurements × time steps)
%    time     - Time vector for simulation
%    inr      - Input matrix (required by cost function)
%    outr     - Output matrix (required by cost function)
%
% Outputs:
%    BestSol     - Structure containing best solution found
%                  .Position: 40 optimized parameters
%                  .Cost: Final cost value
%    BestCosts   - Array of best cost at each iteration
%    pop         - Final population
%
% Usage:
%    load('best.mat')
%    iwo
%
% Original Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% Modified for helicopter system identification
% Author: S. Mostapha Kalami Heris, Modified by Project Team
% Date: 2026-01-06
%

%% Validate Required Workspace Variables
requiredVars = {'in_H', 'out_H', 'time', 'inr', 'outr'};
missingVars = {};

for i = 1:length(requiredVars)
    if ~evalin('base', ['exist(''' requiredVars{i} ''', ''var'')'])
        missingVars{end+1} = requiredVars{i}; %#ok<AGROW>
    end
end

if ~isempty(missingVars)
    error('IWO:MissingVariables', ...
          'Missing required workspace variables: %s\nPlease load best.mat first.', ...
          strjoin(missingVars, ', '));
end

%% Load Configuration
config = config_iwo();

%% Prepare Data (avoiding global variables)
inputData = in_H;
outputData = out_H;
timeVector = time;

%% Problem Definition
CostFunction = @(x) Sphere(x, inputData, outputData, timeVector);

nVar = config.nVar;
varSize = [1 nVar];
varMin = config.varMin;
varMax = config.varMax;

%% Validate Parameter Bounds
if length(varMin) ~= nVar || length(varMax) ~= nVar
    error('IWO:InvalidBounds', ...
          'Parameter bounds must have length %d (nVar)', nVar);
end

%% IWO Algorithm Parameters
maxIterations = config.maxIterations;
initialPopSize = config.initialPopSize;
maxPopSize = config.maxPopSize;
minSeeds = config.minSeeds;
maxSeeds = config.maxSeeds;
varianceExponent = config.varianceExponent;
sigmaInitial = config.sigmaInitial;
sigmaFinal = config.sigmaFinal;

%% Performance Options
useParallel = config.useParallel;
showProgress = config.showProgress;
progressUpdateInterval = config.progressUpdateInterval;
preallocateOffspring = config.preallocateOffspring;

% Check for Parallel Computing Toolbox
if useParallel
    if license('test', 'Distrib_Computing_Toolbox')
        try
            % Try to get current pool, create if doesn't exist
            poolObj = gcp('nocreate');
            if isempty(poolObj)
                fprintf('Starting parallel pool...\n');
                poolObj = parpool;
            end
            fprintf('Using parallel pool with %d workers\n', poolObj.NumWorkers);
        catch ME
            warning('IWO:ParallelError', ...
                    'Could not start parallel pool: %s\nFalling back to serial execution', ...
                    ME.message);
            useParallel = false;
        end
    else
        warning('IWO:NoParallelToolbox', ...
                'Parallel Computing Toolbox not available. Using serial execution.');
        useParallel = false;
    end
end

%% Initialization
fprintf('Starting Invasive Weed Optimization...\n');
fprintf('Parameters: %d variables, %d iterations, population %d-%d\n', ...
        nVar, maxIterations, initialPopSize, maxPopSize);

% Empty plant structure template
emptyPlant.Position = [];
emptyPlant.Cost = [];

% Initialize population array
population = repmat(emptyPlant, initialPopSize, 1);

% Initialize each plant with random position
fprintf('Initializing population...\n');

if useParallel
    % Parallel initialization
    positions = zeros(initialPopSize, nVar);
    costs = zeros(initialPopSize, 1);

    for i = 1:initialPopSize
        positions(i, :) = unifrnd(varMin, varMax);
    end

    parfor i = 1:initialPopSize
        try
            costs(i) = CostFunction(positions(i, :));
            if isnan(costs(i)) || ~isreal(costs(i))
                costs(i) = Inf;
            end
        catch
            costs(i) = Inf;
        end
    end

    % Assign to population
    for i = 1:initialPopSize
        population(i).Position = positions(i, :);
        population(i).Cost = costs(i);
    end
else
    % Serial initialization
    for i = 1:initialPopSize
        % Initialize position with uniform random values within bounds
        population(i).Position = unifrnd(varMin, varMax);

        % Evaluate cost function
        try
            population(i).Cost = CostFunction(population(i).Position);
        catch ME
            warning('IWO:InitializationError', ...
                    'Error evaluating initial population member %d: %s', i, ME.message);
            population(i).Cost = Inf;
        end

        % Check for invalid cost
        if isnan(population(i).Cost) || ~isreal(population(i).Cost)
            warning('IWO:InvalidCost', ...
                    'Invalid cost for population member %d, setting to Inf', i);
            population(i).Cost = Inf;
        end
    end
end

% Initialize best cost history
bestCosts = zeros(maxIterations, 1);

fprintf('Initialization complete. Starting optimization...\n\n');

%% Initialize Progress Tracking
if showProgress
    progressBar = waitbar(0, 'IWO Optimization: Starting...', ...
                          'Name', 'IWO Progress');
end
startTime = tic;

%% IWO Main Loop
for iteration = 1:maxIterations
    % Update standard deviation (decreases over iterations for convergence)
    sigma = ((maxIterations - iteration) / (maxIterations - 1))^varianceExponent ...
            * (sigmaInitial - sigmaFinal) + sigmaFinal;

    % Get best and worst cost values in current population
    costs = [population.Cost];
    bestCost = min(costs);
    worstCost = max(costs);

    % Calculate total number of offspring for preallocation
    totalSeeds = 0;
    seedCounts = zeros(numel(population), 1);

    for i = 1:numel(population)
        % Calculate fitness ratio (better plants produce more seeds)
        if worstCost ~= bestCost
            ratio = (population(i).Cost - worstCost) / (bestCost - worstCost);
        else
            ratio = 0.5;  % All plants equal, use middle value
        end

        % Determine number of seeds (better plants → more seeds)
        numSeeds = floor(minSeeds + (maxSeeds - minSeeds) * ratio);
        seedCounts(i) = numSeeds;
        totalSeeds = totalSeeds + numSeeds;
    end

    % Preallocate offspring array for memory efficiency
    if preallocateOffspring && totalSeeds > 0
        offspringPositions = zeros(totalSeeds, nVar);
        offspringCosts = zeros(totalSeeds, 1);
    else
        offspring = [];
    end

    % Generate all offspring positions first
    offspringIdx = 1;
    for i = 1:numel(population)
        numSeeds = seedCounts(i);

        for j = 1:numSeeds
            % Generate new position: parent position + random variation
            newPosition = population(i).Position + sigma * randn(varSize);

            % Apply bounds
            newPosition = max(newPosition, varMin);
            newPosition = min(newPosition, varMax);

            if preallocateOffspring
                offspringPositions(offspringIdx, :) = newPosition;
                offspringIdx = offspringIdx + 1;
            else
                newSolution = emptyPlant;
                newSolution.Position = newPosition;
                offspring = [offspring; newSolution]; %#ok<AGROW>
            end
        end
    end

    % Evaluate all offspring (parallel or serial)
    if preallocateOffspring && totalSeeds > 0
        if useParallel
            % Parallel evaluation
            parfor i = 1:totalSeeds
                try
                    offspringCosts(i) = CostFunction(offspringPositions(i, :));
                    if isnan(offspringCosts(i)) || ~isreal(offspringCosts(i))
                        offspringCosts(i) = Inf;
                    end
                catch
                    offspringCosts(i) = Inf;
                end
            end
        else
            % Serial evaluation
            for i = 1:totalSeeds
                try
                    offspringCosts(i) = CostFunction(offspringPositions(i, :));
                catch ME
                    warning('IWO:EvaluationError', ...
                            'Error evaluating offspring at iteration %d: %s', ...
                            iteration, ME.message);
                    offspringCosts(i) = Inf;
                end

                if isnan(offspringCosts(i)) || ~isreal(offspringCosts(i))
                    offspringCosts(i) = Inf;
                end
            end
        end

        % Create offspring structure array
        offspring = repmat(emptyPlant, totalSeeds, 1);
        for i = 1:totalSeeds
            offspring(i).Position = offspringPositions(i, :);
            offspring(i).Cost = offspringCosts(i);
        end
    elseif ~preallocateOffspring
        % Evaluate offspring one by one (old method)
        for i = 1:numel(offspring)
            try
                offspring(i).Cost = CostFunction(offspring(i).Position);
            catch ME
                warning('IWO:EvaluationError', ...
                        'Error evaluating offspring at iteration %d: %s', ...
                        iteration, ME.message);
                offspring(i).Cost = Inf;
            end

            if isnan(offspring(i).Cost) || ~isreal(offspring(i).Cost)
                offspring(i).Cost = Inf;
            end
        end
    end

    % Merge parent and offspring populations
    population = [population; offspring];

    % Sort population by cost (ascending - lower is better)
    [~, sortOrder] = sort([population.Cost]);
    population = population(sortOrder);

    % Competitive exclusion: keep only the best maxPopSize plants
    if numel(population) > maxPopSize
        population = population(1:maxPopSize);
    end

    % Store best solution found so far
    bestSolution = population(1);

    % Record best cost for this iteration
    bestCosts(iteration) = bestSolution.Cost;

    % Update progress bar
    if showProgress && (mod(iteration, progressUpdateInterval) == 0 || iteration == 1 || iteration == maxIterations)
        progress = iteration / maxIterations;
        elapsedTime = toc(startTime);
        estimatedTotal = elapsedTime / progress;
        remainingTime = estimatedTotal - elapsedTime;

        if ishandle(progressBar)
            waitbar(progress, progressBar, ...
                    sprintf('Iteration %d/%d | Best: %.6f | ETA: %.1f min', ...
                            iteration, maxIterations, bestCosts(iteration), ...
                            remainingTime / 60));
        end
    end

    % Display progress
    if mod(iteration, config.displayInterval) == 0 || iteration == 1
        fprintf('Iteration %4d/%d: Best Cost = %.6f, Sigma = %.6f\n', ...
                iteration, maxIterations, bestCosts(iteration), sigma);
    end
end

%% Cleanup Progress Bar
if showProgress && ishandle(progressBar)
    close(progressBar);
end

%% Store Results in Base Workspace
fprintf('\nOptimization complete!\n');
fprintf('Final best cost: %.6f\n', bestSolution.Cost);

% Assign results to base workspace for user access
assignin('base', 'BestSol', bestSolution);
assignin('base', 'BestCosts', bestCosts);
assignin('base', 'pop', population);

%% Results Visualization
if config.plotResults
    figure('Name', 'IWO Convergence', 'NumberTitle', 'off');

    if config.useSemilogy
        semilogy(bestCosts, 'LineWidth', 2);
    else
        plot(bestCosts, 'LineWidth', 2);
    end

    xlabel('Iteration');
    ylabel('Best Cost');
    title('IWO Convergence History');
    grid on;

    % Add annotations
    hold on;
    plot(1, bestCosts(1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(maxIterations, bestCosts(end), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    legend('Cost History', 'Start', 'End', 'Location', 'best');
    hold off;
end

fprintf('\nResults stored in workspace:\n');
fprintf('  BestSol    - Best solution structure (.Position contains 40 parameters)\n');
fprintf('  BestCosts  - Convergence history\n');
fprintf('  pop        - Final population\n');
fprintf('\nTo visualize results, run: PopulCheck\n');
