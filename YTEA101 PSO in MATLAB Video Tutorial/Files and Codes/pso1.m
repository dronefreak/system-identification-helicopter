%
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YTEA101
% Project Title: Particle Swarm Optimization Video Tutorial
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer and Instructor: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

%% Problem Definiton
global in ou;
in=inr;
ou=outr;


problem.CostFunction = @(x1) Sphere(x1,in,ou,t);  % Cost Function
problem.nVar = 40;       % Number of Unknown (Decision) Variables
problem.VarMin =  -1000;  % Lower Bound of Decision Variables
problem.VarMax =  1000;   % Upper Bound of Decision Variables

%% Parameters of PSO

params.MaxIt = 10000;        % Maximum Number of Iterations
params.nPop = 20;           % Population Size (Swarm Size)
params.w = 1;               % Intertia Coefficient
params.wdamp = 0.99;        % Damping Ratio of Inertia Coefficient
params.c1 = 5;              % Personal Acceleration Coefficient
params.c2 = 3;              % Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

%% Calling PSO

 out = PSO(problem, params,population);

BestSol = out.BestSol;
BestCosts = out.BestCosts;

%% Results

figure;
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;


