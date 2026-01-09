%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA119
% Project Title: Implementation of Invasive Weed Optimization in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
%
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
%
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

global in ou t;
in=in_H;
ou=out_H;
t=time;
% t=t';
%% Problem Definition

CostFunction = @(x) Sphere(x,in,ou,t);  % Objective Function

nVar = 6;           % Number of Decision Variables
VarSize = [1 nVar]; % Decision Variables Matrix Size

                    VarMin = [ -200 -200 -200 -200 -200 -200];
                    VarMax =[200 200 200 200 200 200]  ;               % Upper Bound of Decision Variables
%% IWO Parameters

MaxIt = 500;    % Maximum Number of Iterations

nPop0 = 5;     % Initial Population Size
nPop = 20;      % Maximum Population Size

Smin = 2;       % Minimum Number of Seeds
Smax = 7;       % Maximum Number of Seeds

Exponent = 0.5;           % Variance Reduction Exponent
sigma_initial = 0.9;    % Initial Value of Standard Deviation
sigma_final = 0.001;	% Final Value of Standard Deviation

%% Initialization

% Empty Plant Structure
empty_plant.Position = [];
empty_plant.Cost = [];

pop = repmat(empty_plant, nPop0, 1);    % Initial Population Array

for i = 1:numel(pop)
%     if i==1
%       pop(i).Position =[-0.318654622287613,-17.8453381426888,-0.993998295501291,30.0546962349369,-0.0433158223659158,-0.334791393471090,211.623545732571,0,0.100000000000000,-0.0974421245437278,99.4487350826012,0,0.120286733634098,-0.0364429419846437,5.87927234911562e-05,0.0469387166324384,-0.507832782032205,-11.3842492471931,79.0519781190210,0.0255219773194581,-0.578720868835278,0.000978027055460892,-1.09512609205338,0,-7.65712691407085,71.8741055136934,-1.58712756290823,1.49453531166143,0.0414907198741265,0,0,-0.100000000000000,0.994234839455727,0.621657054488682,-0.0861540909996352,-81.4739396640063,-37.6390629641359,8.66095448573757,0.961950586660375,-0.649161518199722];
%    pop(i).Cost = CostFunction(pop(i).Position);
%     else
      % Initialize Position
    pop(i).Position = unifrnd(VarMin,VarMax);


    % Evaluation
    pop(i).Cost = CostFunction(pop(i).Position);

end
% pop(1).Position=[-0.0473064690979077,59.4041049445030,-0.296884284761468,-56.2215819489845,0,-0.518282228913252,120,0,0.0907165029098635,0.0979909431328154,96.8053646034147,0,0.287820031408099,-0.103159686663789,0.696835091490538,0.985710497955199,-0.445529600673732,6.84257706236652,12.1818231714891,0.152009890513607,1.79127797360832,0.00324660882419941,0.879634078320579,0,-2.73801255483213,-52.0307324907271,0.159480671534981,-17.0261020023028,0.177668710695164,0,0,-0.0105976859805020,0.0670731474180124,0.123792286988724,-0.0623295610193758,70.8489177449588,37.4610946299896,-1.48779112825753,-0.773922780669298,-0.793421860213022];
% pop(1).Cost = CostFunction(pop(1).Position);
% Initialize Best Cost History
BestCosts = zeros(MaxIt, 1);

%% IWO Main Loop

for it = 1:MaxIt
%     if it==round(MaxIt*0.75)
%     pop(1).Position=[-0.0473064690979077,59.4041049445030,-0.296884284761468,-56.2215819489845,0,-0.518282228913252,120,0,0.0907165029098635,0.0979909431328154,96.8053646034147,0,0.287820031408099,-0.103159686663789,0.696835091490538,0.985710497955199,-0.445529600673732,6.84257706236652,12.1818231714891,0.152009890513607,1.79127797360832,0.00324660882419941,0.879634078320579,0,-2.73801255483213,-52.0307324907271,0.159480671534981,-17.0261020023028,0.177668710695164,0,0,-0.0105976859805020,0.0670731474180124,0.123792286988724,-0.0623295610193758,70.8489177449588,37.4610946299896,-1.48779112825753,-0.773922780669298,-0.793421860213022];
%     pop(1).Cost = CostFunction(pop(1).Position);
%     end
    % Update Standard Deviation
    sigma = ((MaxIt - it)/(MaxIt - 1))^Exponent * (sigma_initial - sigma_final) + sigma_final;

    % Get Best and Worst Cost Values
    Costs = [pop.Cost];
    BestCost = min(Costs);
    WorstCost = max(Costs);

    % Initialize Offsprings Population
    newpop = [];

    % Reproduction
    for i = 1:numel(pop)

        ratio = (pop(i).Cost - WorstCost)/(BestCost - WorstCost);
        S = floor(Smin + (Smax - Smin)*ratio);

        for j = 1:S

            % Initialize Offspring
            newsol = empty_plant;

            % Generate Random Location
            newsol.Position = pop(i).Position + sigma * randn(VarSize);

            % Apply Lower/Upper Bounds
            newsol.Position = max(newsol.Position, VarMin);
            newsol.Position = min(newsol.Position, VarMax);

            % Evaluate Offsring
            newsol.Cost = CostFunction(newsol.Position);

            % Add Offpsring to the Population
            newpop = [newpop
                      newsol];  %#ok

        end

    end

    % Merge Populations
    pop = [pop
           newpop];

    % Sort Population
    [~, SortOrder]=sort([pop.Cost]);
    pop = pop(SortOrder);

    % Competitive Exclusion (Delete Extra Members)
    if numel(pop)>nPop
        pop = pop(1:nPop);
    end

    % Store Best Solution Ever Found
    BestSol = pop(1);

    % Store Best Cost History
    BestCosts(it) = BestSol.Cost;

    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);

end

%% Results

figure;
% plot(BestCosts,'LineWidth',2);
semilogy(BestCosts,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

% %%Validation
% x= BestSol.position;
% x=x';
% Xu=x(1,1);
% Xa=x(2,1);
% Yv=x(3,1);
% Yb=x(4,1);
% Lu=x(5,1);
% Lv=x(6,1);
% Lb=x(7,1);
% Lw=x(8,1);
% Mu=x(9,1);
% Mv=x(10,1);
% Ma=x(11,1);
% Mw=x(12,1);
% Tf=x(13,1);
% Ab=x(14,1);
% Ac=x(15,1);
% Ba=x(16,1);
% Bd=x(17,1);
% Za=x(18,1);
% Zb=x(19,1);
% Zw=x(20,1);
% Zr=x(21,1);
%
% Nv=x(22,1);
% Np=x(23,1);
% Nw=x(24,1);
% Nr=x(25,1);
% Nrfb=x(26,1);
% Kr=x(27,1);
% Krfb=x(28,1);
% Ts=x(29,1);
%
% Yped=x(30,1);
% Mcol=x(31,1);
% Alat=x(32,1);
% Alon=x(33,1);
% Blat=x(34,1);
% Blon=x(35,1);
% Zcol=x(36,1);
% Nped=x(37,1);
% Ncol=x(38,1);
% Clon=x(39,1);
% Dlat=x(40,1);
%
%  A=[Xu,0,0,0,0,-9.81,Xa,0,0,0,0,0,0;
%     0,Yv,0,0,9.81,0,0,Yb,0,0,0,0,0;
%     Lu,Lv,0,0,0,0,0,Lb,Lw,0,0,0,0;
%     Mu,Mv,0,0,0,0,Ma,0,Mw,0,0,0,0;
%     0,0,1,0,0,0,0,0,0,0,0,0,0;
%     0,0,0,1,0,0,0,0,0,0,0,0,0;
%     0,0,0,-Tf,0,0,-1,Ab,0,0,0,Ac,0;
%     0,0,-Tf,0,0,0,Ba,-1,0,0,0,0,Bd;
%     0,0,0,0,0,0,Za,Zb,Zw,Zr,0,0,0;
%     0,Nv,Np,0,0,0,0,0,Nw,Nr,Nrfb,0,0;
%     0,0,0,0,0,0,0,0,0,Kr,Krfb,0,0;
%     0,0,0,-Ts,0,0,0,0,0,0,0,-1,0;
%     0,0,-Ts,0,0,0,0,0,0,0,0,0,-1];
%
%  B=[0,0,0,0;
%     0,0,Yped,0;
%     0,0,0,0;
%     0,0,0,Mcol;
%     0,0,0,0;
%     0,0,0,0;
%     Alat,Alon,0,0;
%     Blat,Blon,0,0;
%     0,0,0,Zcol;
%     0,0,Nped,Ncol;
%     0,0,0,0;
%     0,Clon,0,0;
%     Dlat,0,0,0];
% C=eye(13);
% D=0;
% mod=ss(A,B,C,D);
%
%
%  y(:,:,1)=lsim(mod,in,t);
%    cor1=corrcoef(y(:,1,1),ou(:,1));
%    cor2=corrcoef(y(:,2,1),ou(:,2));
%  cor3=corrcoef(y(:,3,1),ou(:,3));
% cor4=corrcoef(y(:,4,1),ou(:,4));
%   cor5=corrcoef(y(:,5,1),ou(:,5));
%  cor6=corrcoef(y(:,6,1),ou(:,6));
%  cor7=corrcoef(y(:,7,1),ou(:,7));
%  cor8=corrcoef(y(:,8,1),ou(:,8));
%  cor9=corrcoef(y(:,9,1),ou(:,9));
%  cor10=corrcoef(y(:,10,1),ou(:,10));
