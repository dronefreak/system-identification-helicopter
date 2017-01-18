% Demonstration of different neural network training algorithms used
% for curve fitting

close all
%---------- Generate training and test set ----------
clc
PHI=0:0.25:6;
Y=sin(PHI);
PHI1 = PHI(:,1:2:25);
Y1   = Y(:,1:2:25);
PHI2 = PHI(:,2:2:24);
Y2   = Y(:,2:2:24);
figure(1);
plot(PHI1,Y1,'+')
title('Data set for training')
drawnow


% ----------- Training the network ----------
k=5;
while k~=7,

  k = menu('Choose one of the training algorithms below:',...
            'Back Propagation (Batch)',...
            'Back Propagation (Recursive)',...
            'Recursive Prediction Error Method (forgetting factor)',...
            'Recursive Prediction Error Method (Constant Trace)',...
            'Recursive Prediction Error Method (EFRA)',...
            'Levenberg-Marquardt method',...
            'quit');
         
  chvec = get(0,'Children'); % Close all figures except for figure(1)
  chvec = chvec(find(chvec~=1));       
  close(chvec);
  %---------- Define network structure and initialize weights ----------
  rand('seed',0);
  W1         = rand(4,2);  % Weights to hidden layer 
  W2         = rand(1,5);  % Weights to output
  NetDef = ['HHHH'       
            'L---'];       % The top row represents the structure of 
                           % the hidden layer. The bottom row represents
                           % the output structure.
  trparms = settrain;      % Set training parameters to default values

  % ----- Back propagation (Batch) -----
  if k==1,
    trparms=settrain(trparms,'maxiter',1000,'eta',0.01);
    [W1,W2,PI_vector,iter]=batbp(NetDef,W1,W2,PHI1,Y1,trparms);

  
  % ----- Back propagation (Recursive) -----
  elseif k==2,
    trparms=settrain(trparms,'maxiter',1000,'eta',0.01);
    [W1,W2,PI_vector,iter]=incbp(NetDef,W1,W2,PHI1,Y1,trparms);
    

  % ----- RPE algorithm (Forgetting factor) -----
  elseif k==3,
    trparms=settrain(trparms,'maxiter',200,'p0',200,'lambda',0.98);
    [W1,W2,PI_vector,iter]=rpe(NetDef,W1,W2,PHI1,Y1,trparms);

    
  % ----- RPE algorithm (Constant Trace) -----
  elseif k==4,
    trparms=settrain(trparms,'method','ct');
    trparms=settrain(trparms,'maxiter',200,'alpha_max',200,'alpha_min',0.001);
    [W1,W2,PI_vector,iter]=rpe(NetDef,W1,W2,PHI1,Y1,trparms);

  
  % ----- RPE algorithm (EFRA) -----
  elseif k==5,  
    trparms=settrain(trparms,'method','efra','maxiter',200,'alpha',1);
    [W1,W2,PI_vector,iter]=rpe(NetDef,W1,W2,PHI1,Y1,trparms);

  
  % ----- Marquardt algorithm -----
  elseif k==6,
    trparms=settrain(trparms,'maxiter',200);
    [W1,W2,PI_vector,iter,lambda]=marq(NetDef,W1,W2,PHI1,Y1,trparms);
  end

  if k~=7,
    
    % -----------  Validate Network  -----------
    [Y_sim,E,PI] = nneval(NetDef,W1,W2,PHI2,Y2);
    
    % -----------  Plot Cost function  -----------
    figure
    semilogy(PI_vector)
    title('Criterion evaluated after each iteration')
    xlabel('Iteration (epoch)')
    ylabel('Criterion')
    grid
  end
end
