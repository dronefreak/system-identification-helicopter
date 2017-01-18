function [theta_data,PI_vector,FPE_vector,PI_test_vec,deff_vec,pvec]=...
                             obdprune(NetDef,W1,W2,PHI,Y,trparms,prparms,PHI2,Y2)
%  OBDPRUNE
%  --------
%           This function applies the Optimal Brain Damage (OBD) algorithm
%           for pruning ordinary feedforward neural networks.
%
%  CALL:
%   [theta_data,NSSEvec,FPEvec,NSSEtestvec,deff,pvec]=...
%                    obdprune(NetDef,W1,W2,PHI,Y,trparms,prparms,PHI2,Y2)
%
%  INPUT:
%  NetDef, W1, W2,
%  PHI, Y, trparms    : See for example the function MARQ
%  PHI2,Y2 (optional) : Test data. Can be used for pointing out the
%                       optimal network architecture. 
%  prparms (optional) : Parameters associated with the pruning session
%                       prparms = [iter RePercent]
%                       iter      : Max. number of retraining iterations
%                       RePercent : Prune 'RePercent' percent of the
%                                   remaining weights (0 = prune one at a time)
%                       If passed as [], prparms=[50 0] will be used.
% 
%  OUTPUT:
%  theta_data  : Matrix containing all the parameter vectors.
%  NSSEvec     : Vector containing the training error (SSE/2N) after each
%                weight elimination.
%  FPEvec      : Contains the FPE estimate of the average generalization error.
%  NSSEtestvec : Contains the test error.
%  deff        : Contains the "effective" number of weights.
%  pvec        : Index into the above vectors.

%  Programmed by : Magnus Norgaard, IAU/IMM, Technical University of Denmark
%  LastEditDate  : Jan. 7, 2000


%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------
more off
if nargin>7, TestDataFlag = 1;          % Check if test data was given as argument
else TestDataFlag = 0;end
if isempty(prparms),
  prparms=[50 0];
end
iter      = prparms(1);                 % Max. retraining iterations
RePercent = prparms(2);                 % % of remaining weights to prune
[outputs,N] = size(Y);                  % # of outputs and # of data
[hidden,inputs] = size(W1);             % # of hidden units 
inputs=inputs-1;                        % # of inputs
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neuron
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
y1       = zeros(hidden,N);             % Hidden layer outputs
y2       = zeros(outputs,N);            % Network output
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A usefull vector!
index2 = (0:N-1)*outputs;               % Yet another usefull vector
PHI_aug  = [PHI;ones(1,N)];             % Augment PHI with a row containg ones
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters = parameters1 + parameters2; % Total # of weights
ones_h   = ones(hidden+1,1);            % A vector of ones
ones_i   = ones(inputs+1,1);            % Another vector of ones
                                        % Parameter vector containing all weights
theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
reduced0 = reduced;                     % Copy of 'reduced'. Will be constant
theta_data=zeros(parameters,parameters);% Matrix used for collecting theta vectors
theta_data(:,reduced) = theta;          % Insert 'initial' theta
PSI      = zeros(parameters,outputs*N); % Deriv. of each output w.r.t. each weight
p0       = 1e6;                         % Diag. element of H_inv (no weight decay)
H_inv    = p0*eye(reduced);             % Initial inverse Hessian (no weight decay)
Ident    = eye(outputs);                % Identity matrix
PI_vector= zeros(1,reduced);            % A vector containing the collected PI's
FPE_vector= zeros(1,reduced);           % Vector used for collecting FPE estimates
if TestDataFlag,                        % Initializations if a test set exists
  [tmp,N2]    = size(Y2);               % # of data in test set
  ytest1      = zeros(hidden,N2);       % Hidden layer outputs 
  ytest2      = zeros(outputs,N2);      % Network output
  PHI2_aug    = [PHI2;ones(1,N2)];      % Augment PHI with a row containing ones
  PI_test_vec = zeros(1,reduced);       % Collected PI's for the test set
end
deff_vec = zeros(1,reduced);            % The effective number of parameters
minweights = 2;                         % Prune until 'minweights' weights remain
FirstTimeFlag=1;                        % Initialize flag
pr = 0;                                 % Initialize counter
pvec=[];                                % Initialize index vector
if ~exist('trparms') | isempty(trparms) % Default training parameters
  trparms = settrain;
  D0       = trparms.D;
else                                    % User specified values
  if ~isstruct(trparms),
     error('''trparms'' must be a structure variable.');
  end
  if ~isfield(trparms,'infolevel')
     trparms = settrain(trparms,'infolevel','default');
  end
  if ~isfield(trparms,'maxiter')
     trparms = settrain(trparms,'maxiter','default');
  end
  if ~isfield(trparms,'critmin')
     trparms = settrain(trparms,'critmin','default');
  end
  if ~isfield(trparms,'critterm')
     trparms = settrain(trparms,'critterm','default');
  end
  if ~isfield(trparms,'gradterm')
     trparms = settrain(trparms,'gradterm','default');
  end
  if ~isfield(trparms,'paramterm')
     trparms = settrain(trparms,'paramterm','default');
  end
  if ~isfield(trparms,'lambda')
     trparms = settrain(trparms,'lambda','default');
  end
  if ~isfield(trparms,'D')
     trparms = settrain(trparms,'D','default');
     D0 = trparms.D;
  end
end
if length(trparms.D)==1,              % Scalar weight decay parameter
   D0 = trparms.D(ones(1,parameters))';      
elseif length(trparms.D)==2,          % Two weight decay parameters
   D0 = trparms.D([ones(1,parameters2) 2*ones(1,parameters1)])';
elseif length(trparms.D)>2,           % Individual weight decay
   D0 = trparms.D;
end
D0 = D0(:);
D = D0(theta_index);
trparms = settrain(trparms,'maxiter',prparms(1));


%----------------------------------------------------------------------------------
%---------------                    MAIN LOOP                        --------------
%----------------------------------------------------------------------------------
while reduced>=minweights,   
        
  % >>>>>>>>>>>>>>>>>>>>>>>>>      Retrain Network      <<<<<<<<<<<<<<<<<<<<<<<<<<<
  % -- Don't retrain the first time --
  if ~FirstTimeFlag,
    trparms.D = D;
    [W1,W2,dummy1,dummy2,dummy3] = marq(NetDef,W1,W2,PHI,Y,trparms);
    theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
    theta_red = theta(theta_index);       % Vector containing  non-zero parameters
    theta_data(:,reduced) = theta;        % Store parameter vector
  end
      

  % >>>>>>>>>>>>>  COMPUTE NETWORK OUTPUT FROM TEST DATA y2(theta)   <<<<<<<<<<<<<<
  % -- Compute only if a test set is present --
  if TestDataFlag,
    htest1 = W1*PHI2_aug; 
    ytest1(H_hidden,:) = pmntanh(htest1(H_hidden,:));
    ytest1(L_hidden,:) = htest1(L_hidden,:);
    ytest1_aug=[ytest1;ones(1,N2)];
        
    htest2 = W2*ytest1_aug;
    ytest2(H_output,:) = pmntanh(htest2(H_output,:));
    ytest2(L_output,:) = htest2(L_output,:);

    E        = Y2 - ytest2;               % Training error
    E_vector = E(:);                      % Reshape E into a long vector
    SSE      = E_vector'*E_vector;        % Sum of squared errors (SSE)
    PI_test = SSE/(2*N2);                 % Cost function evaluated on test data
    PI_test_vec(reduced) = PI_test;       % Collect PI_test in vector
  end


  % >>>>>>>>>>>  COMPUTE NETWORK OUTPUT FROM TRAINING DATA y2(theta)   <<<<<<<<<<<<
  h1 = W1*PHI_aug;  
  y1(H_hidden,:) = pmntanh(h1(H_hidden,:));
  y1(L_hidden,:) = h1(L_hidden,:);
  y1_aug=[y1; ones(1,N)];

  h2 = W2*y1_aug;
  y2(H_output,:) = pmntanh(h2(H_output,:));
        y2(L_output,:) = h2(L_output,:);
        
  E        = Y - y2;                      % Training error
  E_vector = E(:);                        % Reshape E into a long vector
  SSE      = E_vector'*E_vector;          % Sum of squared errors (SSE)
  PI = SSE/(2*N);                         % Value of cost function
  PI_vector(reduced) = PI;                % Collect PI in vector


  % >>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX    <<<<<<<<<<<<<<<<<<<<<<<<<
  % (The derivative of each network output (y2) with respect to each weight)

  % ============   Elements corresponding to the linear output units   ============
  for i = L_output',
    index1 = (i-1) * (hidden + 1) + 1;

    % -- The part of PSI corresponding to hidden-to-output layer weights --
    PSI(index1:index1+hidden,index2+i) = y1_aug;
    % ---------------------------------------------------------------------
  
    % -- The part of PSI corresponding to input-to-hidden layer weights ---
    for j = L_hidden',
       PSI(index(j):index(j)+inputs,index2+i) = W2(i,j)*PHI_aug;
    end
      
    for j = H_hidden',
      tmp = W2(i,j)*(1-y1(j,:).*y1(j,:)); 
      PSI(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:).*PHI_aug;
    end 
    % ---------------------------------------------------------------------    
  end
  
  % ======= Elements corresponding to the hyperbolic tangent output units   =======
  for i = H_output',
    index1 = (i-1) * (hidden + 1) + 1;

    % -- The part of PSI corresponding to hidden-to-output layer weights --
    tmp = 1 - y2(i,:).*y2(i,:);
    PSI(index1:index1+hidden,index2+i) = y1_aug.*tmp(ones_h,:);
    % ---------------------------------------------------------------------
         
    % -- The part of PSI corresponding to input-to-hidden layer weights ---
    for j = L_hidden',
      tmp = W2(i,j)*(1-y2(i,:).*y2(i,:));
      PSI(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:).* PHI_aug;
    end
      
    for j = H_hidden',
      tmp  = W2(i,j)*(1-y1(j,:).*y1(j,:));
      tmp2 = (1-y2(i,:).*y2(i,:));
      PSI(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:)...
                                                .*tmp2(ones_i,:).* PHI_aug;
    end
    % ---------------------------------------------------------------------
  end
        

  % >>>>>>>>>>>>>>>>>>>>>>>>    COMPUTE THE HESSIAN MATRIX   <<<<<<<<<<<<<<<<<<<<<<
  PSI_red = PSI(theta_index,:);
  R     = PSI_red*PSI_red';
  H     = R;
  index3   = 1:(reduced+1):(reduced^2);            % A third useful vector
  H(index3) = H(index3) + D';                      % Add weight decay to diagonal
  gamma1 = (diag(R)./diag(H))'*(diag(R)./diag(H)); % Effective # of parameters
  FPE = (N+gamma1)*PI/(N-gamma1);                  % FPE estimate
  
  % --- Use the following lines for a more precise FPE estimate ---
  % H_inv = inv(H);      
  % RHinv = R*H_inv;
  % gamma1=trace(RHinv*RHinv);
  % gamma2=trace(RHinv);        
  % FPE = (N+gamma1)*PI/(N+gamma1-2*gamma2);  
  
  FPE_vector(reduced) = FPE;                       % Collect FPE estimate
  deff_vec(reduced)=gamma1;                        % Collect effective # of param.


  % >>>>>>>>>>>>    PLOT THE PI's AND THE CURRENT NETWORK STRUCTURE   <<<<<<<<<<<<<
  % --- Draw PI's ---
  figure(1);
  pvec=[reduced pvec];
  if TestDataFlag
    plot(pvec,PI_vector(pvec),'x',pvec,FPE_vector(pvec),'+',...
    pvec,PI_test_vec(pvec),'o')
    title('x = training error,   + = FPE,   o = test error')
  else
    plot(pvec,PI_vector(pvec),'x',pvec,FPE_vector(pvec),'+')
    title('x = training error,  + = FPE')
  end
  set(gca,'Xlim',[0 reduced0]);
  xlabel('Parameters');
  drawnow
     
  % --- Draw pruned network ---
  figure(2);
  drawnet(W1,W2,eps);
  title(['Network after having pruned ',int2str(pr),' weights']);
  drawnow


  % >>>>>>>>>>>>>  ELIMINATE THE WEIGHT HAVING THE SMALLEST SALIENCY  <<<<<<<<<<<<<
  nopruned = floor(max(1,reduced*RePercent/100));       % No of parms to prune
  zeta = (D+diag(R)/2).*theta_red.*theta_red;           % Saliencies
  [zeta_sorted,min_index] = sort(zeta);                 % Sort in ascending order
  theta_red(min_index(1:nopruned)) = zeros(nopruned,1); % Eliminate weights
  theta(theta_index) =theta_red;
       
  theta_index = theta_index(sort(min_index(nopruned+1:reduced)));
  theta_red = theta(theta_index);                       % Non-zero weights
  reduced  = reduced - nopruned;                        % Remaining weights
  pr = pr + nopruned;                                   % Total # of pruned weights
  theta_data(:,reduced) = theta;                        % Store parameter vector
  D = D0(theta_index);   
   
  % -- Put the parameters back into the weight matrices --
  W1 = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
  W2 = reshape(theta(1:parameters2),hidden+1,outputs)';
  FirstTimeFlag=0;
end

%----------------------------------------------------------------------------------
%-------------                END OF NETWORK PRUNING                  -------------
%----------------------------------------------------------------------------------

fprintf('\n\n\n  -->  Pruning session terminated  <--\n\n\n');
