function [theta_data,PI_vector,FPE_vector,PI_test_vec,deff_vec,pvec]=...
                             obsprune(NetDef,W1,W2,PHI,Y,trparms,prparms,PHI2,Y2)
%  OBSPRUNE
%  --------
%           This function applies the Optimal Brain Surgeon (OBS) algorithm
%           for pruning ordinary feedforward neural networks
%
%  CALL:
%   [theta_data,NSSEvec,FPEvec,NSSEtestvec,deff,pvec]=...
%                    obsprune(NetDef,W1,W2,PHI,Y,trparms,prparms,PHI2,Y2)
%
%  INPUT:
%  NetDef, W1, W2,
%  PHI, Y, trparms    : See for example the function MARQ
%  PHI2,Y2 (optional) : Test data. Can be used for pointing out the
%                       optimal network architecture. 
%  prparms            : Parameters assocoated with the pruning session
%                       prparms = [iter RePercent]
%                       iter      : Max. number of retraining iterations
%                       RePercent : Prune 'RePercent' percent of the
%                                   remaining weights (0 = prune one at a time)
%                       If passed as [], prparms=[50 0] will be used.
%  
%  OUTPUT:
%  theta_data  : Matrix containing all the parameter vectors
%  NSSEvec     : Vector containing the training error (SSE/2N) after each
%                weight elimination
%  FPEvec      : Contains the FPE estimate of the average generalization error
%  NSSEtestvec : Contains the test error
%  deff        : Contains the "effective" number of parameters
%  pvec        : Index to the above vectors
 
%  Programmed by : Magnus Norgaard, IAU/IMM, Technical University of Denmark
%  LastEditDate  : Jan. 8, 2000


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
PHI_aug  = [PHI;ones(1,N)];             % Augment PHI with a row containing ones
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
theta_data=zeros(parameters,reduced);   % Matrix used for collecting theta vectors
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
minweights = 2*outputs;                 % Prune until 'minweights'(>=2) weights remain
FirstTimeFlag=1;                        % Initialize flag
pr = 0;                                 % Initialize counter
pvec=[];                                % Initialize index vector
HiddenIndex = [];                       % Connection to hidden number X
for k=1:outputs,
  HiddenIndex = [HiddenIndex;(1:(hidden+1))'];
end
for k=1:hidden,
  HiddenIndex = [HiddenIndex;k*ones(inputs+1,1)];
end
ConnectToHidden = (inputs+1)*ones(hidden,1); % Connections to each hidden unit
ConnectFromHidden = outputs*ones(hidden,1);  % Connections from each hidden unit
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
    if ElimWeights==1,                    % Store parameter vector
      theta_data(:,reduced) = theta;
    else
      theta_data(:,[reduced reduced+LEidx-1]) = theta(:,ones(1,LEidx));
    end
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


  % >>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<  
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

  % --- Inverse Hessian if no weight decay ---
  if D==0,
    H_inv = p0*eye(reduced);
    for k=1:outputs:N,                        % Iterative solution
      psi=PSI_red(:,(k-1)*outputs+1:k*outputs);
      H_inv = H_inv - H_inv*psi*inv(Ident + psi'*H_inv*psi)*psi'*H_inv;
    end
    FPE = (N+reduced)*PI/(N-reduced);         % FPE estimate
    gamma1 = reduced;
    
  % --- Inverse Hessian if weight decay is being used ---
  else
    R     = PSI_red*PSI_red';
    H     = R;
    index3   = 1:(reduced+1):(reduced^2);       % A third useful vector
    H(index3) = H(index3) + D';                 % Add weight deacy to diagonal
    H_inv = inv(H);                             % Inverse Hessian
    RHinv = R*H_inv;
    gamma1=trace(RHinv*RHinv);                  % Effective # of parameters
    gamma2=trace(RHinv);        
    FPE = (N+gamma1)*PI/(N+gamma1-2*gamma2);    % FPE estimate
  end
  FPE_vector(reduced) = FPE;                    % Collect FPE estimate
  deff_vec(reduced)=gamma1;                     % Collect effective # of param.


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
  figure(2); drawnow


  % >>>>>>>>>>>>>  ELIMINATE THE WEIGHT HAVING THE SMALLEST SALIENCY  <<<<<<<<<<<<<
  nopruned = floor(max(1,reduced*RePercent/100)); % No of parms to prune
  if reduced<=minweights, break; end
  bpr = 1;
  while bpr<=nopruned
  
    % ----- Calculate all saliences -----
    if D==0,
      gamma = theta_red./diag(H_inv);             % Lagrange multipliers
      zeta = theta_red.*gamma;                    % Salincies if D=0
    else
      gamma = theta_red./diag(H_inv);             % Lagrange multipliers
      HRH   = H_inv*R*H_inv;
                                                  % Saliencies
      zeta  = gamma.*(H_inv*(D.*theta_red))+gamma.*gamma.*diag(HRH)/2; 
    end
    
  
   % ----- Add a large number to uninteresting saliences -----
   index5 = find(ConnectFromHidden==1);
   [zeta_max,max_index] = max(zeta);              % Find largest saliency
   IndexInZetaVec=[];   
   for hidno=index5',
     % Find the weight's location in theta
     index6 = hidno:(hidden+1):(hidden+1)*outputs;
     IndexInTheta = find(theta(index6)~=0);
     IndexInTheta = index6(IndexInTheta);
     
     % Find the location in zeta
     IndexInZeta = find(theta_index == IndexInTheta);
     zeta(IndexInZeta) = 2*abs(zeta_max);
     IndexInZetaVec=[IndexInZetaVec;IndexInZeta];
   end
   
    Critical=[]; WhileFlag=1; ElimWeights=1;
    while WhileFlag,
      WhileFlag = 0;
      [zeta_min,min_index] = min(zeta);               % Find smallest saliency
      
      % -- Check if weight in question belongs to W1 or W2 --
      if (theta_index(min_index)>=1 & theta_index(min_index)<=(hidden+1)*outputs),
        W=2;     % Belongs to W2
      else
        W=1;     % Belongs to W1
      end

      HiddenNo = HiddenIndex(theta_index(min_index)); % Connection to hidden no?

      % -- If not a bias for an output unit --
      if ~( any([hidden+1:hidden+1:outputs*(hidden+1)]==theta_index(min_index)) ),
      
        % -- If weight has only one connection leading to it --
        if W==1 & ConnectToHidden(HiddenNo)==1,
          index6 = HiddenNo:(hidden+1):(hidden+1)*outputs;
          HiddenNoRed = [];
          for k=index6,
            HiddenNoRed = [HiddenNoRed,find(theta_index==k)];% Index in zeta vector
          end
          Eidx = [HiddenNoRed,min_index];     % Indices to weights to and from unit
          gamma_new = inv(H_inv(Eidx,Eidx))*theta_red(Eidx);
          if ~isempty(Critical)
            if find(Critical==min_index);
              ElimWeights=2;
              break;
            end
	  end
          if D==0,
            zeta(min_index) = gamma_new'*theta_red(Eidx);
          else
            zeta(min_index) = gamma_new'*H_inv(Eidx,:)*(D.*theta_red)...
                              + gamma_new'*HRH(Eidx,Eidx)*gamma_new/2;
          end
          [zeta_max,max_index] = max(zeta);              % Find largest saliency	  
          zeta(IndexInZetaVec) = 2*abs(zeta_max);	  
          Critical = [Critical;min_index];
          WhileFlag = 1;
        end
      end      
    end
    
    % If not a bias for an output unit
    if ~(any([hidden+1:hidden+1:outputs*(hidden+1)]==theta_index(min_index))),
    
      % -- If a W2 weight is eliminated --
      if W==2,
        ConnectFromHidden(HiddenNo) = ConnectFromHidden(HiddenNo) - 1;
        
      % --if a W1 weight is eliminated -- 
      else
        ConnectToHidden(HiddenNo) = ConnectToHidden(HiddenNo) - 1;
      end
      
      if ElimWeights==2,
        ConnectFromHidden(HiddenNo) = 0;
      end
    end
    
    % ----- Eliminate one weight -----
    if ElimWeights==1, 
      theta_red = theta_red - gamma(min_index)*H_inv(:,min_index);
      theta_red(min_index) = 0;                     % Non-zero weights
      theta(theta_index) = theta_red;
      tmpindex = [1:min_index-1 min_index+1:length(theta_index)];
      theta_index = theta_index(tmpindex);
      theta_red = theta(theta_index);
      reduced  = reduced-1;                         % Remaining weights
      theta_data(:,reduced) = theta;                % Store parameter vector
      D=D0(theta_index);
    
      % --- Update inverse Hessian ---
      H_inv = H_inv(tmpindex,tmpindex)-H_inv(tmpindex,min_index)...
               *H_inv(min_index,tmpindex)/H_inv(min_index,min_index);
      if D~=0, R = R(tmpindex,tmpindex); end
      
     % ----- Eliminate more than one weight -----
    elseif ElimWeights==2,
      LEidx = length(Eidx);                  % # of weights to be pruned
      theta_red = theta_red - H_inv(:,Eidx)*gamma_new;
      theta_red(Eidx) = zeros(LEidx,1);      % Non-zero weights
      theta(theta_index) = theta_red;
      tmpindex = [1:Eidx(1)-1];
      for k=2:LEidx,
        tmpindex = [tmpindex Eidx(k-1)+1:Eidx(k)-1];
      end
      tmpindex = [tmpindex Eidx(LEidx)+1:length(theta_index)];
      theta_index = theta_index(tmpindex);
      theta_red = theta(theta_index);
      reduced  = reduced-length(Eidx);                   % Remaining weights
      theta_data(:,[reduced reduced+LEidx-1]) = theta(:,ones(1,LEidx)); % Store theta
      D = D0(tmpindex);
    
      % --- Update inverse Hessian ---
      H_inv = H_inv(tmpindex,tmpindex)-H_inv(tmpindex,Eidx)...
               *inv(H_inv(Eidx,Eidx))*H_inv(Eidx,tmpindex);
      if D~=0, R = R(tmpindex,tmpindex); end
      bpr = bpr + LEidx-1;
    end
    bpr = bpr + 1;
  end    
  pr = pr + bpr-1;                             % Total # of pruned weights
  
  % -- Put the parameters back into the weight matrices --
  W1 = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
  W2 = reshape(theta(1:parameters2),hidden+1,outputs)';
  FirstTimeFlag=0;
end
%----------------------------------------------------------------------------------
%-------------                END OF NETWORK PRUNING                  -------------
%----------------------------------------------------------------------------------
fprintf('\n\n\n  -->  Pruning session terminated  <--\n\n\n');
