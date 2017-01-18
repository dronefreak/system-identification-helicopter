function [theta_data,PI_vector,FPE_vector,PI_test_vec,deff_vec,pvec]=...
             nnprune(method,NetDef,W1,W2,U,Y,NN,trparms,prparms,U2,Y2,Chat)
%  NNPRUNE
%  -------
%         This function applies the Optimal Brain Surgeon (OBS) strategy for
%         pruning neural network models of dynamic systems. That is networks
%         trained by NNARX, NNOE, NNARMAX1, NNARMAX2, or their recursive
%         counterparts.
%
%
%  CALL:
% [theta_data,NSSEvec,FPEvec,NSSEtestvec,deff,pvec]=...
%          nnprune(method,NetDef,W1,W2,U,Y,NN,trparms,prparms,U2,Y2,Chat)
%
%  INPUT:
%  method             : The function applied for generating the model. For
%                       example method='nnarx' or method='nnoe'
%  NetDef, W1, W2,
%  U, Y, trparms      : See for example the function MARQ
%  U2,Y2              : Test data. This can be used for pointing out
%                       the optimal network architecture is achieved. Pass
%                       two []'s if a test set is not available.
%  Chat               : See NNARMAX1
%  prparms (optional) : Parameters associated with the pruning session
%                       prparms = [iter RePercent]
%                       iter      : Max. number of retraining iterations
%                       RePercent : Prune 'RePercent' percent of the
%                                   remaining weights (0 = prune one at a time)
%                       If passed as [], prparms=[50 0] will be used.
%  
%  OUTPUT:
%  theta_data  : Matrix containing the parameter vectors saved after each
%                weight elimination round.
%  NSSEvec     : Vector containing the training error (SSE/2N) after each
%                weight elimination.
%  FPEvec      : Contains the FPE estimate of the average generalization error
%  NSSEtestvec : Contains the normalized SSE evaluated on the test set
%  deff        : Contains the "effective" number of weights
%  pvec        : Index to the above vectors
%
%  SEE ALSO: OBSPRUNE and OBDPRUNE on how to prune ordinary feedforward
%            networks. See also the function NETSTRUC on how to extract
%            the weight matrices from the matrix theta_data (notice that for
%            NNARMAX1-models one must remove the bottom deg(C) rows first).
 
%  Written by : Magnus Norgaard, IAU/IMM, Technical Univ. of Denmark
%  LastEditDate  : Jan. 6, 2000

%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------
more off
if ~isempty(Y2), TestDataFlag = 1;      % Check if test data was given as argument
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
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters = parameters1 + parameters2; % Total # of weights
                                        % Parameter vector containing all weights
if strcmp(method,'nnarmax1') | strcmp(method,'nnrarmx1'),
  mflag=2;
  parameters12 = parameters;
  nc = length(Chat)-1;
  parameters = parameters12+nc;
  theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1) ; Chat(2:nc+1)'];
else
  theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
end
if nargin<12, Chat=[]; end;
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
reduced0 = reduced;                     % Copy of 'reduced'. Will be constant
theta_data=zeros(parameters,parameters);% Matrix used for collecting theta vectors
theta_data(:,reduced) = theta;          % Insert 'initial' theta
p0       = 1e6;                         % Diag. element of H_inv (no weight decay)
H_inv    = p0*eye(reduced);             % Initial inverse Hessian (no weight decay)
Ident    = eye(outputs);                % Identity matrix
PI_vector= zeros(1,reduced);            % A vector containing the collected PI's
FPE_vector= zeros(1,reduced);           % Vector used for collecting FPE estimates
deff_vec = zeros(1,reduced);            % The effective number of parameters
minweights = 2;                         % Prune until 'minweights' weights remain
FirstTimeFlag=1;                        % Initialize flag
pr = 0;                                 % Initialize counter
pvec=[];                                % Initialize index vector
HiddenIndex = ones((hidden+1),1);       % Connection to hidden no.
for k=1:hidden,
  HiddenIndex = [HiddenIndex;k*ones(inputs+1,1)];
end
ConnectToHidden = (inputs+1)*ones(hidden,1); % Connections to each hidden unit
if ~exist('trparms') | isempty(trparms) % Default training parameters
  trparms = settrain;
  lambda  = trparms.lambda;
  skip=trparms.skip+1;
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
  lambda    = trparms.lambda;
  if ~isfield(trparms,'skip')
     trparms= settrain(trparms,'skip','default');
  end
  skip=trparms.skip+1;
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
D0=D0(:);
D = D0(theta_index);


% ---------- NNARX model ----------
if strcmp(method,'nnarx') | strcmp(method,'nnrarx'),
  mflag=1;
    if length(NN)==1                      % nnar model
    nb = 0;
    nk = 0;
    nu = 0;
  else                                  % nnarx or nnoe model
    [nu,N] = size(U); 
    nb = NN(2:1+nu); 
    nk = NN(2+nu:1+2*nu);
  end
  nc = 0;
  
% --------- NNARMAX1 model --------
elseif strcmp(method,'nnarmax1') | strcmp(method,'nnrarmx1'),
  mflag=2;

% --------- NNARMAX2 model --------
elseif strcmp(method,'nnarmax2') | strcmp(method,'nnrarmx2'),
  mflag=3;

% --------- NNOE model --------
elseif strcmp(method,'nnoe'),
 mflag=4;
  if length(NN)==1                      % nnar model
    nb = 0;
    nk = 0;
    nu = 0;
  else                                  % nnarx or nnoe model
    [nu,N] = size(U); 
    nb = NN(2:1+nu); 
    nk = NN(2+nu:1+2*nu);
  end
  nc = 0;
else
  disp('Unknown method!!!!!!!!');
  break
end

if mflag==2 | mflag==3,
  if length(NN)==2                      % nnarma model
    nc     = NN(2);
    nb     = 0;
    nk     = 0;
    nu     = 0;
  else                                  % nnarmax model
    [nu,Ndat]= size(U); 
    nb     = NN(2:1+nu);
    nc     = NN(2+nu);
    nk     = NN(2+nu+1:2+2*nu);
  end
end

% --------- Common initializations --------
Ndat     = length(Y);                   % # of data
na = NN(1);
nab = na+sum(nb);
nabc = nab+nc;
nmax     = max([na,nb+nk-1,nc]);        % 'Oldest' signal used as input to the model
N        = Ndat - nmax;                 % Size of training set
NN2      = N-skip+1;
if TestDataFlag,                        % Initializations if a test set exists
  Ndat2 = length(Y2);                   % Total # of data in test set
  N2    = Ndat2 - nmax;                 % Size of test set
  N2tot = N2 - skip+1;
  ytest1 = zeros(hidden,N2);            % Hidden layer outputs 
  ytest1 = [ytest1;ones(1,N2)];         % Hidden layer outputs 
  ytest2 = zeros(outputs,N2);           % Network output

  %------  CONSTRUCT THE REGRESSION MATRIX PHI ------
  if mflag~=3,
    PHI2 = zeros(nab,N2);
  else
    PHI2 = zeros(nabc,N2);
  end
  jj  = nmax+1:Ndat2;
  for k = 1:na, PHI2(k,:)    = Y2(jj-k); end
  index4 = na;
  for kk = 1:nu,
    for k = 1:nb(kk), PHI2(k+index4,:) = U2(kk,jj-k-nk(kk)+1); end
    index4 = index4 + nb(kk);
  end
  PHI2_aug    = [PHI2;ones(1,N2)];      % Augment PHI with a row containing ones
  Y2    = Y2(nmax+1:Ndat2);
  PI_test_vec = zeros(1,reduced);       % Collected PI's for the test set
end
trparmsp = settrain(trparms,'maxiter',iter);

%----------------------------------------------------------------------------------
%---------------                    MAIN LOOP                        --------------
%----------------------------------------------------------------------------------
while reduced>=minweights,   
  % >>>>>>>>>>>>>>>>>>>>>>>>>      Retrain Network      <<<<<<<<<<<<<<<<<<<<<<<<<<< 
  % -- Don't retrain the first time --
  
  if ~FirstTimeFlag,
    trparmsp.D = D;
    if mflag==1,
      [W1,W2,dummy1,dummy2,dummy3] = nnarx(NetDef,NN,W1,W2,trparmsp,Y,U);
    elseif mflag==2,
      [W1,W2,Chat] = nnarmax1(NetDef,NN,W1,W2,Chat,trparmsp,Y,U);
    elseif mflag==3,
      [W1,W2,dummy1,dummy2,dummy3] = nnarmax2(NetDef,NN,W1,W2,trparmsp,Y,U);
    elseif mflag==4,
      [W1,W2,dummy1,dummy2,dummy3] = nnoe(NetDef,NN,W1,W2,trparmsp,Y,U);
    end
    if mflag==2,
      theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1) ; Chat(2:nc+1)'];
    else
      theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
    end
    theta_red = theta(theta_index);       % Vector containing non-zero parameters
    if ElimWeights==1,                    % Store parameter vector
      theta_data(:,reduced) = theta;
    else
      theta_data(:,[reduced reduced+1]) = [theta theta];
    end
  end  


  % >>>>>>>>>>>>>  COMPUTE NETWORK OUTPUT FROM TEST DATA y2(theta)   <<<<<<<<<<<<<<
  % -- Compute only if a test set is present -- 
  if TestDataFlag,
% ---------- NNARX model ----------
if mflag==1,
  htest1 = W1*PHI2_aug;  
  ytest1(H_hidden,:) = pmntanh(htest1(H_hidden,:));
  ytest1(L_hidden,:) = htest1(L_hidden,:);
    
  htest2 = W2*ytest1;
  ytest2(H_output,:) = pmntanh(htest2(H_output,:));
  ytest2(L_output,:) = htest2(L_output,:);

  E     = Y2 - ytest2;                    % Error between Y and deterministic part


% --------- NNARMAX1 model --------
elseif mflag==2,
  htest1 = W1*PHI2_aug;  
  ytest1(H_hidden,:) = pmntanh(htest1(H_hidden,:));
  ytest1(L_hidden,:) = htest1(L_hidden,:);
    
  htest2 = W2*ytest1;
  ytest2(H_output,:) = pmntanh(htest2(H_output,:));
  ytest2(L_output,:) = htest2(L_output,:);

  Ebar     = Y2 - ytest2;                    % Error between Y and deterministic part
  E        = filter(1,Chat,Ebar);            % Prediction error
  ytest2   = ytest2 - E;                       % One step ahead prediction


% --------- NNARMAX2 model --------
elseif mflag==3,
  for t=1:N2,
    htest1 = W1*PHI2_aug(:,t);  
    ytest1(H_hidden,t) = pmntanh(htest1(H_hidden));
    ytest1(L_hidden,t) = htest1(L_hidden);    

    htest2 = W2*ytest1(:,t);
    ytest2(H_output,t) = pmntanh(htest2(H_output,:));
    ytest2(L_output,t) = htest2(L_output,:);

    E(:,t) = Y2(:,t) - ytest2(:,t);          % Prediction error
    for d=1:min(nc,N2-t),
      PHI2_aug(nab+d,t+d) = E(:,t);
    end
  end


% ---------- NNOE model ----------
elseif mflag==4,
  for t=1:N2,
    htest1 = W1*PHI2_aug(:,t);;  
    ytest1(H_hidden,t) = pmntanh(htest1(H_hidden));
    ytest1(L_hidden,t) = htest1(L_hidden);    

    htest2 = W2*ytest1(:,t);
    ytest2(H_output,t) = pmntanh(htest2(H_output,:));
    ytest2(L_output,t) = htest2(L_output,:);

    for d=1:min(na,N2-t),
      PHI2_aug(d,t+d) = ytest2(:,t);
    end
  end
  E     = Y2 - ytest2;                    % Error between Y and deterministic part
end

    SSE     = E(skip:N2)*E(skip:N2)';     % Sum of squared errors (SSE)
    PI_test = SSE/(2*N2tot);              % Cost function evaluated on test data
    PI_test_vec(reduced) = PI_test;       % Collect PI_test in vector
end


  % >>>>>>>>>>>>>>>>>>>>>>  GET NETWORK OUTPUT AND GRADIENT   <<<<<<<<<<<<<<<<<<<<<<
  [PSI,E] = getgrad(method,NetDef,NN,W1,W2,Chat,Y,U);
  PI = E(skip:N)*E(skip:N)'/(2*NN2);      % Performance index
  PI_vector(reduced) = PI;                % Collect PI in vector

        
        
  % >>>>>>>>>>>>>>>>>>>>>>>>    COMPUTE THE HESSIAN MATRIX   <<<<<<<<<<<<<<<<<<<<<<
  PSI_red = PSI(theta_index,skip:N);

  % --- Inverse Hessian if no weight decay ---
  if D==0,
    H_inv = p0*eye(reduced);
    for k=1:outputs:NN2,                       % Iterative solution
      psi=PSI_red(:,(k-1)*outputs+1:k*outputs);
      H_inv = H_inv - H_inv*psi*inv(Ident + psi'*H_inv*psi)*psi'*H_inv;
    end
    FPE = (NN2+reduced)*PI/(NN2-reduced);     % FPE estimate
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
    FPE = (NN2+gamma1)*PI/(NN2+gamma1-2*gamma2);% FPE estimate
  end 
  FPE_vector(reduced) = FPE;                  % Collect FPE estimate
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
    if D==0,
      gamma = theta_red./diag(H_inv);             % Lagrange multipliers
      zeta = theta_red.*gamma;                    % Salincies if D=0
    else
      gamma = theta_red./diag(H_inv);             % Lagrange multipliers
      HRH   = H_inv*R*H_inv;
                                                  % Saliencies
      zeta  = gamma.*(H_inv*(D.*theta_red))+gamma.*gamma.*diag(HRH)/2; 
    end
    Critical=[]; WhileFlag=1; ElimWeights=1;
    while WhileFlag,
      WhileFlag = 0;
      HiddenLeft = hidden-length(find(ConnectToHidden==0)); % Hidden units;
      [zeta_min,min_index] = min(zeta(HiddenLeft+1:reduced));% Find smallest saliency
      min_index = min_index+HiddenLeft;
      HiddenNo = HiddenIndex(theta_index(min_index)); % Connection to hidden no?
      if theta_index(min_index)~=(hidden+1),  % Not the bias
        if ConnectToHidden(HiddenNo)==1
          HiddenNoRed = find(theta_index(1:HiddenLeft)==HiddenNo);
          Eidx = [HiddenNoRed;min_index];
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
          Critical = [Critical;min_index];
          WhileFlag = 1;
        end
      end
    end
    
    if theta_index(min_index)~=(hidden+1),
      ConnectToHidden(HiddenNo) = ConnectToHidden(HiddenNo) - 1;
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
      
     % ----- Eliminate two weights -----
    elseif ElimWeights==2,
      theta_red = theta_red - H_inv(:,Eidx)*gamma_new;
      theta_red(Eidx) = [0;0];                      % Non-zero weights
      theta(theta_index) =theta_red;
      tmpindex = [1:Eidx(1)-1 Eidx(1)+1:Eidx(2)-1 Eidx(2)+1:length(theta_index)];
      theta_index = theta_index(tmpindex);
      theta_red = theta(theta_index);
      reduced  = reduced-2;                         % Remaining weights
      theta_data(:,[reduced reduced+1]) = [theta theta]; % Store parameter vector
      D=D0(theta_index);
    
      % --- Update inverse Hessian ---
      H_inv = H_inv(tmpindex,tmpindex)-H_inv(tmpindex,Eidx)...
               *inv(H_inv(Eidx,Eidx))*H_inv(Eidx,tmpindex);
      if D~=0, R = R(tmpindex,tmpindex); end
      bpr = bpr + 1;
    end
    bpr = bpr + 1;
  end    
  pr = pr + bpr-1;                             % Total # of pruned weights
  
 
  % -- Put the parameters back into the weight matrices --  
  if mflag==2,
    W1 = reshape(theta(parameters2+1:parameters12),inputs+1,hidden)';
    W2 = reshape(theta(1:parameters2),hidden+1,outputs)';
    Chat   = [1 theta(parameters12+1:parameters)'];
  else
    W1 = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
    W2 = reshape(theta(1:parameters2),hidden+1,outputs)';
  end
  FirstTimeFlag=0;
end
%----------------------------------------------------------------------------------
%-------------                END OF NETWORK PRUNING                  -------------
%----------------------------------------------------------------------------------
fprintf('\n\n\n  -->  Pruning session terminated  <--\n\n\n');
