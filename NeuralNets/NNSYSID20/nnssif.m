function [W1,W2,obsidx,PI_vector,iteration,lambda]=nnssif(NetDef,nx,W1,W2,obsidx,trparms,Y,U)
%  NNSSIF
%  ------
%          Determine a neural network model on state space innovations form:
%          x(t) = f(x(t-1),U(t-1),E(t-1))
%          y(t) = Cx(t)
%
%          The can function handle Multi-Input Multi-Outputs systems (MIMO).
%
%  CALL:
%  [W1,W2,obsidx,critvec,iteration,lambda]=...
%                      nnssif(NetDef,nx,W1,W2,obsidx,trparms,Y,U)
%
%  INPUTS:
%  U       : Input signal (= control signal) (left out in the nnarma case)
%            dim(U) = [(inputs) * (# of data)]
%  Y       : Output signal. dim(Y) = [(outputs) * (# of data)]
%  nx      : # of states (= system order)
%  NetDef  : Network definition. Remember to define 'nx' output units.
%  W1,W2   : Input-to-hidden layer and hidden-to-output layer weights.
%            dim(W1)=[(# of hidden units) * (# nx+inputs+outputs)]
%            dim(W2)=[nx * (# of hidden units)]
%            If they are passed as [], they are initialized automatically
%  obsidx  : Pseudo-observability indices. Sum must equal nx.
%            If passed as [] a particular set of indices is selected.
%  trparms:  Data structure with parameters associated with the
%            training algorithm (optional). Use the function SETTRAIN if
%            you do not want to use the default values.
% 
%  See the function MARQ for an explanation of the remaining input arguments
%  as well as of the returned variables.

%  Programmed by : Magnus Norgaard, IAU/IMM, technical University of Denmark
%  LastEditDate  : January 15, 2000

%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------
[ny,Ndat]= size(Y);                     % # of outputs and # of data
[nu,Ndat]= size(U);                     % # of inputs 
N        = Ndat - 1;                    % Size of complete training set
hidden   = length(NetDef(1,:));         % Number of hidden neurons
inputs   = nx+nu+ny;                    % Number of inputs to the network
outputs  = nx;                          % Network outputs equal to the number of states 
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
y1       = zeros(hidden,N);             % Hidden layer outputs
y1       = [y1;ones(1,N)];
y2       = zeros(outputs,N);            % Network output
Yhat     = zeros(ny,N);                 % Initialize prediction matrix
E        = zeros(ny,N);                 % Initialize prediction error matrix
E_new    = zeros(ny,N);                 % Initialize prediction error matrix
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A useful vector!
index2  = (0:N-1)*outputs;              % Yet another useful vector
iteration= 1;                           % Counter variable
dw       = 1;                           % Flag telling that the weights are new
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters=parameters1 + parameters2;   % Total # of weights
ones_h   = ones(hidden+1,1);            % A vector of ones
ones_i   = ones(inputs+1,1);            % Another vector of ones
if isempty(obsidx),                     % Choose observability indices if necessary
  rowidx   = [1:ny-1 nx];
  obsidx = [rowidx(1) rowidx(2:ny)-rowidx(1:ny-1)]; 
else                                    % Otherwise compute the row indices
  obsidx=obsidx(:)';
  rowidx=obsidx;
  for k=2:ny,
    rowidx(k)=obsidx(k)+rowidx(k-1);
  end
end
nrowidx = 1:nx;                         % Not row indices
nrowidx(rowidx)=[];
hid1=floor(hidden/2+0.5);
hid2=hidden-hid1;
if isempty(W1) | isempty(W2),           % Initialize weights if nescessary
  W1 = 0.05*(rand(hidden,inputs+1)-0.5);
  W2 = 0.05*(rand(nx,hidden+1)-0.5);                       
end
W1(hid1+1:hidden,1:nx) = zeros(hid2,nx);
W2(rowidx,hid1+1:hidden) = zeros(ny,hid2);
W2(nrowidx,1:hid1) = zeros(nx-ny,hid1);
Cidx=[1 rowidx(1:ny-1)+1];
C = zeros(ny,nx);
C(1:ny,Cidx)=eye(ny);
                                        % Parameter vector containing all weights
theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
index3   = 1:(reduced+1):(reduced^2);   % A third useful vector
PSIx    = zeros(reduced,N*nx);          % Deriv. of net output w.r.t. each weight
PSI      = zeros(reduced,N*ny);         % Deriv. of model output w.r.t. each weight
RHO      = zeros(parameters,N*ny);      % Partial -"-  -"-
txindex  = reshape(1:nx*N,nx,N);        % Index vector for PSIx
tyindex  = reshape(1:ny*N,ny,N);        % Index vector for PSI
NNy      = N*ny;

lambda_old = 0;
if ~exist('trparms') | isempty(trparms) % Default training parameters
  trparms = settrain;
  lambda  = trparms.lambda;
  skip    = trparms.skip+1;
  D       = trparms.D;
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
     D = trparms.D;
  else
    if length(trparms.D)==1,              % Scalar weight decay parameter
      D = trparms.D(ones(1,reduced))';      
    elseif length(trparms.D)==2,          % Two weight decay parameters
      D = trparms.D([ones(1,parameters2) 2*ones(1,parameters1)])';
      D = D(theta_index);
    elseif length(trparms.D)>2,           % Individual weight decay
      D = trparms.D(:);
    end
  end
end
D = D(:);
skipidx  = ny*(skip-1)+1:NNy;
N2       = N-skip+1;
critdif  = trparms.critterm+1;            % Initialize stopping variables
gradmax  = trparms.gradterm+1;
paramdif = trparms.paramterm+1;
PI_vector = zeros(trparms.maxiter,1);     % Vector for storing criterion values


% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
PHI = zeros(inputs,N);
PHI(nx+1:nx+nu,:) = U(:,1:N);
PHI_aug = [PHI;ones(1,N)];                % Augment PHI with a row containing ones
Y       = Y(:,2:Ndat);                    % Extract the 'target' part of Y

%----------------------------------------------------------------------------------
%--------------                   TRAIN NETWORK                       -------------
%----------------------------------------------------------------------------------
clc;
c=fix(clock);
fprintf('Network training started at %2i.%2i.%2i\n\n',c(4),c(5),c(6));


% >>>>>>>>>>>>>>>>>>>>>  COMPUTE NETWORK OUTPUT  y2(theta)   <<<<<<<<<<<<<<<<<<<<<<
for t=1:N,
  h1 = W1*PHI_aug(:,t);                 % Hidden neuron outputs
  y1(H_hidden,t) = pmntanh(h1(H_hidden));
  y1(L_hidden,t) = h1(L_hidden);    

  h2 = W2*y1(:,t);                      % Predicted states
  y2(H_output,t) = pmntanh(h2(H_output,:));
  y2(L_output,t) = h2(L_output,:);
  y2(nrowidx,t)  = y2(nrowidx,t) + PHI_aug(nrowidx+1,t);
  Yhat(:,t)      = C*y2(:,t);

  E(:,t) = Y(:,t) - Yhat(:,t);          % Prediction error
  for d=1:min(1,N-t),
    PHI_aug(1:nx,t+1) = y2(:,t);
    PHI_aug(nx+nu+1:inputs,t+1) = E(:,t);
  end
end
Evec   = E(:);
SSE=Evec(skipidx)'*Evec(skipidx);                  % Sum of squared errors (SSE)
PI       = (SSE+theta_red'*(D.*theta_red))/(2*N2); % Performance index

% Iterate until stopping criterion is satisfied
while (iteration<=trparms.maxiter & PI>trparms.critmin & lambda<1e7 & ...
       (critdif>trparms.critterm | gradmax>trparms.gradterm | ...
       paramdif>trparms.paramterm))
if dw==1,
% >>>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE RHO MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<<
% Partial derivative of output (y2) with respect to each weight and neglecting
% that the model inputs (the residuals) depends on the weights

    % ==========   Elements corresponding to the linear output units   ============
    for i = L_output'
      index1 = (i-1) * (hidden + 1) + 1;

      % -- The part of RHO corresponding to hidden-to-output layer weights --
      RHO(index1:index1+hidden,index2+i) = y1;
      % ---------------------------------------------------------------------
 
      % -- The part of RHO corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        RHO(index(j):index(j)+inputs,index2+i) = W2(i,j)*PHI_aug;
      end
     
      for j = H_hidden',
        tmp = W2(i,j)*(1-y1(j,:).*y1(j,:)); 
        RHO(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:).*PHI_aug;
      end 
      % ---------------------------------------------------------------------    
    end

    % ============  Elements corresponding to the tanh output units   =============
    for i = H_output',
      index1 = (i-1) * (hidden + 1) + 1;

      % -- The part of RHO corresponding to hidden-to-output layer weights --
      tmp = 1 - y2(i,:).*y2(i,:);
      RHO(index1:index1+hidden,index2+i) = y1.*tmp(ones_h,:);
      % ---------------------------------------------------------------------
         
      % -- The part of RHO corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        tmp = W2(i,j)*(1-y2(i,:).*y2(i,:));
        RHO(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:).* PHI_aug;
      end
      
      for j = H_hidden',
        tmp  = W2(i,j)*(1-y1(j,:).*y1(j,:));
        tmp2 = (1-y2(i,:).*y2(i,:));
        RHO(index(j):index(j)+inputs,index2+i) = tmp(ones_i,:)...
                                                  .*tmp2(ones_i,:).* PHI_aug;
      end
      % ---------------------------------------------------------------------
    end
    RHO_red = RHO(theta_index(1:reduced),:);


% >>>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<<
    % ---------- Find derivative of output wrt. the past residuals ----------
    for t=1:N,
      dxdy1 = W2(:,1:hidden);
      for j = H_output',
        dxdy1(j,:) = W2(j,1:hidden)*(1-y2(j,t).*y2(j,t));
      end

      % Matrix with partial derivatives of the output from each hidden neurons with
      % respect to each input:
      dy1dx = W1(:,1:nx);
      dy1de = W1(:,nx+nu+1:inputs);
      for j = H_hidden',
        dy1dx(j,:) = W1(j,1:nx)*(1-y1(j,t).*y1(j,t));
        dy1de(j,:) = W1(j,nx+nu+1:inputs)*(1-y1(j,t).*y1(j,t));
      end

      % Matrix with partial derivative of each output with respect to each input
      Ahat = dxdy1 * dy1dx;
      Ahat(nrowidx,nrowidx+1) = eye(nx-ny);
      Khat = dxdy1 * dy1de;


    % ---------- Determine PSI by "filtering" ----------
      PSIx(:,txindex(:,t))=RHO_red(:,txindex(:,t));
      for t1=1:min(1,t-1),
        PSIx(:,txindex(:,t))  = PSIx(:,txindex(:,t)) + PSIx(:,txindex(:,t-t1))...
                                                                 * (Ahat-Khat*C)';
      end
      PSI(:,tyindex(:,t)) = PSIx(:,txindex(:,t))*C';
    end

   
% >>>>>>>>>>>>>>>>>>>>>>>>>>>        COMPUTE h_k        <<<<<<<<<<<<<<<<<<<<<<<<<<<
    % -- Gradient --
    G = PSI(:,skipidx)*Evec(skipidx)-D.*theta_red;
    
    % -- Hessian  --
    H = PSI(:,skipidx)*PSI(:,skipidx)';
    H(index3) = H(index3)'+D;                       % Add diagonal matrix  
    dw = 0;
  end
   
  % -- L-M Hessian  --
  H(index3) = H(index3)'+(lambda-lambda_old);       % Add diagonal matrix

  % -- Search direction --
  h = H\G;                                          % Solve for search direction

  % -- Compute 'apriori' iterate --
  theta_red_new = theta_red + h;                    % Update parameter vector
  theta(theta_index) = theta_red_new;

  % -- Put the parameters back into the weight matrices --
  W1_new = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
  W2_new = reshape(theta(1:parameters2),hidden+1,outputs)';

    
    
% >>>>>>>>>>>>>>>>>>>>   COMPUTE NETWORK OUTPUT  y2(theta+h)   <<<<<<<<<<<<<<<<<<<<
for t=1:N,
  h1 = W1_new*PHI_aug(:,t);                 % Hidden neuron outputs
  y1(H_hidden,t) = pmntanh(h1(H_hidden));
  y1(L_hidden,t) = h1(L_hidden);    

  h2 = W2_new*y1(:,t);                      % Predicted states
  y2(H_output,t) = pmntanh(h2(H_output,:));
  y2(L_output,t) = h2(L_output,:);
  y2(nrowidx,t)  = y2(nrowidx,t) + PHI_aug(nrowidx+1,t);
  Yhat(:,t)      = C*y2(:,t);

  E_new(:,t) = Y(:,t) - Yhat(:,t);          % Prediction error
  Evec_new   = E_new(:);
  for d=1:min(1,N-t),
    PHI_aug(1:nx,t+1) = y2(:,t);
    PHI_aug(nx+nu+1:inputs,t+1) = E_new(:,t);
  end
end
SSE_new=Evec_new(skipidx)'*Evec_new(skipidx);     % Sum of squared errors (SSE)
PI_new   = (SSE_new + theta_red_new'*(D.*theta_red_new))/(2*N2); % PI
    

% >>>>>>>>>>>>>>>>>>>>>>>>>>>       UPDATE  lambda     <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  L = h'*G + h'*(h.*(D+lambda));
  lambda_old = lambda;

  % Decrease lambda if SSE has fallen 'sufficiently'
  if 2*N2*(PI - PI_new) > (0.75*L),
    lambda = lambda/2;
  
  % Increase lambda if SSE has grown 'sufficiently'
  elseif 2*N2*(PI-PI_new) <= (0.25*L),
    lambda = 2*lambda;
  end


% >>>>>>>>>>>>>>>>>>>>       UPDATES FOR NEXT ITERATION        <<<<<<<<<<<<<<<<<<<<
  % Update only if criterion has decreased
  if PI_new < PI,                      
    critdif  = PI-PI_new;                           % Criterion difference
    gradmax  = max(abs(G))/N2;                      % Maximum gradient
    paramdif = max(abs(theta_red_new - theta_red)); % Maximum parameter dif.
    W1 = W1_new;
    W2 = W2_new;
    theta_red = theta_red_new;
    E = E_new;
    Evec = Evec_new;
    PI = PI_new;
    dw = 1;
    lambda_old = 0;
    iteration = iteration + 1;
    PI_vector(iteration-1) = PI;                             % Collect PI in vector
    switch(trparms.infolevel)                                % Print on-line inform
       case 1
          fprintf('# %i   W=%4.3e  critdif=%3.2e  maxgrad=%3.2e  paramdif=%3.2e\n',...
                                                  iteration-1,PI,critdif,gradmax,paramdif);
       otherwise
          fprintf('iteration # %i   W = %4.3e\r',iteration-1,PI);
    end
  end
end
%----------------------------------------------------------------------------------
%--------------              END OF NETWORK TRAINING                  -------------
%----------------------------------------------------------------------------------
iteration = iteration-1;
PI_vector = PI_vector(1:iteration);
c=fix(clock);
fprintf('\n\nNetwork training ended at %2i.%2i.%2i\n',c(4),c(5),c(6));

