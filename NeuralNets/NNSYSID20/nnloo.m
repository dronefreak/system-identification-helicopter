function [Eloo] = nnloo(NetDef,W1,W2,U,Y,NN,trparms)
%  NNLOO
%  ----- 
%  Leave-one-out estimate of the average generalization error for NNARX models.
%
%  The leave-one-out cross-validation scheme is a method for estimating
%  the average generalization error. When calling
%  Eloo = nnloo(NetDef,W1,W2,U,Y,NN,trparms) with trparms.maxiter>0, the network
%  will be retrained a maximum of trparms.maxiter iterations for each input-output
%  pair in the data set, starting from the initial weights (W1,W2). If
%  trparms.maxiter=0 an approximation to the loo-estimate based on "linear 
%  unlearning" is produced. This is in general less accurate, but it is much
%  faster to compute.
%
%  Unless trparms.maxiter=0 it is recommended that trparms.maxiter=20-40
%
%  INPUT:
%           See the function 'NNARX'
%           If trparms=[], trparms.maxiter will be set to 0.
%  
%  OUTPUT:
%  Eloo   : The leave-one-out estimate of the average generalization error 
%
%  SEE the function "nnfpe" for the final prediction error estimate
%
%  REFERENCE:
%             L.K. Hansen and J. Larsen (1995):
%             "Linear Unlearning for Cross-Validation"
%             Submitted for Advances in Computational Mathematics, April 1995
% 
%  Programmed by : Magnus Norgaard, IAU/IMM, Technical University of Denmark 
%  LastEditDate  : Jan. 15, 2000


%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------
[outputs,N] = size(Y);                  % # of outputs and # of data
[hidden,inputs] = size(W1);             % # of hidden units
inputs=inputs-1;                        % # of inputs
na = NN(1);
if length(NN)==1
  nb = 0;                               % nnar model
  nk = 0;
  nu = 0; 
else
  [nu,N]      = size(U); 
  nb = NN(2:1+nu);                      % nnarx model
  nk = NN(2+nu:1+2*nu);
end
nmax        = max([na,nb+nk-1]);        % Maximum delay
nab         = na+sum(nb);

PHI_aug = [zeros(nab,N-nmax);ones(1,N-nmax)]; % Construct the regression matrix PHI
jj  = nmax+1:N;
for k = 1:na, PHI_aug(k,:)    = Y(jj-k); end
index = na;
for kk = 1:nu,
  for k = 1:nb(kk), PHI_aug(k+index,:) = U(kk,jj-k-nk(kk)+1); end
  index = index + nb(kk);
end
Y = Y(nmax+1:N);                        % Remove first outputs in Y vector     
N=N-nmax;

L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neuron
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
y1       = zeros(hidden,N);             % Hidden layer outputs
y2       = zeros(outputs,N);            % Network output
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A usefull vector!
index2 = (0:N-1)*outputs;               % Yet another usefull vector
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
PSI      = zeros(parameters,outputs*N); % Der. of each output w.r.t. each weight
p0       = 1e8;                         % Diag. element of H_inv (no weight decay)

if nargin<7 | isempty(trparms) % Default training parameters
  trparms = settrain;
  trparms = settrain(trparms,'maxiter',0);
  D = trparms.D;
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

% Determine LOO-estimate by retraining
if trparms.maxiter>0,
  E        = zeros(outputs,1);
  y1       = zeros(hidden,1);
  y2       = zeros(outputs,1);
  ELOOvector = zeros(N,1);
  for i=1:N,
    notbeta=[1:i-1 i+1:N];
    [W1x,W2x,PI_vector,iteration,lambda]=...
                 marq(NetDef,W1,W2,PHI_aug(1:nab,notbeta),Y(:,notbeta),trparms);
                 
    % ---  Compute network output ---
    h1           = W1x*PHI_aug(:,i);  
    y1(H_hidden) = pmntanh(h1(H_hidden));
    y1(L_hidden) = h1(L_hidden);
    
    h2 = W2x*[y1;1];
    y2(H_output) = pmntanh(h2(H_output));
    y2(L_output) = h2(L_output);

    E            = Y(:,i) - y2;                % Test error
    PI           = (E'*E)/(2*N);               % Sum of squared errors 
  
    ELOOvector(i)=PI;
  end
  Eloo=sum(ELOOvector);


% Determine LOO-estimate using linear unlearning
else
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
  Ident    = eye(outputs);                % Identity matrix
    H_inv = p0*eye(reduced);
    for k=1:outputs:N,                        % Iterative solution
      psi=PSI_red(:,(k-1)*outputs+1:k*outputs);
      H_inv = H_inv - (H_inv*psi*inv(Ident + psi'*H_inv*psi)*psi'*H_inv);
    end

  % --- Inverse Hessian if weight decay is being used ---
  else
    R     = PSI_red*PSI_red';
    H     = R;
    index3   = 1:(reduced+1):(reduced^2);       % A third useful vector
    H(index3) = H(index3) + D';                 % Add weight deacy to diagonal
    H_inv = inv(H);                             % Inverse Hessian
  end


  % LOO average generalization error estimate
  Eloo=0;
  for beta=1:N,
    hjh=(PSI_red(:,beta)'*H_inv*PSI_red(:,beta));
    Eloo=Eloo+E(beta)*E(beta)*(1+hjh)/(1-hjh);
  end
  Eloo=Eloo/(2*N);
end
