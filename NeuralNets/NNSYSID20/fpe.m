function [FPE,deff,varest,H] = fpe(NetDef,W1,W2,PHI,Y,trparms)
% 
%  FPE
%  --- 
%           This function calculates Akaike's final prediction error
%           estimate of the average generalization error.
%
%  [FPE,deff,varest,H] = fpe(NetDef,W1,W2,PHI,Y,trparms) produces the
%  final prediction error estimate (fpe), the effective number of
%  weights in the network if the network has been trained with
%  weight decay, an estimate of the noise variance, and the Gauss-Newton
%  Hessian.
%  
%  INPUT:
%           See for example the function MARQ 
%  
%  OUTPUT:
%  FPE    : The Final prediction error estimate 
%  deff   : The effective number of weights
%  varest : Estimate of the noise variance
%  H      : The Gauss-Newton Hessian
%
%  REFERENCE:
%       J. Larsen & L.K. Hansen:
%       "Generalization Performance of Regularized Neural Network Models"
%        Proc. of the IEEE Workshop on Neural networks for Signal Proc. IV,
%        Piscataway, New Jersey, pp.42-51, 1994
%
%  SEE ALSO:  NNFPE, LOO

%  Programmed by : Magnus Norgaard, IAU/IMM, Technical Univ. of Denmark
%  LastEditDate  : Jan 5, 2000


%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------
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
theta_data=zeros(parameters,parameters);% Matrix used for collecting theta vectors
theta_data(:,reduced) = theta;          % Insert 'initial' theta
PSI      = zeros(parameters,outputs*N); % Deriv. of each output w.r.t. each weight
if nargin<6 | isempty(trparms) % Default training parameters
  trparms = settrain;
  D       = trparms.D;
else                                    % User specified values
  if ~isstruct(trparms),
     error('''trparms'' must be a structure variable.');
  end
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
 
  % --- Calculate the HEssian matrix ---
  PSI_red = PSI(theta_index,:);
  R     = PSI_red*PSI_red';
  H     = R;
  index3   = 1:(reduced+1):(reduced^2);       % A third useful vector
  H(index3) = H(index3) + D';                 % Add weight deacy to diagonal

  % --- FPE in case of no weight decay ---
if D==0,
  FPE  = PI*(N + reduced) / (N - reduced);
  deff = reduced;
  varest = 2*N*PI/(N-reduced);
else

  % --- FPE in case of weight decay ---
  H_inv  = inv(H);                            % Inverse Hessian
  RHinv  = R*H_inv;
  Dmat   = diag(D);
  gamma1 = trace(RHinv*RHinv);                % Effective # of parameters
  gamma2 = trace(RHinv);        
  gamma3 = theta(theta_index)'*Dmat*H_inv*RHinv*Dmat*theta(theta_index)/N;
  varest = (2*N*PI-N*gamma3) / (N + gamma1 - 2*gamma2);
  FPE    = (varest*(1+gamma1/N) + gamma3)/2;  % FPE estimate
  deff = gamma1;                              % Effective # of parameters
end


