function [W1,W2,PI_vector,iteration]=rpe(NetDef,W1,W2,PHI,Y,trparms)
%  RPE
%  --- 
%          Train a two layer neural network with a recursive prediction error
%          algorithm ("recursive Gauss-Newton"). Also pruned (i.e., not fully
%          connected) networks can be trained.
%
%          The activation functions can either be linear or tanh. The network
%          architecture is defined by the matrix 'NetDef', which has of two
%          rows. The first row specifies the hidden layer while the second
%          specifies the output layer.
%
%          E.g.:    NetDef = ['LHHHH'
%                             'LL---']
% 
%          (L = linear, H = tanh)
%
%          A weight is pruned by setting it to zero.
%
%          The algorithm is described in:
%          L. Ljung: "System Identification - Theory for the User"
%          (Prentice-Hall, 1987)
%
%          Notice that the bias is included as the last column
%          in the weight matrices.
%
%  CALL:
%            [W1,W2,critvec,iter]=rpe(NetDef,W1,W2,PHI,Y,trparms)
%
%  INPUT:
%  NetDef: Network definition 
%  W1    : Input-to-hidden layer weights. The matrix dimension is
%          dim(W1) = [(# of hidden units) * (inputs + 1)] (the 1 is due to the bias)
%          Use [] for a random initialization.
%  W2    : hidden-to-output layer weights.
%          dim(W2) = [(outputs)  *  (# of hidden units + 1)]
%          Use [] for a random initialization.
%  PHI   : Input vector. dim(PHI) = [(inputs)  *  (# of data)]
%  Y     : Output data. dim(Y) = [(outputs)  * (# of data)]
%  trparms: Data structure with parameters associated with the
%           training algorithm (optional). Use the function SETTRAIN if
%           you do not want to use the default values. 
%
%  OUTPUT:
%  W1, W2    : Weight matrices after training
%  critvec   : Vector containing the criterion of fit after each iteration
%  iter      : # of iterations
%  
%  Programmed by : Magnus Norgaard, IAU/IMM, Technical University of Denmark 
%  LastEditDate  : January 15, 2000


%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------
[outputs,N] = size(Y);                  % # of outputs and # of data
[inputs,N] = size(PHI);                 % # of hidden units
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
hidden = length(L_hidden)+length(H_hidden);
if isempty(W1) | isempty(W2),           % Initialize weights if nescessary
  W1 = rand(hidden,inputs+1)-0.5;
  W2 = rand(outputs,hidden+1)-0.5;
end
if (size(W1,2)~=inputs+1 | size(W1,1)~=hidden |... % Check dimensions
      size(W2,2)~=hidden+1 | size(W2,1)~=outputs)
   error('Dimension mismatch in weights, data, or NetDef');
end
y1       = zeros(hidden,1);             % Hidden layer outputs
y2       = zeros(outputs,1);            % Network output
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A useful vector!
PHI_aug  = [PHI;ones(1,N)];             % Augment PHI with a row containg ones
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters = parameters1 + parameters2; % Total # of weights
PSI        = zeros(parameters,outputs); % Deriv. of each output w.r.t. each weight
                                        % Parametervector containing all weights
theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
index3= 1:(reduced+1):(reduced^2);      % Yet another useful vector
I        = eye(outputs);                % (outputs|outputs) unity matrix

if nargin<6 | isempty(trparms) % Default training parameters
  trparms = settrain;
else                                    % User specified values
  trparmsdef = settrain;                
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
  if ~isfield(trparms,'method')
     trparms = settrain(trparms,'method','default');
  end  
end
if strcmp(trparms.method,'ff'),         % Forgetting factor method
  mflag     = 1;                        % Method flag
  if isfield(trparms,'lambda')          % Forgetting factor
     lambda = trparms.lambda;
  else lambda = trparmsdef.lambda; end
  if isfield(trparms,'p0'),             % Diagonal of covariance matrix
     p0 = trparms.p0;
  else p0 = trparmsdef.p0; end
  P         = p0 * eye(reduced);        % Initialize covariance matrix
  lambdaI  = lambda*I;                  % Useful diagonal matrix
  
elseif strcmp(trparms.method,'ct'),     % Constant trace method
  mflag     = 2;                        % Method flag
  if isfield(trparms,'alpha_max')       % Max eigenvalue
     alpha_max = trparms.alpha_max;
  else alpha_max = trparmsdef.alpha_max; end
  if isfield(trparms,'alpha_min'),      % Min eigenvalue
     alpha_min = trparms.alpha_min;
  else alpha_min = trparmsdef.alpha_min; end
  P      = alpha_max * eye(reduced);    % Initialize covariance matrix
  
elseif strcmp(trparms.method,'efra'),   % EFRA method
  mflag     = 3;                        % Method flag
  if isfield(trparms,'alpha')           % EFRA parameter 'alpha'
     alpha = trparms.alpha;
  else alpha = trparmsdef.alpha; end
  if isfield(trparms,'beta')            % EFRA parameter 'beta'
     beta = trparms.beta;
  else beta = trparmsdef.beta; end
  if isfield(trparms,'delta')           % EFRA parameter 'delta'
     delta = trparms.delta;
  else delta = trparmsdef.delta; end
  if isfield(trparms,'eflambda')        % EFRA parameter 'eflambda'
     lambda = trparms.eflambda;
  else lambda = trparmsdef.eflambda; end
  gamma     = (1-lambda)/lambda;
                                        % Max. eigenvalue
  maxeig = gamma/(2*delta)*(1+sqrt(1+4*beta*delta/(gamma*gamma)));
  P      = maxeig * eye(reduced);       % Initialize covariance matrix
  betaI  = beta*eye(reduced);           % Useful diagonal matrix
end
PI_vector= zeros(trparms.maxiter,1);    % Vector containing the normalized SSE for each iterat.
critdif  = trparms.critterm+1;          % Initialize stopping variables


%----------------------------------------------------------------------------------
%-------------                    TRAIN NETWORK                       -------------
%----------------------------------------------------------------------------------
clc;
c=fix(clock);
fprintf('Network training started at %2i.%2i.%2i\n\n',c(4),c(5),c(6));

for iteration=1:trparms.maxiter,
  SSE=0;
  for t=1:N,

% >>>>>>>>>>>>>>>>>>>>>>>>  COMPUTE NETWORK OUTPUT y2(theta) <<<<<<<<<<<<<<<<<<<<<<
    h1 = W1(:,1:inputs)*PHI(:,t) + W1(:,inputs+1);  
    y1(H_hidden) = pmntanh(h1(H_hidden)); 
    y1(L_hidden) = h1(L_hidden);
    
    h2 = W2(:,1:hidden)*y1 + W2(:,hidden+1);
    y2(H_output) = pmntanh(h2(H_output));
    y2(L_output) = h2(L_output);

    y1_aug=[y1;1];
    E = Y(:,t) - y2;                      % Training error


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  COMPUTE THE PSI MATRIX  <<<<<<<<<<<<<<<<<<<<<<<<<<<
% (The derivative of each y2(t) with respect to each weight)

    % ==========   Elements corresponding to the linear output units   ============
    for i = L_output'

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      index1 = (i-1) * (hidden + 1) + 1;
      PSI(index1:index1+hidden,i) = y1_aug;
      % ---------------------------------------------------------------------
 
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        PSI(index(j):index(j)+inputs,i) = W2(i,j)*PHI_aug(:,t);
      end
      
      for j = H_hidden',
        PSI(index(j):index(j)+inputs,i) = W2(i,j)*(1-y1(j)*y1(j))*PHI_aug(:,t);
      end 
      % ---------------------------------------------------------------------    
    end

    % ============  Elements corresponding to the tanh output units   =============
    for i = H_output',
      % -- The part of PSI corresponding to hidden-to-output layer weights --
      index1 = (i-1) * (hidden + 1) + 1;
      PSI(index1:index1+hidden,i) = y1_aug * (1 - y2(i)*y2(i));
      % ---------------------------------------------------------------------
       
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        PSI(index(j):index(j)+inputs,i) = W2(i,j)*(1-y2(i)*y2(i))...
                                                              * PHI_aug(:,t);
      end
      
      for j = H_hidden',
        PSI(index(j):index(j)+inputs,i) = W2(i,j)*(1-y1(j)*y1(j))...
                                             *(1-y2(i)*y2(i)) * PHI_aug(:,t);
      end
      % ---------------------------------------------------------------------
    end


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    UPDATE THE WEIGHTS    <<<<<<<<<<<<<<<<<<<<<<<<<<<
    PSI_red = PSI(theta_index,:);
    
    % ---------- Forgetting factor method ----------
    if mflag == 1,
      % -- Update P matrix --
      P = (P - P*PSI_red*inv(lambdaI + PSI_red'*P*PSI_red)*PSI_red'*P ) / lambda;

      % -- Update Parameters --
      theta_red = theta_red + P*PSI_red*E;
      
    % ----------  Constant trace method   ---------- 
    elseif mflag == 2,
      % -- Measurement update of P matrix --
      P = (P - P*PSI_red * inv(I + PSI_red'*P*PSI_red) * PSI_red'*P );

      % -- Update Parameters --
      theta_red = theta_red + P*PSI_red*E;

      % -- Time update of P matrix --
      P         = ((alpha_max-alpha_min)/trace(P))*P;
      P(index3) = P(index3)+alpha_min;
      
    % ----------       EFRA method        ---------- 
    else 
      % -- Correction factor --
      K = P*PSI_red * (alpha*inv(I + PSI_red'*P*PSI_red));

      % -- Update Parameters --
      theta_red = theta_red + K*E;
      
      % -- Update P --
      P = P/lambda - K*PSI_red'*P + betaI-delta*P*P;
    end
    theta(theta_index) = theta_red;       % Put estimated weights back into theta


    % -- Put the parameters back into the weight matrices --
    W1 = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
    W2 = reshape(theta(1:parameters2),hidden+1,outputs)';

    % -- Accumulate SSE --
    SSE = SSE + E'*E;
  end
  
  
%>>>>>>>>>>>>>>>>>>>>>>       UPDATES FOR NEXT ITERATION       <<<<<<<<<<<<<<<<<<<<
  PI = SSE/(2*N);
  PI_vector(iteration) = PI;            % Collect PI
  if iteration>1, 
     critdif  = abs(PI_vector(iteration-1)-PI);    % Criterion difference
  end   
  switch(trparms.infolevel)                                % Print on-line inform
       case 1
          fprintf('iteration # %i   W=%4.3e  critdif=%3.2e\n',iteration,PI,critdif);
       otherwise
          fprintf('iteration # %i   W = %4.3e\r',iteration,PI);
  end
  if (PI < trparms.critmin | critdif<trparms.critterm) % Check if stop condition is satisfied
     break
  end
end


%----------------------------------------------------------------------------------
%-------------              END OF NETWORK TRAINING                  --------------
%----------------------------------------------------------------------------------
c=fix(clock);
PI_vector = PI_vector(1:iteration);
fprintf('\n\nNetwork training ended at %2i.%2i.%2i\n',c(4),c(5),c(6));
