function [W1,W2,PI_vector,iteration,lambda]=marqlm(NetDef,W1,W2,PHI,Y,trparms)
%  MARQLM
%  ------
%            Levenberg-Marquardt training algorithm that uses less memory
%            than MARQ but is slower. The difference in speed occurs because
%            the function is less 'vectorized' (which is a MATLAB problem)
%            but also because some calculations are made more than once.
%
%  Written by : Magnus Norgaard, IAU/IMM, Technical Univ. of Denmark
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
   error('Dimension mismatch in weights, data, or NetDef.');
end
y1       = [zeros(hidden,1);1];         % Hidden layer outputs
y2       = zeros(outputs,1);            % Network output
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A useful vector!
iteration = 1;                          % Counter variable
dw       = 1;                           % Flag telling that the weights are new
PHI      = [PHI;ones(1,N)];             % Augment PHI with a row containg ones
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters = parameters1 + parameters2; % Total # of weights
PSI      = zeros(parameters,outputs);   % Deriv. of each output w.r.t. each weight                                     % Parameter vector containing all weights
theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
index3   = 1:(reduced+1):(reduced^2);   % A third useful vector
lambda_old = 0;
if nargin<6 | isempty(trparms) % Default training parameters
  trparms = settrain;
  lambda  = trparms.lambda;
  D       = trparms.D;
else                                    % User specified values
  if ~isstruct(trparms),
     error('''trparms'' must be a structure variable.');
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
      D = trparms.D(ones(1,reduced));      
    elseif length(trparms.D)==2,          % Two weight decay parameters
      D = trparms.D([ones(1,parameters2) 2*ones(1,parameters1)])';
      D = D(theta_index);
    elseif length(trparms.D)>2,           % Individual weight decay
      D = trparms.D(:);
    end
  end
end
D = D(:);
critdif  = trparms.critterm+1;            % Initialize stopping variables
gradmax  = trparms.gradterm+1;
paramdif = trparms.paramterm+1;
PI_vector = zeros(trparms.maxiter,1);     % Vector for storing criterion values

%----------------------------------------------------------------------------------
%--------------                   TRAIN NETWORK                       -------------
%----------------------------------------------------------------------------------
clc;
c=fix(clock);
fprintf('Network training started at %2i.%2i.%2i\n\n',c(4),c(5),c(6));


% >>>>>>>>>>>>>>>>>>>>>  COMPUTE NETWORK OUTPUT  y2(theta)   <<<<<<<<<<<<<<<<<<<<<<
SSE = 0;
for t=1:N,
  h1 = W1*PHI(:,t);  
  y1(H_hidden) = pmntanh(h1(H_hidden));
  y1(L_hidden) = h1(L_hidden);

  h2 = W2*y1;
  y2(H_output) = pmntanh(h2(H_output));
  y2(L_output) = h2(L_output);

  E        = Y(:,t) - y2;                   % Training error
  SSE      = SSE + E'*E;                  % Sum of squared errors (SSE)
end
PI       = (SSE+theta_red'*(D.*theta_red))/(2*N); % Performance index

while (iteration<=trparms.maxiter & PI>trparms.critmin & lambda<1e7 & ...
       (critdif>trparms.critterm | gradmax>trparms.gradterm | ...
       paramdif>trparms.paramterm))   
if dw==1,
% >>>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<<
% (The derivative of each network output (y2) with respect to each weight)
  G = zeros(reduced,1);
  H = diag(D);
  for t=1:N,
    h1 = W1*PHI(:,t);  
    y1(H_hidden) = pmntanh(h1(H_hidden));
    y1(L_hidden) = h1(L_hidden);    
    h2 = W2*y1;
    y2(H_output) = pmntanh(h2(H_output));
    y2(L_output) = h2(L_output);
    E        = Y(:,t) - y2;                      % Training error

    % ==========   Elements corresponding to the linear output units   ============
    for i = L_output'
      index1 = (i-1) * (hidden + 1) + 1;

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      PSI(index1:index1+hidden,i) = y1;
      % ---------------------------------------------------------------------
 
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        PSI(index(j):index(j)+inputs,i) = W2(i,j)*PHI(:,t);
      end

      for j = H_hidden',
        PSI(index(j):index(j)+inputs,i) = W2(i,j)*(1-y1(j).*y1(j))*PHI(:,t);
      end 
      % ---------------------------------------------------------------------    
    end

    % ============  Elements corresponding to the tanh output units   =============
    for i = H_output',
      index1 = (i-1) * (hidden + 1) + 1;

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      PSI(index1:index1+hidden,i) = y1*(1 - y2(i,:).*y2(i,:));
      % ---------------------------------------------------------------------
     
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        PSI(index(j):index(j)+inputs,i) = W2(i,j)*(1-y2(i).*y2(i))* PHI(:,t);
      end
      
      for j = H_hidden',
        PSI(index(j):index(j)+inputs,i) = W2(i,j)*(1-y1(j).*y1(j))...
                                           *(1-y2(i).*y2(i))* PHI(:,t);
      end
      % ---------------------------------------------------------------------
    end
    PSI_red = PSI(theta_index,:);
    G = G + PSI_red*E;             % -- Gradient --
    H = H + PSI_red*PSI_red';      % -- Hessian  --
  end
  G = G - D.*theta_red;
  H(index3) = H(index3)'+ D;  
  dw = 0;
  end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>        COMPUTE h_k        <<<<<<<<<<<<<<<<<<<<<<<<<<<
  % -- Hessian  --
  H(index3) = H(index3)'+(lambda-lambda_old);  % Add diagonal matrix

  % -- Search direction --
  U = chol(H);
  x = U'\G;
  h = U\x;
%  h = H\G;                                     % Solve for search direction

  % -- Compute 'apriori' iterate --
  theta_red_new = theta_red + h;               % Update parameter vector
  theta(theta_index) = theta_red_new;

  % -- Put the parameters back into the weight matrices --
  W1_new = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
  W2_new = reshape(theta(1:parameters2),hidden+1,outputs)';
    
    
% >>>>>>>>>>>>>>>>>>>>   COMPUTE NETWORK OUTPUT  y2(theta+h)   <<<<<<<<<<<<<<<<<<<<
  SSE_new = 0;
  for t=1:N,
    h1 = W1_new*PHI(:,t);  
    y1(H_hidden) = pmntanh(h1(H_hidden));
    y1(L_hidden) = h1(L_hidden);
    
    h2 = W2_new*y1;
    y2(H_output) = pmntanh(h2(H_output));
    y2(L_output) = h2(L_output);

    E_new    = Y(:,t) - y2;                     % Training error
    SSE_new  = SSE_new + E_new'*E_new;     % Sum of squared errors (SSE)
  end
  PI_new   = (SSE_new + theta_red_new'*(D.*theta_red_new))/(2*N); % PI


% >>>>>>>>>>>>>>>>>>>>>>>>>>>       UPDATE  lambda     <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  L = h'*G + h'*(h.*(D+lambda));
  lambda_old = lambda;

  % Decrease lambda if SSE has fallen 'sufficiently'
  if 2*N*(PI - PI_new) > (0.75*L),
    lambda = lambda/2;
  
  % Increase lambda if SSE has grown 'sufficiently'
  elseif 2*N*(PI-PI_new) <= (0.25*L),
    lambda = 2*lambda;
  end


% >>>>>>>>>>>>>>>>>>>>       UPDATES FOR NEXT ITERATION        <<<<<<<<<<<<<<<<<<<<
  % Update only if criterion has decreased
  if PI_new < PI,                      
    critdif  = PI-PI_new;                           % Criterion difference
    gradmax  = max(abs(G))/N;                       % Maximum gradient
    paramdif = max(abs(theta_red_new - theta_red)); % Maximum parameter difference
    W1 = W1_new;
    W2 = W2_new;
    theta_red = theta_red_new;
    PI = PI_new;
    lambda_old=0;
    dw = 1;
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
