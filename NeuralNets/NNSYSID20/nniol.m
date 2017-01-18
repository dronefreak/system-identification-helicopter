function [W1f,W2f,W1g,W2g,PI_vector,iteration,lambda]=nniol(NetDeff,NetDefg,...
                                                  NN,W1f,W2f,W1g,W2g,trparms,Y,U)
%  NNIOL
%  -----
%       Train a neural network to model a dynamic system, assuming the 
%       following structure:
%         y(t) = f(y(t-1),..,y(t-na),u(t-nk-1),...,u(t-nk-nb)) 
%                   +  g(y(t-1),..,y(t-na),u(t-nk-1),...,u(t-nk-nb))*u(t-nk)
%       with the Levenberg-Marquardt method. This type of model description
%       is particularly relevant in control by discrete input-output
%       linearization.
%
%  CALL:
%  [W1f,W2f,W1g,W2g,critvec,iteration,lambda]=...
%                     nniol(NetDeff,NetDefg,NN,W1f,W2f,W1g,W2g,trparms,Y,U)
%
%  INPUTS:
%  U       : Input signal (= control signal) (left out in the nnarma case)
%            dim(U) = [(inputs) * (# of data)]
%  Y       : Output signal. dim(Y) = [1 * # of data]
%  NN      : NN=[na nb nk].
%            na = # of past outputs used for determining the prediction
%            nb = # of past inputs used for determining prediction
%            nk = time delay (usually 1)
%  NetDeff : Architecture of network used for modelling the function f
%  NetDefg : Archtecture of network used for modelling the function g
%  W1f,W2f : Input-to-hidden layer and hidden-to-output layer weights for 
%  W1g,W2g   the "f" and "g" nets, respectively.
%            If they are passed as [], they will be initialized automatically
%  trparms:  Data structure with parameters associated with the
%            training algorithm (optional). Use the function SETTRAIN if
%            you do not want to use the default values.
%
%            For time series (NNARMA models), use only NN=[na nc].
% 
%  See the function MARQ for an explanation of the remaining input arguments
%  as well as of the returned variables.

%  Programmed by : Magnus Norgaard, IAU/IMM, technical University of Denmark
%  LastEditDate  : January 2, 2000

%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------

% >>>>>>>>>>>>>>>>>>>>>>>>>>>     REGRESSOR STRUCTURE     <<<<<<<<<<<<<<<<<<<<<<<<<   
[nu,N]      = size(U);                   % # of inputs and # of data
na = NN(1);                              % # of past y's to be used in TDL
nb = NN(2:1+nu);                         % # of past u's to be used in TDL
nk = NN(2+nu:1+2*nu);                    % Time delay in system
nmax        = max([na nb+nk-1]);
nab         = na+sum(nb);                % Number of "signal inputs" to each net
outputs  = 1;                            % # of outputs
inputs   = nab-1;                        % # of inputs


% >>>>>>>>>>>>>>>>>>>>>>>>>>    OPTIONAL INITIALIZATIONS    <<<<<<<<<<<<<<<<<<<<<<<
% -- Initialize weights if nescessary --
if isempty(W1f) | isempty(W2f) | isempty(W1g) | isempty(W2g)
  hiddenf = length(NetDeff(1,:));        % Number of hidden neurons in f-net
  W1f = rand(hiddenf,nab)-0.5;
  W2f = rand(1,hiddenf+1)-0.5;
  hiddeng = length(NetDefg(1,:));        % Number of hidden neurons in g-net
  W1g = rand(hiddeng,nab)-0.5;
  W2g = rand(1,hiddeng+1)-0.5;
end


% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
PHI = zeros(nab-1,N-nmax);
jj  = nmax+1:N;
for k = 1:na, PHI(k,:)    = Y(jj-k); end
index = na;
for kk = 1:nu,                          % Here nu=1, allways!
for k = 1:nb(kk)-1, PHI(k+index,:) = U(kk,jj-k-nk(kk)); end
  index = index + nb(kk);
end
Y = Y(nmax+1:N);
U = U(nmax+1-nk:N-nk);
N = N-nmax;


% >>>>>>>>>>>>>>>>>>>>>     DETERMINE NETWORK STRUCTURE      <<<<<<<<<<<<<<<<<<<<<<
% ---------- f-net structure ----------
L_hiddenf = find(NetDeff(1,:)=='L')';   % Location of linear hidden neurons
H_hiddenf = find(NetDeff(1,:)=='H')';   % Location of tanh hidden neurons
L_outputf = find(NetDeff(2,:)=='L')';   % Location of linear output neurons
H_outputf = find(NetDeff(2,:)=='H')';   % Location of tanh output neurons
y1f       =[zeros(hiddenf,N);ones(1,N)];% Hidden layer outputs
y2f       = zeros(outputs,N);           % Network output          
indexf = outputs*(hiddenf+1) + 1 + [0:hiddenf-1]*(inputs+1); % A useful vector!
params1f  = hiddenf*(inputs+1);         % # of input-to-hidden weights
params2f  = outputs*(hiddenf+1);        % # of hidden-to-output weights
paramsf   = params1f + params2f;        % Total # of weights in net f
PSIf      = zeros(paramsf,outputs*N);   % Deriv. of fhat w.r.t. each weight
onesf_h   = ones(hiddenf+1,1);          % A vector of ones
onesf_i   = ones(inputs+1,1);           % Another vector of ones

% ---------- g-net structure ----------
L_hiddeng = find(NetDefg(1,:)=='L')';   % Location of linear hidden neurons
H_hiddeng = find(NetDefg(1,:)=='H')';   % Location of tanh hidden neurons
L_outputg = find(NetDefg(2,:)=='L')';   % Location of linear output neurons
H_outputg = find(NetDefg(2,:)=='H')';   % Location of tanh output neurons
y1g       =[zeros(hiddeng,N);ones(1,N)];% Hidden layer outputs
y2g       = zeros(outputs,N);           % Network output
indexg = outputs*(hiddeng+1) + 1 + [0:hiddeng-1]*(inputs+1); % A useful vector!
params1g  = hiddeng*(inputs+1);         % # of input-to-hidden weights
params2g  = outputs*(hiddeng+1);        % # of hidden-to-output weights
paramsg   = params1g + params2g;        % Total # of weights in net g
PSIg      = zeros(paramsg,outputs*N);   % Deriv. of ghat w.r.t. each weight
onesg_h   = ones(hiddeng+1,1);          % A vector of ones
onesg_i   = ones(inputs+1,1);           % Another vector of ones
ones_u    = ones(paramsg,1);            % Yet another vector of ones


% >>>>>>>>>>>>>>>>>>>>     MISCELLANOUS INITIALIZATIONS      <<<<<<<<<<<<<<<<<<<<<
index2 = (0:N-1)*outputs;               % Yet another useful vector
iteration = 1;                          % Counter variable
dw        = 1;                          % Flag telling that the weights are new
PHI_aug   = [PHI;ones(1,N)];            % Augment PHI with a row containg ones
parameters = paramsf + paramsg;         % Total # of weights
PSI      = zeros(parameters,outputs*N); % Deriv. of each output w.r.t. each weight
identity = eye(parameters);             % Identity matrix
                                        % Parameter vector containing all weights
theta = [reshape(W2f',params2f,1) ; reshape(W1f',params1f,1)];
theta = [theta ; reshape(W2g',params2g,1) ; reshape(W1g',params1g,1)];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
index3   = 1:(reduced+1):(reduced^2);   % A third useful vector


% >>>>>>>>>>>>>>>>>>>>>>>>>    TRAINING PARAMETERS     <<<<<<<<<<<<<<<<<<<<<<<<<<<
% -- Initialize 'trparms' if nescessary --
lambda_old = 0;
if ~exist('trparms') | isempty(trparms) % Default training parameters
  trparms = settrain;
  lambda  = trparms.lambda;
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
h1f = W1f*PHI_aug;  
y1f(H_hiddenf,:) = pmntanh(h1f(H_hiddenf,:));
y1f(L_hiddenf,:) = h1f(L_hiddenf,:);    
h2f = W2f*y1f;
y2f(H_outputf,:) = pmntanh(h2f(H_outputf,:));
y2f(L_outputf,:) = h2f(L_outputf,:);

h1g = W1g*PHI_aug;  
y1g(H_hiddeng,:) = pmntanh(h1g(H_hiddeng,:));
y1g(L_hiddeng,:) = h1g(L_hiddeng,:);    
h2g = W2g*y1g;
y2g(H_outputg,:) = pmntanh(h2g(H_outputg,:));
y2g(L_outputg,:) = h2g(L_outputg,:);

y2       = y2f + y2g.*U;                % Network output
E        = Y - y2;                      % Training error
E_vector = E(:);                        % Reshape E into a long vector
SSE      = E_vector'*E_vector;          % Sum of squared errors (SSE)
PI       = (SSE+theta_red'*(D.*theta_red))/(2*N); % Performance index

% Iterate until stopping criterion is satisfied
while (iteration<=trparms.maxiter & PI>trparms.critmin & lambda<1e7 & ...
       (critdif>trparms.critterm | gradmax>trparms.gradterm | ...
       paramdif>trparms.paramterm))   
if dw==1,
% >>>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<<
% (The derivative of fhat with respect to each weight in the f net)

    % ==========   Elements corresponding to the linear output units   ============
    for i = L_outputf'
      index1 = (i-1) * (hiddenf + 1) + 1;

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      PSIf(index1:index1+hiddenf,index2+i) = y1f;
      % ---------------------------------------------------------------------
 
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hiddenf',
        PSIf(indexf(j):indexf(j)+inputs,index2+i) = W2f(i,j)*PHI_aug;
      end
     
      for j = H_hiddenf',
        tmp = W2f(i,j)*(1-y1f(j,:).*y1f(j,:)); 
        PSIf(indexf(j):indexf(j)+inputs,index2+i) = tmp(onesf_i,:).*PHI_aug;
      end 
      % ---------------------------------------------------------------------    
    end

    % ============  Elements corresponding to the tanh output units   =============
    for i = H_outputf',
      index1 = (i-1) * (hiddenf + 1) + 1;

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      tmp = 1 - y2f(i,:).*y2f(i,:);
      PSIf(index1:index1+hiddenf,index2+i) = y1f.*tmp(onesf_h,:);
      % ---------------------------------------------------------------------
         
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hiddenf',
        tmp = W2f(i,j)*(1-y2f(i,:).*y2f(i,:));
        PSIf(indexf(j):indexf(j)+inputs,index2+i) = tmp(onesf_i,:).* PHI_aug;
      end
      
      for j = H_hiddenf',
        tmp  = W2f(i,j)*(1-y1f(j,:).*y1f(j,:));
        tmp2 = (1-y2f(i,:).*y2f(i,:));
        PSIf(indexf(j):indexf(j)+inputs,index2+i) = tmp(onesf_i,:)...
                                                  .*tmp2(onesf_i,:).* PHI_aug;
      end
      % ---------------------------------------------------------------------
    end


    % (The derivative of ghat with respect to each weight in the g net)

    % ==========   Elements corresponding to the linear output units   ============
    for i = L_outputg'
      index1 = (i-1) * (hiddeng + 1) + 1;

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      PSIg(index1:index1+hiddeng,index2+i) = y1g;
      % ---------------------------------------------------------------------
 
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hiddeng',
        PSIg(indexg(j):indexg(j)+inputs,index2+i) = W2g(i,j)*PHI_aug;
      end
     
      for j = H_hiddeng',
        tmp = W2g(i,j)*(1-y1g(j,:).*y1g(j,:)); 
        PSIg(indexg(j):indexg(j)+inputs,index2+i) = tmp(onesg_i,:).*PHI_aug;
      end 
      % ---------------------------------------------------------------------    
    end

    % ============  Elements corresponding to the tanh output units   =============
    for i = H_outputg',
      index1 = (i-1) * (hiddeng + 1) + 1;

      % -- The part of PSI corresponding to hidden-to-output layer weights --
      tmp = 1 - y2g(i,:).*y2g(i,:);
      PSIg(index1:index1+hiddeng,index2+i) = y1g.*tmp(onesg_h,:);
      % ---------------------------------------------------------------------
         
      % -- The part of PSI corresponding to input-to-hidden layer weights ---
      for j = L_hiddeng',
        tmp = W2g(i,j)*(1-y2g(i,:).*y2g(i,:));
        PSIg(indexg(j):indexg(j)+inputs,index2+i) = tmp(onesg_i,:).* PHI_aug;
      end
      
      for j = H_hiddeng',
        tmp  = W2g(i,j)*(1-y1g(j,:).*y1g(j,:));
        tmp2 = (1-y2g(i,:).*y2g(i,:));
        PSIg(indexg(j):indexg(j)+inputs,index2+i) = tmp(onesg_i,:)...
                                                  .*tmp2(onesg_i,:).* PHI_aug;
      end
      % ---------------------------------------------------------------------
    end
    PSI=[PSIf;PSIg.*U(ones_u,:)];
    PSI_red = PSI(theta_index,:);
    
    % -- Gradient --
    G = PSI_red*E_vector-D.*theta_red;

    % -- Hessian  --
    H = PSI_red*PSI_red';
    H(index3) = H(index3)'+D;                  % Add diagonal matrix
    dw = 0;
  end
  
   
% >>>>>>>>>>>>>>>>>>>>>>>>>>>        COMPUTE h_k        <<<<<<<<<<<<<<<<<<<<<<<<<<<
  % -- L-M Hessian  --
  H(index3) = H(index3)'+(lambda-lambda_old);       % Add diagonal matrix

  % -- Search direction --
  h = H\G;                                          % Solve for search direction

  % -- Compute 'apriori' iterate --
  theta_red_new = theta_red + h;                    % Update parameter vector
  theta(theta_index) = theta_red_new;

  % -- Put the parameters back into the weight matrices --
  W1f_new = reshape(theta(params2f+1:paramsf),inputs+1,hiddenf)';
  W2f_new = reshape(theta(1:params2f),hiddenf+1,outputs)';
  W1g_new = reshape(theta(paramsf+params2g+1:parameters),inputs+1,hiddeng)';
  W2g_new = reshape(theta(paramsf+1:paramsf+params2g),hiddeng+1,outputs)';
    
% >>>>>>>>>>>>>>>>>>>>   COMPUTE NETWORK OUTPUT  y2(theta+h)   <<<<<<<<<<<<<<<<<<<<
  h1f = W1f_new*PHI_aug;  
  y1f(H_hiddenf,:) = pmntanh(h1f(H_hiddenf,:));
  y1f(L_hiddenf,:) = h1f(L_hiddenf,:);
  h2f = W2f_new*y1f;
  y2f(H_outputf,:) = pmntanh(h2f(H_outputf,:));
  y2f(L_outputf,:) = h2f(L_outputf,:);

  h1g = W1g_new*PHI_aug;  
  y1g(H_hiddeng,:) = pmntanh(h1g(H_hiddeng,:));
  y1g(L_hiddeng,:) = h1g(L_hiddeng,:);   
  h2g = W2g_new*y1g;
  y2g(H_outputg,:) = pmntanh(h2g(H_outputg,:));
  y2g(L_outputg,:) = h2g(L_outputg,:);

  y2           = y2f + y2g.*U;
  E_new        = Y - y2;                 % Training error
  E_new_vector = E_new(:);               % Reshape E into a long vector
  SSE_new  = E_new_vector'*E_new_vector; % Sum of squared errors (SSE)
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
    paramdif = max(abs(theta_red_new - theta_red)); % Maximum parameter dif.
    W1f = W1f_new;
    W2f = W2f_new;
    W1g = W1g_new;
    W2g = W2g_new;
    theta_red = theta_red_new;
    E_vector = E_new_vector;
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

