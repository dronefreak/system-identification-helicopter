function [W1,W2,Chat,PI_vector,iteration,lambda]=nnarmax1(NetDef,NN,W1,W2,Chat,...
                                                                     trparms,Y,U)
%  NNARMAX1
%  --------
%          [W1,W2,Chat,critvec,iteration,lambda]=...
%                              nnarmax1(NetDef,NN,W1,W2,Chat,trparms,Y,U)
%          Determines a nonlinear ARMAX model of a dynamic system by training
%          a two layer neural network with the Marquardt method. The function
%          can handle multi-input systems (MISO). It is assumed that the noise 
%          is filtered by a linear MA-filter. y(t)=f(Y(t-1),U(t-1)) + Ce(t)
%
%  INPUT:
%  U       : Input signal (= control signal) (left out in the nnarma case)
%            dim(U) = [(inputs) * (# of data)]
%  Y       : Output signal. dim(Y) = [1 * # of data]
%  NN      : NN=[na nb nc nk].
%            na = # of past outputs used for determining the prediction
%            nb = # of past inputs
%            nc = # of past residuals (= order of C)
%            nk = time delay (usually 1)
%            For multi-input systems, nb and nk contain as many columns as
%            there are inputs.
%  W1,W2   : Input-to-hidden layer and hidden-to-output layer weights.
%            If they are passed as [], they are initialized automatically
%  Chat    : Initial MA-filter estimate (initialized automatically if Chat=[])
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
Ndat     = length(Y);                   % # of data
na = NN(1);
if length(NN)==2                        % nnarma model
  nc     = NN(2);
  nb     = 0;
  nk     = 0;
  nu     = 0;
else                                    % nnarmax model
  [nu,Ndat]= size(U); 
  nb     = NN(2:1+nu);
  nc     = NN(2+nu);
  nk     = NN(2+nu+1:2+2*nu);
end
nmax     = max([na,nb+nk-1,nc]);        % 'Oldest' signal used as input to the model
N        = Ndat - nmax;                 % Size of training set
nab      = na+sum(nb);                  % na+nb
nabc     = nab+nc;                      % na+nb+nc
hidden   = length(NetDef(1,:));         % Number of hidden neurons
inputs   = nab;                         % Number of inputs to the network
outputs  = 1;                           % Only one output 
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
y1       = [zeros(hidden,N);ones(1,N)]; % Hidden layer outputs
y2       = zeros(outputs,N);            % Network output
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A useful vector!
index2 = (0:N-1)*outputs;               % Yet another useful vector
iteration= 1;                           % Counter variable
dw       = 1;                           % Flag telling that the weights are new
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters12=parameters1 + parameters2; % Total # of weights
parameters = parameters12+nc;           % Total # of parameters
ones_h   = ones(hidden+1,1);            % A vector of ones
ones_i   = ones(inputs+1,1);            % Another vector of ones
if isempty(W1) | isempty(W2),           % Initialize weights if nescessary
  trparmsi = settrain;
  trparmsi = settrain(trparmsi,'maxiter',100);
  if nb==0,
    [W1,W2]=nnarx(NetDef,na,[],[],trparmsi,Y);
  else
    [W1,W2]=nnarx(NetDef,[na nb nk],[],[],trparmsi,Y,U);
  end                            
end
if isempty(Chat),                       % Initialize a stable Chat if nescessary
  while sum(abs(roots(Chat))<1)<nc
    Chat = [1 rand(1,nc)-0.5];
  end
end
                                        % Parameter vector containing all weights
theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1) ; Chat(2:nc+1)'];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
index3   = 1:(reduced+1):(reduced^2);   % A third useful vector
PSI_red  = zeros(reduced,N);            % Deriv. of output w.r.t. each weight
RHO      = zeros(parameters12,N);       % Partial -"-  -"-
RHO2     = zeros(nc,N);                 % Partial deriv. of output wrt. each C-par.
lambda_old = 0;
if ~exist('trparms') | isempty(trparms) % Default training parameters
  trparms = settrain;
  lambda  = trparms.lambda;
  skip=trparms.skip+1;
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
      D = trparms.D([ones(1,parameters2) 2*ones(1,parameters1) ones(1,nc)])';
      D = D(theta_index);
    elseif length(trparms.D)>2,           % Individual weight decay
      D = trparms.D(:);
    end
  end
end
D = D(:);
N2       = N-skip+1;
critdif  = trparms.critterm+1;            % Initialize stopping variables
gradmax  = trparms.gradterm+1;
paramdif = trparms.paramterm+1;
PI_vector = zeros(trparms.maxiter,1);     % Vector for storing criterion values


% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
PHI = zeros(nab,N);
jj  = nmax+1:Ndat;
for k = 1:na, PHI(k,:)    = Y(jj-k); end
index4 = na;
for kk = 1:nu,
  for k = 1:nb(kk), PHI(k+index4,:) = U(kk,jj-k-nk(kk)+1); end
  index4 = index4 + nb(kk);
end
PHI_aug = [PHI;ones(1,N)];              % Augment PHI with a row containg ones
Y       = Y(nmax+1:Ndat);               % Extract the 'target' part of Y


%----------------------------------------------------------------------------------
%--------------                   TRAIN NETWORK                       -------------
%----------------------------------------------------------------------------------
clc;
c=fix(clock);
fprintf('Network training started at %2i.%2i.%2i\n\n',c(4),c(5),c(6));


% >>>>>>>>>>>>>>>>>>>>>  COMPUTE NETWORK OUTPUT  y2(theta)   <<<<<<<<<<<<<<<<<<<<<<
h1 = W1*PHI_aug;  
y1(H_hidden,:) = pmntanh(h1(H_hidden,:));
y1(L_hidden,:) = h1(L_hidden,:);

h2 = W2*y1;
y2(H_output,:) = pmntanh(h2(H_output,:));
y2(L_output,:) = h2(L_output,:);

Ebar     = Y - y2;                      % Error between Y and deterministic part
E        = filter(1,Chat,Ebar);         % Prediction error
SSE      = E(skip:N)*E(skip:N)';        % Sum of squared errors (SSE)
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
    RHO_red = RHO(theta_index(1:reduced-nc),:);


    % Partial deriv. of output wrt. each C-par.  
    for i=1:nc,
      RHO2(i,nmax+i-1:N) = E(1:N-i-nmax+2);
    end
  
  

% >>>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<<
    for i=1:reduced-nc,
      PSI_red(i,:) = filter(1,Chat,RHO_red(i,:));
    end
    for i=1:nc,
      PSI_red(reduced-nc+i,:) = filter(1,Chat,RHO2(i,:));
    end
  
   
% >>>>>>>>>>>>>>>>>>>>>>>>>>>        COMPUTE h_k        <<<<<<<<<<<<<<<<<<<<<<<<<<<
    % -- Gradient --
    G = PSI_red(:,skip:N)*E(skip:N)'-D.*theta_red;
    
    % -- Hessian  --
    H = PSI_red(:,skip:N)*PSI_red(:,skip:N)';
    H(index3) = H(index3)'+(lambda-lambda_old);     % Add diagonal matrix
    dw = 0;
  end
  
  % -- Hessian  --
  H(index3) = H(index3)'+(lambda-lambda_old);       % Add diagonal matrix

  % -- Search direction --
  h = H\G;                                          % Solve for search direction

  % -- Compute 'apriori' iterate --
  theta_red_new = theta_red + h;                    % Update parameter vector
  theta(theta_index) = theta_red_new;

  % -- Put the parameters back into the weight matrices and Chat --
  W1_new = reshape(theta(parameters2+1:parameters12),inputs+1,hidden)';
  W2_new = reshape(theta(1:parameters2),hidden+1,outputs)';
  Chat   = [1 theta(parameters12+1:parameters)'];
  croots = roots(Chat);
  for i=1:length(croots),
    if abs(croots(i))>1, croots(i)=1/croots(i); end
  end
  Chat=real(poly(croots));
  theta(parameters-nc+1:parameters)=Chat(2:nc+1)';
  theta_red_new = theta(theta_index);
  
    
% >>>>>>>>>>>>>>>>>>>>   COMPUTE NETWORK OUTPUT  y2(theta+h)   <<<<<<<<<<<<<<<<<<<<
  h1 = W1_new*PHI_aug;  
  y1(H_hidden,:) = pmntanh(h1(H_hidden,:));
  y1(L_hidden,:) = h1(L_hidden,:);
    
  h2 = W2_new*y1;
  y2(H_output,:) = pmntanh(h2(H_output,:));
  y2(L_output,:) = h2(L_output,:);
  Ebar           = Y - y2;                % Error between Y and deterministic part
  E_new    = filter(1,Chat,Ebar);         % Prediction error
  SSE_new  = E_new(skip:N)*E_new(skip:N)';% Sum of squared errors (SSE)
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
