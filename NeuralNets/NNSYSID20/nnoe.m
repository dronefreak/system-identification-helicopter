function [W1,W2,PI_vector,iteration,lambda]=nnoe(NetDef,NN,W1,W2,trparms,Y,U)
%  NNOE  
%  ---- 
%          Determine a nonlinear output error model of a dynamic system
%          by training a two layer neural network with the Marquardt method.
%
%            yhat(t)=f(yhat(t-1),...,yhat(t-na),u(t-nk),...u(t-nk-nb+1))
%          The function can handle multi-input systems (MISO).
%
%  CALL:
%  [W1,W2,critvec,iteration,lambda]=nnoe(NetDef,NN,W1,W2,trparms,Y,U)
%
%  INPUTS:
%  U       : Input signal (= control signal).
%            dim(U) = [(inputs) * (# of data)]
%  Y       : Output signal. dim(Y) = [1 * # of data]
%  NN      : NN=[na nb nk]
%            na = # of past outputs used for determining prediction
%            nb = # of past inputs used for determining prediction
%            nk = time delay (usually 1)
%            For multi-input systems nb and nk contain as many columns as
%            there are inputs.
%  W1,W2   : Input-to-hidden layer and hidden-to-output layer weights.
%            If they are passed as [] they are initialized automatically 
%  trparms:  Data structure with parameters associated with the
%            training algorithm (optional). Use the function SETTRAIN if
%            you do not want to use the default values.
% 
%  See the function MARQ for an explanation of the remaining input arguments
%  as well as of the returned variables.

%  Programmed by : Magnus Norgaard, IAU/IMM, technical University of Denmark
%  LastEditDate  : June 1, 2001

%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------
Ndat     = length(Y);                   % # of data
na       = NN(1);                       % Order of polynomials
[nu,Ndat]= size(U); 
nb       = NN(2:1+nu);
nk       = NN(1+nu+1:1+2*nu);
nmax     = max([na nb+nk-1]);           % Oldest signal used as input to the model
N        = Ndat - nmax;                 % Size of training set
nab      = na+sum(nb);                  % na+nb
hidden   = length(NetDef(1,:));         % Number of hidden neurons
inputs   = nab;                         % Number of inputs to the network
outputs  = 1;                           % Only one output 
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
y1       = zeros(hidden,N);             % Hidden layer outputs
y1       = [y1;ones(1,N)];
y2       = zeros(outputs,N);            % Network output
E        = zeros(outputs,N);            % Initialize prediction error vector
E_new    = zeros(outputs,N);            % Initialize prediction error vector
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A useful vector!
index2  = (0:N-1)*outputs;              % Yet another useful vector
iteration= 1;                           % Counter variable
dw       = 1;                           % Flag telling that the weights are new
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters=parameters1 + parameters2;   % Total # of weights
ones_h   = ones(hidden+1,1);            % A vector of ones
ones_i   = ones(inputs+1,1);            % Another vector of ones
if isempty(W1) | isempty(W2),           % Initialize weights if nescessary
   trparmsi = settrain;
   trparmsi = settrain(trparmsi,'maxiter',100);
   [W1,W2]=nnarx(NetDef,[na nb nk],[],[],trparmsi,Y,U);                   
end
                                        % Parameter vector containing all weights
theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
index3   = 1:(reduced+1):(reduced^2);   % A third useful vector
dy2dy    = zeros(na,N);                 % Der. of output wrt. the past outputs
dy1dy    = zeros(hidden,na);            % Der. of hidden unit outp. wrt. past outputs
index4   = 1:na;                        % And a fourth
PSI_red  = zeros(reduced,N);            % Deriv. of output w.r.t. each weight
RHO      = zeros(parameters,N);         % Partial -"-  -"-

lambda_old = 0;
if ~exist('trparms') | isempty(trparms) % Default training parameters
  trparms = settrain;
  lambda  = trparms.lambda;
  D       = trparms.D;
  skip    = trparms.skip+1;
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
N2       = N-skip+1;
critdif  = trparms.critterm+1;            % Initialize stopping variables
gradmax  = trparms.gradterm+1;
paramdif = trparms.paramterm+1;
PI_vector = zeros(trparms.maxiter,1);     % Vector for storing criterion values

% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
PHI = zeros(nab,N);
jj  = nmax+1:Ndat;
for k = 1:na, PHI(k,:)    = Y(jj-k); end
index5 = na;
for kk = 1:nu,
  for k = 1:nb(kk), PHI(k+index5,:) = U(kk,jj-k-nk(kk)+1); end
  index5 = index5 + nb(kk);
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
for t=1:N,
  h1 = W1*PHI_aug(:,t);  
  y1(H_hidden,t) = pmntanh(h1(H_hidden));
  y1(L_hidden,t) = h1(L_hidden);    

  h2 = W2*y1(:,t);
  y2(H_output,t) = pmntanh(h2(H_output,:));
  y2(L_output,t) = h2(L_output,:);

  for d=1:min(na,N-t),
    PHI_aug(d,t+d) = y2(:,t);
  end
end
E = Y - y2;                                        % Prediction error
SSE      = E(skip:N)*E(skip:N)';                   % Sum of squared errors (SSE)
PI       = (SSE+theta_red'*(D.*theta_red))/(2*N2); % Performance index

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
    % ---------- Find derivative of output wrt. the past outputs ----------
    for t=1:N,
      dy2dy1 = W2(:,1:hidden);
      for j = H_output',
        dy2dy1(j,:) = W2(j,1:hidden)*(1-y2(j,t).*y2(j,t));
      end

      % Matrix of partial derivatives of the output from each hidden unit with
      % respect to each input:
      dy1dy(L_hidden,:) = W1(L_hidden,index4);
      for j = H_hidden',
        dy1dy(j,:) = W1(j,index4)*(1-y1(j,t).*y1(j,t));
      end

      % Matrix of partial derivatives of each output with respect to each input
      dy2dy(:,t)= (dy2dy1 * dy1dy)';
    end


    % ---------- Determine PSI by "filtering" ----------
    for t=1:N,
      PSI_red(:,t)=RHO_red(:,t);
      for t1=1:min(na,t-1),
        PSI_red(:,t)  = PSI_red(:,t)+dy2dy(t1,t)*PSI_red(:,t-t1);
      end
    end
 
   
% >>>>>>>>>>>>>>>>>>>>>>>>>>>        COMPUTE h_k        <<<<<<<<<<<<<<<<<<<<<<<<<<<
    % -- Gradient --
    G = PSI_red(:,skip:N)*E(skip:N)'-D.*theta_red;

    % -- Hessian  --
    H = PSI_red(:,skip:N)*PSI_red(:,skip:N)';
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
  h1 = W1_new*PHI_aug(:,t);  
  y1(H_hidden,t) = pmntanh(h1(H_hidden));
  y1(L_hidden,t) = h1(L_hidden);    

  h2 = W2_new*y1(:,t);
  y2(H_output,t) = pmntanh(h2(H_output,:));
  y2(L_output,t) = h2(L_output,:);
  
  for d=1:min(na,N-t),
    PHI_aug(d,t+d) = y2(:,t);
  end
end
E_new    = Y - y2;                      % Prediction error
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
