function [W1,W2,PI_vector,iteration]=nnrarmx2(NetDef,NN,W1,W2,trparms,Y,U)
%  NNRARMX2 
%  ---------
%          Determines a nonlinear ARMAX model
%                          y(t)=f(y(t-1),...,u(t-k),...,e(t-1),...)
%          of a dynamic system by training a two layer neural network with
%          a recursive Gauss-Newton method. The function can handle multi-
%          input systems (MISO).
%
%  CALL:
%   [W1,W2,NSSEvec,iteration]=nnrarmx2(NetDef,NN,W1,W2,trparms,Y,U)
%
%  INPUTS:
%  U       : Input signal (= control signal) (left out in the nnarma case)
%            dim(U) = [(inputs) * (# of data)]
%  Y       : Output signal. dim(Y) = [1 * # of data]
%  NN      : NN=[na nb nc nk].
%            na = # of past outputs used for determining the prediction
%            nb = # of past inputs
%            nc = # of past residuals (= order of C)
%            nk = time delay (usually 1)
%            For multi-input systems nb and nk contain as many columns as
%            there are inputs.
%  W1,W2   : Input-to-hidden layer and hidden-to-output layer weights.
%            If they are passed as [] they are initialized automatically
%
%           For time series (NNARMA models), NN=[na nc] only.
% 
%  See the function RPE for an explanation of the remaining inputs as 
%  well as of the returned variables.
%                                                                                
%  Programmed by : Magnus Norgaard, IAU/IMM, Technical University of Denmark
%  LastEditDate  : Jan. 4, 2000

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   INITIALIZATIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<   
Ndat = length(Y);
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
inputs   = nabc;                        % Number of inputs to the network
outputs  = 1;                           % Only one output 
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
y1       = zeros(hidden,1);             % Hidden layer outputs
y2       = zeros(outputs,1);            % Network output
Eold     = zeros(nc,outputs);           % The nc past residuals
index = outputs*(hidden+1) + 1 + [0:hidden-1]*(inputs+1); % A useful vector!
parameters1= hidden*(inputs+1);         % # of input-to-hidden weights
parameters2= outputs*(hidden+1);        % # of hidden-to-output weights
parameters = parameters1 + parameters2; % Total # of weights
PSI        = zeros(parameters,outputs); % Deriv. of each output w.r.t. each weight
RHO        = zeros(parameters,outputs); % Partial -"-  -"-
RHO2       = zeros(nc,outputs);         % Partial deriv. of output wrt. each C-par.
if isempty(W1) | isempty(W2),           % Initialize weights if nescessary
  trparmsi = settrain;
  trparmsi = settrain(trparmsi,'maxiter',100); 
  if nb==0,
    [W1,W2]=nnarx(NetDef,na,[],[],trparmsi,Y);
  else
    [W1,W2]=nnarx(NetDef,[na nb nk],[],[],trparmsi,Y,U);
  end
  W1=[W1(:,1:nab) , 0.05*rand(hidden,nc)-0.025 W1(:,nab+1)];                           
end
                                  % Parametervector containing all weights
theta = [reshape(W2',parameters2,1) ; reshape(W1',parameters1,1)];
theta_index = find(theta);              % Index to weights<>0
theta_red = theta(theta_index);         % Reduced parameter vector
reduced  = length(theta_index);         % The # of parameters in theta_red
index3 = 1:(reduced+1):(reduced^2);     % Yet another useful vector
dy2de    = zeros(nc,1);                 % Der. of outputs wrt. the past residuals
index4   = nab+1:nabc;                  % And a fourth
I        = eye(outputs);                % (outputs|outputs) unity matrix
PSIold = zeros(reduced,nc);             % Past PSI vectors
if ~exist('trparms') | isempty(trparms) % Default training parameters
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
critdif  = trparms.critterm+1;          % Initialize stopping variables
PI_vector = zeros(trparms.maxiter,1);   % Vector for storing criterion values


% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
PHI = zeros(nabc,N);
jj  = nmax+1:Ndat;
for k = 1:na, PHI(k,:)    = Y(jj-k); end
index5 = na;
for kk = 1:nu,
  for k = 1:nb(kk), PHI(k+index5,:) = U(kk,jj-k-nk(kk)+1); end
  index5 = index5 + nb(kk);
end
PHI_aug  = [PHI;ones(1,N)];             % Augment PHI with a row containg ones
Y        = Y(nmax+1:Ndat);              % Extract the 'target' part of Y



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

    y1_aug = [y1;1];
    E   = Y(:,t) - y2;
    PHI_aug(nab+1:nabc,t) = Eold;


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  COMPUTE THE RHO MATRIX  <<<<<<<<<<<<<<<<<<<<<<<<<<<
% Partial derivative of output (y2) with respect to each weight and neglecting
% that the model inputs (the residuals) depends on the weights

    % ==========   Elements corresponding to the linear output units   ============
    for i = L_output'
      % -- The part of RHO corresponding to hidden-to-output layer weights --
      index1 = (i-1) * (hidden + 1) + 1;
      RHO(index1:index1+hidden,i) = y1_aug;
      % ---------------------------------------------------------------------
 
      % -- The part of RHO corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        RHO(index(j):index(j)+inputs,i) = W2(i,j)*PHI_aug(:,t);
      end
      
      for j = H_hidden',
        RHO(index(j):index(j)+inputs,i) = W2(i,j)*(1-y1(j)*y1(j))*PHI_aug(:,t);
      end 
      % ---------------------------------------------------------------------    
    end

    % ============  Elements corresponding to the tanh output units   =============
    for i = H_output',
      % -- The part of RHO corresponding to hidden-to-output layer weights --
      index1 = (i-1) * (hidden + 1) + 1;
      RHO(index1:index1+hidden,i) = y1_aug * (1 - y2(i)*y2(i));
      % ---------------------------------------------------------------------
       
      % -- The part of RHO corresponding to input-to-hidden layer weights ---
      for j = L_hidden',
        RHO(index(j):index(j)+inputs,i) = W2(i,j)*(1-y2(i)*y2(i))...
                                                              * PHI_aug(:,t);
      end
      
      for j = H_hidden',
        RHO(index(j):index(j)+inputs,i) = W2(i,j)*(1-y1(j)*y1(j))...
                                             *(1-y2(i)*y2(i)) * PHI_aug(:,t);
      end
      % ---------------------------------------------------------------------
    end
    RHO_red = RHO(theta_index(1:reduced));



% >>>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE THE PSI MATRIX   <<<<<<<<<<<<<<<<<<<<<<<<<<
    % ---------- Find derivative of output wrt. the past residuals ----------
    dy2dy1 = W2(:,1:hidden);
    for j = H_output',
        dy2dy1(j,:) = W2(j,1:hidden)*(1-y2(j,t).*y2(j,t));
    end
      
    % Matrix with partial derivatives of the output from each hidden neurons with
    % respect to each input:
    dy1de = W1(:,index4);
    for j = H_hidden',
      dy1de(j,:) = W1(j,index4)*(1-y1(j)*y1(j));
    end

    % Matrix with partial derivative of each output with respect to each input
    Chat = [1 (dy2dy1 * dy1de)];

    PSI_red=RHO_red;
    for i=1:nc,
      PSI_red = PSI_red-Chat(i+1)*PSIold(:,i);
    end

  
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    UPDATE THE WEIGHTS    <<<<<<<<<<<<<<<<<<<<<<<<<<<
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
    W1   = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
    W2   = reshape(theta(1:parameters2),hidden+1,outputs)';
 
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
  Eold = [E;Eold(1:nc-1)];              % Past residuals
  PSIold = [PSI_red,PSIold(:,1:nc-1)];  % Past gradients
end


%----------------------------------------------------------------------------------
%-------------              END OF NETWORK TRAINING                  --------------
%----------------------------------------------------------------------------------
c=fix(clock);
PI_vector = PI_vector(1:iteration);
fprintf('\n\nNetwork training ended at %2i.%2i.%2i\n',c(4),c(5),c(6));
