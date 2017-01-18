function [W1,W2,PI_vector,iteration,lambda]=...
                                       nnarxm(NetDef,NN,W1,W2,trparms,GAM,Y,U)
%  NNARXM
%  ------
%          Determine a nonlinear ARX model of a dynamic system by training a
%          two-layer neural network with the Marquardt method. The function
%          can handle multi-input multi-output systems (MIMO).
%
%  CALL:
%   [W1,W2,critvec,iteration,lambda]=nnarxm(NetDef,NN,W1,W2,trparms,Gamma,Y,U)
%
%  INPUTS: 
%  U       : Input to system  (left out in the nnar case)
%            A matrix with the structure: [(# of inputs) | (# of data)]
%  Y       : Output of the system. ((# of outputs)  | # of data)
%  NN      : NN=[na1 nb1 nk1;na2 nb2 nk2;...].
%            naX = # of past output X used for determining the prediction
%            nbX = # of past input X used to determine prediction
%            nkX = time delay (usually 1)
%            For multi-input systems, nb and nk contains as many columns as
%            there are inputs. Read the manual for more details. 
%  W1,W2   : Input-to-hidden layer and hidden-to-output layer weights.
%            If they are passed as [], they are automatically initialized
%  trparms: Data structure with parameters associated with the
%           training algorithm (optional). Use the function SETTRAIN if
%           you do not want to use the default values.
%  Gamma   : Covariance matrix.
%            If Gamma=[], Gamma is set to the identity matrix
%  
%           For time series (NNAR models), use NN=na.
%
%  See the function "marq" for an explanation of the remaining input arguments
%  as well as of the returned variables.
                                                                                 
%  Programmed by : Magnus Norgaard, IAU/IMM Technical University of Denmark
%  LastEditDate  : Jan. 8, 2000

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   INITIALIZATIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[ny,N]  = size(Y);
[ny,NNn]= size(NN);
na = NN(:,1);
if NNn==1
  nb = 0;                          % nnar model
  nk = 0;
  nu = 0;
  nab=na;
else
  [nu,N] = size(U); 
  nb     = NN(:,2:1+nu);           % nnarx model
  nk     = NN(:,2+nu:1+2*nu);
  if nu>1,
    nab  = na + sum(nb')';
  else
    nab    = na+nb;
  end
end

nmax        = max(max([na nb+nk-1]));


% -- Initialize weights if nescessary --
if isempty(W1) | isempty(W2),
  hidden = length(NetDef(1,:));    % Number of hidden neurons
  W1 = rand(hidden,sum(nab)+1)-0.5;
  W2 = rand(ny,hidden+1)-0.5;
end


% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI  <<<<<<<<<<<<<<<<<<<
PHI = zeros(sum(nab),N-nmax);
jj  = nmax+1:N;
index = 0;
for o=1:ny,
  for k = 1:na(o), PHI(k+index,:)    = Y(o,jj-k); end
  index = index+na(o);
  for kk = 1:nu,
    for k = 1:nb(o,kk), PHI(k+index,:) = U(kk,jj-k-nk(o,kk)+1); end
    index = index + nb(o,kk);
  end
end


% >>>>>>>>>>>>>>>>>>>>        CALL TRAINING FUNCTION        <<<<<<<<<<<<<<<<<<<<
if isempty(GAM),
  Ys=Y;
else
  S = sqrtm(inv(GAM));
  Ys = S*Y;
end
[W1,W2,PI_vector,iteration,lambda]=marq(NetDef,W1,W2,PHI,Ys(:,nmax+1:N),trparms);

