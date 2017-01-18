function [W1,W2,GAMMA,lambda]=nnigls(NetDef,NN,W1,W2,trparms,GAMMA,Y,U)
%  NNIGLS
%  ------
%          Train a multi-output NNARX model using an iterated generalized
%          least squares method (IGLS).
%
%  CALL:
%     [W1,W2,Gamma,lambda] = nnigls(NetDef,NN,W1,W2,trparms,Gamma0,Y,U)
%
%  INPUTS: 
%    U,Y,NN,W1,W2 : See NNARXM
%    trparms: Data structure with parameters associated with the
%             training algorithm (optional). Use the function SETTRAIN if
%             you do not want to use the default values.
%    Gamma0 : Covariance matrix. If passed as [] it is initialized to 
%             the identity matrix.
%
%  OUTPUTS:
%    W1, W2, lambda: See MARQ
%    Gamma         : Estimated covariance matrix
                                                                                  
%  Written by : Magnus Norgaard, IAU/IMM technical University of Denmark
%  LastEditDate  : Jan 8, 2000

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   INITIALIZATIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[ny,N]  = size(Y);                 % Size of data set
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
if isempty(trparms) % Default training parameters
   trd = settrain(trparms,'repeat','default');
   repeat=trd.repeat;
else
  if ~isfield(trparms,'repeat')
     trparms = settrain(trparms,'repeat','default');
  end
  repeat    = trparms.repeat;
end

nmax        = max(max([na nb+nk-1]));


% -- Initialize weights if nescessary --
if isempty(W1) | isempty(W2),
  hidden = length(NetDef(1,:));    % Number of hidden neurons
  W1 = rand(hidden,sum(nab)+1)-0.5;
  W2 = rand(ny,hidden+1)-0.5;
end


% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
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


% >>>>>>>>>>>>>>>>>>>>         CALL TRAINING FUNCTION         <<<<<<<<<<<<<<<<<<<<<
if isempty(GAMMA), GAMMA=eye(ny); end
GAMMAi = inv(GAMMA);
for iglsiter=1:repeat,
  S = sqrtm(GAMMAi);
  YS = S*Y;
  [W1,W2,PIvec,iteration,lambda]=marq(NetDef,W1,W2,PHI,YS(:,nmax+1:N),trparms);
  [Yhat,E] = nneval(NetDef,W1,W2,PHI,YS(:,nmax+1:N),1);
  E=(GAMMA*S')*E;
  GAMMA = (E*E')/(N-nmax);
  GAMMAi= inv(GAMMA);
end
