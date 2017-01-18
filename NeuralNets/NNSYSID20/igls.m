function [W1,W2,GAMMA,lambda]=igls(NetDef,W1,W2,trparms,GAMMA,PHI,Y);
%  IGLS
%  ----
%          Train a multi-output network and estimate the covariance matrix
%          simultaneously using the IGLS-procedure (iterated Generalized
%          Least Squares). The network is trained with a Levenberg-Marquardt
%          method. So far, the function is restriced to work for NNARX and
%          NNSSIF models.
%
%  CALL:
%       [W1,W2,Gamma,lambda]=igls(NetDef,W1,W2,trparms,Gamma0,PHI,Y);
%
%  INPUTS: 
%    NN,W1,W2,PHI,Y : See MARQ
%    trparms:  Data structure with parameters associated with the
%              training algorithm (optional). Use the function SETTRAIN if
%              you do not want to use the default values.
%    Gamma0  : Covariance matrix. If passed as [] it is initialized to 
%              the identity matrix.
%
%  OUTPUTS:
%    W1, W2, lambda: See MARQ
%    Gamma         : Estimated covariance matrix

%  Programmed by : Magnus Norgaard, IAU/IMM Technical University of Denmark
%  LastEditDate  : Jan 8, 2000

if isempty(trparms) % Default training parameters
   trd = settrain(trparms,'repeat','default');
   repeat=trd.repeat;
else
  if ~isfield(trparms,'repeat')
     trparms = settrain(trparms,'repeat','default');
  end
  repeat    = trparms.repeat;
end
[ny,N] = size(Y);
if isempty(GAMMA), GAMMA=eye(ny); end
GAMMAi = inv(GAMMA);
for iglsiter=1:repeat,
  S = sqrtm(GAMMAi);
  YS= S*Y;
  [W1,W2,PIvec,iteration,lambda]=marq(NetDef,W1,W2,PHI,YS,trparms);
  [Yhat,E] = nneval(NetDef,W1,W2,PHI,YS,1);
  E=(GAMMA*S')*E;
  GAMMA = (E*E')/N;
  GAMMAi= inv(GAMMA);
end
