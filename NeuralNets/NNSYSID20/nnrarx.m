function [W1,W2,PI_vector,iteration]=nnrarx(NetDef,NN,W1,W2,trparms,Y,U)
%  NNRARX
%  ------
%          Determine a nonlinear ARX model of a dynamic system by training a
%          two layer neural network with a recursive Gauss-Newton method.
%          The function can handle multi-input systems (MISO).
%
%  CALL:
%          [W1,W2,critvec,iteration,lambda]=...
%                            nnrarx(NetDef,NN,W1,W2,trparms,Y,U)
%
%  INPUTS: 
%  U       : Input signal (= control signal) (left out in the nnarma case)
%            dim(U) = [(inputs) * (# of data)]
%  Y       : Output signal. dim(Y) = [1 * # of data]
%  NN      : NN=[na nb nk].
%            na = # of past outputs used for determining the prediction
%            nb = # of past inputs used for determining prediction
%            nk = time delay (usually 1)
%            For multi-input systems nb and nk contain as many columns as
%            there are inputs. 
%  W1,W2   : Input-to-hidden layer and hidden-to-output layer weights.
%            If they are passed as [] they are initialized automatically
%  
%           For time series (NNAR models), NN=na only.
%
%  See the function RPE for an explanation of the remaining inputs
%  as well as of the returned variables.
%                                                                              
%  Programmed by : Magnus Norgaard, IAU/IMM, Technical University of Denmark
%  LastEditDate  : Jan. 4, 2000

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   INITIALIZATIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<   
N  = length(Y);
na = NN(1);
if length(NN)==1
  nb = 0;                          % nnar model
  nk = 0;
  nu = 0; 
else
  [nu,N]      = size(U); 
  nb = NN(2:1+nu);                 % nnarx model
  nk = NN(2+nu:1+2*nu);
end

nmax        = max(na,nb+nk-1);
nab         = na+sum(nb); 

% -- Initialize weights if nescessary --
if isempty(W1) | isempty(W2),
  hidden = length(NetDef(1,:));    % Number of hidden neurons
  W1 = rand(hidden,nab+1)-0.5;
  W2 = rand(1,hidden+1)-0.5;
end


% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
PHI = zeros(nab,N-nmax);
jj  = nmax+1:N;
for k = 1:na, PHI(k,:)    = Y(jj-k); end
index = na;
for kk = 1:nu,
  for k = 1:nb(kk), PHI(k+index,:) = U(kk,jj-k-nk(kk)+1); end
  index = index + nb(kk);
end


% >>>>>>>>>>>>>>>>>>>>         CALL TRAINING FUNCTION         <<<<<<<<<<<<<<<<<<<<<
[W1,W2,PI_vector,iteration]=rpe(NetDef,W1,W2,PHI,Y(nmax+1:N),trparms);









