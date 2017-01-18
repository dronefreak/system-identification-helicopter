function [Yhat,PI]=ioleval(NetDeff,NetDefg,NN,W1f,W2f,W1g,W2g,Y,U)
%  NNIOL
%  -----
%          Evaluate a neural network trained by 'nniol'.
%
%  The function produces the following plots:
%    o Comparison of output and predicted output.
%    o Prediction error.
%    o Histogram showing the distribution of the prediction errors.
%
%  CALL:
%  [Yhat,NSSE]=ioleval(NetDeff,NetDefg,NN,W1f,W2f,W1g,W2g,Y,U)
%
%  INPUTS:
%  See function 'nniol' for an explanation of the inputs.
% 
%  OUTPUTS: 
%  Yhat    : One-step ahead prediction of the output.
%  NSSE    : Normalized sum of squared error (SSE/2N).

%  Written by : Magnus Norgaard, IAU/IMM, DTU
%  LastEditDate  : Jan 22, 2000

%----------------------------------------------------------------------------------
%--------------             NETWORK INITIALIZATIONS                   -------------
%----------------------------------------------------------------------------------

% >>>>>>>>>>>>>>>>>>>>>>>>   DETERMINE REGRESSOR STRUCTURE   <<<<<<<<<<<<<<<<<<<<<<   
[nu,N]      = size(U);             % # of inputs and # of data
na = NN(1);                        % # of past y's to be used in TDL
nb = NN(2:1+nu);                   % # of past u's to be used in TDL
nk = NN(2+nu:1+2*nu);              % Time delay in system
nmax        = max(na,nb+nk-1);
nab         = na+sum(nb);          % Number of "signal inputs" to each net
outputs  = 1;                      % # of outputs
inputs   = nab-1;                  % # of inputs


% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
PHI = zeros(nab-1,N-nmax);
jj  = nmax+1:N;
for k = 1:na, PHI(k,:)    = Y(jj-k); end
index = na;
for kk = 1:nu,                       % Here nu=1 allways!
for k = 1:nb(kk)-1, PHI(k+index,:) = U(kk,jj-k-nk(kk)); end
  index = index + nb(kk);
end
Y = Y(nmax+1:N);
U = U(nmax+1-nk:N-nk);
N = N-nmax;
PHI_aug   = [PHI;ones(1,N)];         % Augment PHI with a row containg ones


% >>>>>>>>>>>>>>>>>>>>      DETERMINE NETWORK STRUCTURE      <<<<<<<<<<<<<<<<<<<<<
% ---------- f-net structure ----------
hiddenf = length(NetDeff(1,:));         % Number of hidden neurons in f-net
L_hiddenf = find(NetDeff(1,:)=='L')';   % Location of linear hidden neurons
H_hiddenf = find(NetDeff(1,:)=='H')';   % Location of tanh hidden neurons
L_outputf = find(NetDeff(2,:)=='L')';   % Location of linear output neurons
H_outputf = find(NetDeff(2,:)=='H')';   % Location of tanh output neurons
y1f       =[zeros(hiddenf,N);ones(1,N)];% Hidden layer outputs
y2f       = zeros(outputs,N);           % Network output

% ---------- g-net structure ----------
hiddeng = length(NetDefg(1,:));         % Number of hidden neurons in g-net
L_hiddeng = find(NetDefg(1,:)=='L')';   % Location of linear hidden neurons
H_hiddeng = find(NetDefg(1,:)=='H')';   % Location of tanh hidden neurons
L_outputg = find(NetDefg(2,:)=='L')';   % Location of linear output neurons
H_outputg = find(NetDefg(2,:)=='H')';   % Location of tanh output neurons
y1g       =[zeros(hiddeng,N);ones(1,N)];% Hidden layer outputs
y2g       = zeros(outputs,N);           % Network output


%----------------------------------------------------------------------------------
%--------------                  SIMULATE NETWORK                     -------------
%----------------------------------------------------------------------------------

% >>>>>>>>>>>>>>>>>>>>> COMPUTE NETWORK OUTPUT  Yhat(theta)  <<<<<<<<<<<<<<<<<<<<<<
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

Yhat     = y2f + y2g.*U;                % Network output
E        = Y - Yhat;                    % Training error
E_vector = E(:);                        % Reshape E into a long vector
SSE      = E_vector'*E_vector;          % Sum of squared errors (SSE)
PI       = SSE/(2*N);                   % Performance index

si=figure-1;
plot(Y,'b-'); hold on
plot(Yhat,'r--');hold off
xlabel('time (samples)')
title('Output (solid) and one-step ahead prediction (dashed)')
grid;drawnow

figure(si+2)
subplot(211)
plot(E);
title('Prediction error (y-yhat)')
xlabel('time (samples)')
grid;

% Histogram over errors
subplot(212)
hist(E,20)
title('Histogram over prediction errors')
subplot(111)
drawnow
