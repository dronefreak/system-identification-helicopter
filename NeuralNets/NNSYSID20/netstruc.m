function [W1,W2]=netstruc(NetDef,thd,index)
% NETSTRUC
% --------
%          [W1,W2]=netstruc(NetDef,thd,index) extracts the
%          weight matrices from the matrix of parameter vectors
%          produced by the pruning functions OBDPRUNE, OBSPRUNE
%          and NNPRUNE. The variable 'index' specifies the location 
%          in 'thd' where the optimal parameter vector is located.
%
%  Programmed by : Magnus Norgaard, IAU/IMM
%  LastEditDate  : July 16, 1997
theta=thd(:,index);
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of hyperbolic tangent hidden neurons
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of hyperbolic tangent output neurons
outputs  = max([L_output;H_output]);
hidden   = max([L_hidden;H_hidden]);
parameters2 = outputs*(hidden+1);       % # of hidden-to-output weights
parameters1 = length(theta)-parameters2;% # of input-to-hidden weights
parameters  = parameters1 + parameters2;% Total # of weights
inputs      = parameters1/hidden - 1;
W1 = reshape(theta(parameters2+1:parameters),inputs+1,hidden)';
W2 = reshape(theta(1:parameters2),hidden+1,outputs)';
drawnet(W1,W2,eps);
