function [X,Xscale]=dscale(X,Xscale)
% DSCALE
% ------
%         [Xs,Xscale]=dscale(X) scales data to zero mean and variance 1.
%
%         Xs=dscale(X,Xscale) scales data using the scaling parameters in 
%         Xscale: Xs(k,:) = [X(k,:) - Xscale(k,1)]/Xscale(k,2)
%         Typically used to scale a test data set with the same scaling parameters
%         that were used for scaling the training data to zero mean and variance 1.
%
% INPUTS:
%      X      - Data matrix.
%               (dimension is # of data vectors in matrix * # of data points)
%      Xscale - See below.
%
% OUTPUTS:
%      Xs     - Scaled data matrix.
%      Xscale - Matrix containing sample mean (column 1) and standard 
%               deviation (column 2) for each data vector in X.
%
% See the function WRESCALE on how to rescale the weights of the
% trained network.

% Written by Magnus Norgaard, IAU/IMM, Technical University of Denmark
% LastEditDate: Jan. 8 2000

if nargin==0,
    error('DSCALE called with no arguments.');
end
[r,N]  = size(X);              % r = # of data vectors, N = # of data
if nargin==1,
   Xscale = [mean(X')' std(X')']; % Col. 1 contains mean values, col 2 the std's
else
   if (size(Xscale,1)~=r | size(Xscale,2)~=2)
      error('Dimension mismatch between "X" and "Xscale"');
   end
end     
for k=1:r,
  X(k,:) = (X(k,:) - Xscale(k,1))/Xscale(k,2);
end
