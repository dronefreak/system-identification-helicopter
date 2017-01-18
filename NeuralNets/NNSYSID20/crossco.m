function coefs = crossco(v,w,maxlag)
% crossco    Calculate correlation coefficients.
%
% Cross-correlation coefficients:
%   coefs = crossco(v,w)
%   coefs = crossco(v,w,maxlag);
%
% Autocorrelation coefficients:
%   coefs = crossco(v,v,maxlag);
%
% v and w are two signals contained in vectors of equal length.
% Default max. lag is 25 or the vector length-1.

% Written by Magnus Norgaard, IAU/IMM, DTU
% LastEditDate Jan. 23, 2000
if nargin>3 | nargin<2
   error('Wrong number of arguments');
elseif nargin==2,
   maxlag=25;
end

% Extract mean and make sure the vectors are column vectors
v = v(:)-mean(v);
w = w(:)-mean(w);
if length(v) ~= length(w)
   error('v and w must have the same length');
end
maxlag = min(maxlag,length(v)-1);   % Reduce maxlag if vectors are too short

normcoef=sqrt(sum(v.*v)*sum(w.*w));
coefs = zeros(maxlag+1,1);            % Allocate vector for correlation function
for k=0:maxlag,
   coefs(k+1) = sum(v(1:end-k).*w(k+1:end))/normcoef;
end

coefs2 = zeros(maxlag,1);            % Allocate vector for correlation function
for k=1:maxlag,
   coefs2(k) = sum(w(1:end-k).*v(k+1:end))/normcoef;
end

coefs = [flipud(coefs2);coefs];

