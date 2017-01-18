function [OrderIndexMat]=lipschit(U,Y0,mvec,nvec)
% LIPSCHIT
% --------
%            Determine the lag space.
%
%            Given a set of corresponding inputs and outputs the
%            function calculates a matrix of indices, which can be
%            helpful when trying to determine a proper lag space structure
%            (m and n) before identifying a model of a dynamic system:
%                 y(t) = f(y(t-1),...,y(t-n), u(t-1),..., u(t-m))
%
%            An insufficient lag space structure leads to a large index.
%            While increasing the lag space, the index will decrease until
%            a sufficiently large lag space structure is reached. Increasing
%            the lag space beyond this will not reduce the index significantly.
%            In other words: look for theknee-point where the order index 
%            flattens out.
%
%            NB: So far, the function works for SISO systems only.
%
%            Please ignore the message "Divide by zero"
%
% CALL:
%  [OrderIndexMat]=lipschit(U,Y,m,n)
%  
%  m is a vector specifying which input lag spaces to investigate and n is
%  ditto for the output. If one is only interested in the order index for one
%  particular choice of lag structure, n and m are specified as scalars, and
%  only the order index is returned. In the more general case, where one or
%  both are vectors, the function will also produce one or two plots.
%
%  Examples of some typical special cases:
%  o  NNFIR model structure expected:
%     m=[1:20]; n=0;
%
%  o  Time series:
%     U=[]; m=0;
%
%  o  Check unly n=m:
%     m=[1:5]; n=m;
% 
% INPUTS:
% U  - Sequence of inputs  (row vector)
% Y  - Sequence of outputs (row vector)
% m  - Vector specifying the input lag spaces to investigate
% n  - Vector specifying the ouput lag spaces to investigate
%
% OUPUTS:
% OrderIndexMat - A matrix containing the order indices for each combination
%                 of elements in the vectors m and n. The number of rows
%                 corresponds to the number of elements in m, while the
%                 number of columns corresponds to the number of elements in n
%
% REFERENCE:
% X. He & H. Asada: "A New Method for Identifying Orders of Input-Output
%                    Models for Nonlinear Dynamic Systems"
%            Proc. of the American Control Conf., S.F., California, 1993
%
% SEE ALSO: Use function 'dscale' to scale the data.
%
% Programmed by Magnus Norgaard, IAU/IMM, Technical Univ. of Denmark
% LastEditDate: July 16, 1996

% ---------- Generate matrix of regressors, PHI -----------
Ndat = length(Y0);
mvec=sort(mvec); nvec=sort(nvec);
nu   = 1;                      % Inputs (only one is considered here)
nk   = 1;                      % Time delay (only one is considered here)
OrderIndexMat = zeros(length(mvec),length(nvec));
mi=1;
for m=mvec,
  ni=1;
  for n=nvec,
l    = n+m;                    % Rows in phi
nmax = max(n,m+nk-1);
N    = Ndat-nmax;              % # of regressor vector - output pairs
p    = max(floor(0.02*N),min(N,10)); % # of coefficients used to determine order index
np   = N*p;
PHI  = zeros(l,N);
jj   = nmax+1:Ndat;
for k = 1:n, PHI(k,:) = Y0(jj-k); end
index = n;
for kk = 1:nu,
  for k = 1:m(kk), PHI(k+index,:) = U(kk,jj-k-nk(kk)+1); end
  index = index + m(kk);
end


% ---------- Eliminate first elements from Y ----------
Y = Y0(nmax+1:Ndat);


% ---------- Initialize different matrices ----------
DY   = zeros(1,N-1);   % For saving differenced outputs
DPHI = zeros(l,N-1);   % For saving differenced regressor vectors
Q    = zeros(N,p);     % For saving a large number of Lipschitz coef.
q    = zeros(1,N-1);
q2   = q;
DPHI2=q;
onesvec = ones(1,N-1); % Vector of ones


% ---------- Compute the Lipschitz coefficients ----------
for i=2:N,
  DY(1:i-1)     = Y(i)-Y(1:i-1);
  DPHI(:,1:i-1) = PHI(:,i(onesvec(1:i-1)))-PHI(:,1:i-1);
  if l>1,
    DPHI2         = sum(DPHI(:,1:i-1).*DPHI(:,1:i-1));
  else
    DPHI2         = DPHI(:,1:i-1).*DPHI(:,1:i-1);
  end
  q2(1:i-1)     = DY(1:i-1).*DY(1:i-1)./DPHI2(1:i-1);
  q(1:i-1)   = -sort(-q2(1:i-1));
  maxindex      = min(i-1,p);             % Max index to non-zero elements
  if q(maxindex)==inf | q(maxindex)==NaN, % Remove possible NaN's and Inf's i
    finitevec = find(finite(q(1:maxindex)));
    maxindex  = finitevec(length(finitevec));
  else
    Q(i,1:maxindex) = sqrt(q(1:maxindex));% Compute square root and save p largest
  end
end


% ---------- Determine the p largests coefficients ----------
PMaxLipCoef = (sort(Q(:)));
PMaxLipCoef = PMaxLipCoef(np-p+1:np);


% ---------- Finally determine the order index ----------
OrderIndex = prod(sqrt(l)*PMaxLipCoef)^(1/p);

OrderIndexMat(mi,ni)=OrderIndex;
    ni = ni+1;
  end
  mi = mi+1;
end


% ---------- Plot the order indices ----------
if length(mvec)>1 & length(nvec)>1,
  figure(1)
  surf(mvec,nvec,OrderIndexMat');
  set(gca,'Zscale','log');
  set(gca,'YTickLabels',nvec)
  set(gca,'YTick',nvec)
  set(gca,'XTickLabels',mvec)
  set(gca,'XTick',mvec)
  view([-600 30])
  xlabel('Number of past inputs')
  ylabel('Number of past outputs')
  zlabel('Order index')
  title('Order index vs. lag space')
  colormap jet
  grid
  if length(mvec)==length(nvec)
    if (mvec-nvec)==0,
      figure(2)
      semilogy(nvec,diag(OrderIndexMat));
      set(gca,'XTickLabels',nvec)
      set(gca,'XTick',nvec)
      grid
      xlabel('Number of past inputs and outputs');
      ylabel('Order index');
      title('Order index vs. lag space')
    end
  end
  
elseif length(mvec)>1 & length(nvec)==1,
  semilogy(mvec,OrderIndexMat);
  set(gca,'XTickLabels',mvec)
  set(gca,'XTick',mvec)
  grid
  xlabel('Number of past inputs');
  ylabel('Order index');
  title('Order index vs. lag space')
  
elseif length(mvec)==1 & length(nvec)>1,
  semilogy(nvec,OrderIndexMat);
  set(gca,'XTickLabels',nvec)
  set(gca,'XTick',nvec)
  grid
  xlabel('Number of past outputs');
  ylabel('Order index');
  title('Order index vs. lag space')
end


