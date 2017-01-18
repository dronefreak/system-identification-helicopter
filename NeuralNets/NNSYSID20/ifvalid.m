function [Yhat,PI]=ifvalid(NetDef,nx,W1,W2,obsidx,Y,U)
%  IFVALID
%  -------
%          Validate a neural network based state space innovations form model
%          of a (possibly multivariable) dynamic system. I.e., a network model
%          trained with the function NNSSIF.
%
%          The following plots are produced: 
%          - Output(s) together with predicted output(s)
%          - Prediction error
%          - Autocorrelation coefficients for prediction error and cross-
%            correlation between prediction error(s) and input(s)
%          - Histogram(s) showing the distribution of the prediction errors
%          - Coefficients of extracted linear models.
%
%           Call: 
%           [Yhat,NSSE] = ifvalid(NetDef,nx,W1,W2,obsidx,Y,U)               

%  Written by : Magnus Norgaard, IAU/IMM
%  LastEditDate  : Jan 23, 2000  



% >>>>>>>>>>>>>>>>>>>>>>>>>>>>     INITIALIZATIONS     <<<<<<<<<<<<<<<<<<<<<<<<<<<< 
skip=1;
[ny,Ndat]= size(Y);                     % # of outputs and # of data
[nu,Ndat]= size(U);                     % # of inputs 
N        = Ndat - 1;                    % Size of complete training set
N2       = N-skip+1;                    % "Actual" size of training set
hidden   = length(NetDef(1,:));         % Number of hidden neurons
inputs   = nx+nu+ny;                    % Number of inputs to the network
outputs  = ny;                          % Network outputs equal to the number of states 
L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
y1       = zeros(hidden,N);             % Hidden layer outputs
y1       = [y1;ones(1,N)];
y2       = zeros(outputs,N);            % Network output
Yhat     = zeros(ny,N);                 % Initialize prediction matrix
E        = zeros(ny,N);                 % Initialize prediction error matrix
E_new    = zeros(ny,N);                 % Initialize prediction error matrix
if isempty(obsidx),                     % Choose observability indices if necessary
  rowidx   = [1:ny-1 nx];
  obsidx = [rowidx(1) rowidx(2:ny)-rowidx(1:ny-1)]; 
else                                    % Otherwise compute the row indices
  obsidx=obsidx(:)';
  rowidx=obsidx;
  for k=2:ny,
    rowidx(k)=obsidx(k)+rowidx(k-1);
  end
end
nrowidx = 1:nx;                         % Not row indices
nrowidx(rowidx)=[];
Cidx=[1 rowidx(1:ny-1)+1];
C = zeros(ny,nx);
C(1:ny,Cidx)=eye(ny);
NNy      = N*ny;
skipidx  = ny*(skip-1)+1:NNy;



% >>>>>>>>>>>>>>>>>>>>  CONSTRUCT THE REGRESSION MATRIX PHI   <<<<<<<<<<<<<<<<<<<<<
PHI = zeros(inputs,N);
PHI(nx+1:nx+nu,:) = U(:,1:N);
PHI_aug = [PHI;ones(1,N)];                % Augment PHI with a row containing ones
Y       = Y(:,2:Ndat);                    % Extract the 'target' part of Y



% >>>>>>>>>>>>>>>>>>>>>>>>>>   COMPUTE NETWORK OUTPUT   <<<<<<<<<<<<<<<<<<<<<<<<<<<
for t=1:N,
  h1 = W1*PHI_aug(:,t);                 % Hidden neuron outputs
  y1(H_hidden,t) = pmntanh(h1(H_hidden));
  y1(L_hidden,t) = h1(L_hidden);    

  h2 = W2*y1(:,t);                      % Predicted states
  y2(H_output,t) = pmntanh(h2(H_output,:));
  y2(L_output,t) = h2(L_output,:);
  y2(nrowidx,t)  = y2(nrowidx,t) + PHI_aug(nrowidx+1,t);
  Yhat(:,t)      = C*y2(:,t);

  E(:,t) = Y(:,t) - Yhat(:,t);          % Prediction error
  Evec   = E(:);
  for d=1:min(1,N-t),
    PHI_aug(1:nx,t+1) = y2(:,t);
    PHI_aug(nx+nu+1:inputs,t+1) = E(:,t);
  end
end
SSE=Evec(skipidx)'*Evec(skipidx);       % Sum of squared errors (SSE)
PI       = SSE/(2*N2);                  % Performance index



% >>>>>>>>>>>>>>>>>>>>>>>>>>      PLOT THE RESULTS      <<<<<<<<<<<<<<<<<<<<<<<<<<<
si = figure-1;

% ---------- Output, Prediction and Prediction error ----------
for i=1:outputs
    figure(si+i)
    subplot(211)
    plot(Y(i,:),'b-'); hold on
    plot(Yhat(i,:),'r--');hold off
    xlabel('time (samples)')
    if outputs==1,
      title('Output (solid) and one-step ahead prediction (dashed)')
    else
      title([' Output (solid) and one-step ahead prediction (dashed) (output # ',num2str(i) ')'])
    end
    grid

    subplot(212)
    plot(E(i,:)');
    title('Prediction error (y-yhat)')
    xlabel('time (samples)')
    grid
    subplot(111)
    drawnow
end


% --------- Correlation functions ----------
for ii=1:outputs,
  figure(si+outputs+ii)
  subplot(nu+1,1,1);
  M=min(25,N-1);
  Eauto=crossco(E(ii,:),E(ii,:),M);
  Eauto=Eauto(M+1:2*M+1);
  conf=1.96/sqrt(N);
  plot([0:M],Eauto(1:M+1),'b-'); hold on
  plot([0 M],[conf -conf;conf -conf],'r--');hold off
  set(gca,'xlim',[0 M]);
  xlabel('lag')
  title(['Autocorrelation coefficients for prediction error #' num2str(ii)])
  grid
  
  for i=1:nu,
    subplot(nu+1,1,i+1);
    UEcross=crossco(E(ii,:),U(i,1:N),M);
    plot([-M:M], UEcross,'b-'); hold on
    plot([-M M],[conf -conf;conf -conf],'r--');hold off
    xlabel('lag')
    title(['Cross-correlation coefficients of u' num2str(i) ' and prediction error'])
    ymax=min(5*conf,max([abs(UEcross)]));
    ymax=max(1.2*conf,ymax);
    axis([-M M -ymax ymax]);
    grid
  end
  subplot(111)
  drawnow
end

% --------- Histograms over the prediction errors ---------
figure(si+2*outputs+1)
for ii=1:outputs,
  subplot(outputs,1,ii);
  hist(E(ii,:),20)
  title(['Histogram over prediction errors on output #' num2str(ii)])
end
subplot(111)

% ---------- Extract linear model from network ----------
    % ---------- Find derivative of output wrt. the past residuals ----------
    Amat=zeros(ny*nx,N);
    Kmat=zeros(nx*ny,N);
    Bmat=zeros(nx*nu,N);
    for t=1:N,
      dxdy1 = W2(:,1:hidden);
      for j = H_output',
        dxdy1(j,:) = W2(j,1:hidden)*(1-Yhat(j,t).*Yhat(j,t));
      end

      % Matrix with partial derivatives of the output from each hidden neurons with
      % respect to each input:
      dy1dx = W1(:,1:nx);
      dy1de = W1(:,nx+nu+1:inputs);
      dy1du = W1(:,nx+1:nx+nu);
      for j = H_hidden',
        dy1dx(j,:) = W1(j,1:nx)*(1-y1(j,t).*y1(j,t));
        dy1de(j,:) = W1(j,nx+nu+1:inputs)*(1-y1(j,t).*y1(j,t));
        dy1du(j,:) = W1(j,nx+1:nx+nu)*(1-y1(j,t).*y1(j,t));
      end

      % Matrix with partial derivative of each output with respect to each input
      Ahat = dxdy1 * dy1dx;
      Acoef = Ahat(rowidx,:);
      Bhat = dxdy1 * dy1du;
      Khat = dxdy1 * dy1de;
      Amat(:,t) = Acoef(:);
      Bmat(:,t) = Bhat(:);
      Kmat(:,t) = Khat(:);
    end

  figure(si+2*outputs+2)
  subplot(311)
  plot(Amat')
  title('Extracted A matrix coefficients')
  grid
  
  subplot(312)
  plot(Bmat')
  title('Extracted B matrix coefficients')
  grid
  
  subplot(313)
  plot(Kmat')
  title('Extracted K matrix coefficients')
  xlabel('time (samples)')
  grid
  subplot(111)
figure(si+1)

