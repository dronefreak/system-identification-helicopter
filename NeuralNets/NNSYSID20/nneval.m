function [y2,E,PI] = nneval(NetDef,W1,W2,PHI,Y,noplots)
% NNEVAL
% ------
%         Validation of ordinary feedforward neural networks.
%
% The predictions are compared to the true outputs, a histogram is shown
% for the prediction errors, and the autocorrelation coefficients for 
% the prediction error is plotted.
%
% CALL:
%        [Yhat,E,NSSE] = nneval(NetDef,W1,W2,PHI,Y)
%
%
% INPUTS:
%        See for example one of the functions MARQ, RPE, BATBP, INCBP
%
% OUTPUTS:
%         Yhat - Network predictions.
%         E    - Prediction errors.
%         NSSE - Normalized sum of squared errors.

% Written by : Magnus Norgaard, IAU/IMM Technical University of Denmark
% LastEditDate  : Jan. 23, 2000
[outputs,N] = size(Y);
[layers,dummy] = size(NetDef);        % Number of hidden layers
L_hidden = find(NetDef(1,:)=='L')';   % Location of linear hidden units
H_hidden = find(NetDef(1,:)=='H')';   % Location of tanh hidden units
L_output = find(NetDef(2,:)=='L')';   % Location of linear output units
H_output = find(NetDef(2,:)=='H')';   % Location of tanh output units
[hidden,inputs] = size(W1);
inputs   = inputs-1;
E        = zeros(outputs,N);
y1       = zeros(hidden,N);
y2       = zeros(outputs,N);


% ---  Compute network output ---
    h1 = W1*[PHI;ones(1,N)];  
    y1(H_hidden,:) = pmntanh(h1(H_hidden,:));
    y1(L_hidden,:) = h1(L_hidden,:);
    
    h2 = W2*[y1;ones(1,N)];
    y2(H_output,:) = pmntanh(h2(H_output,:));
    y2(L_output,:) = h2(L_output,:);

    E            = Y - y2;                     % Test error
    PI           = sum(sum(E.*E))/(2*N);       % Sum of squared errors 


if nargin~=6,
  % ---------- Output, Prediction and Prediction error ----------
  si=figure-1;
  for i=1:outputs
    figure(si+i);
    subplot(211)
    plot(Y(i,:),'bx'); hold on
    plot(y2(i,:),'ro');hold off
    if outputs==1,
      title('Observations = x     Network output=o')
    else
      title(['Observations = x     Network output=o      (output # ',num2str(i) ')'])
    end
    grid

    subplot(212)
    plot(E(i,:));
    title('Prediction error')
    grid
    subplot(111)
    drawnow
  end

  % Auto correlation function of error
  for i=1:outputs
    figure(si+outputs+i)
    subplot(211)
    M=min(25,N-1);
    Eauto=crossco(E(i,:),E(i,:),M);
    Eauto=Eauto(M+1:2*M+1);
    conf=1.96/sqrt(N);
    plot([0:M],Eauto(1:M+1),'b-'); hold on
    plot([0 M],[conf -conf;conf -conf],'r--');hold off
    set(gca,'xlim',[0 M]);
    xlabel('lag')
    if outputs==1,
      title('Autocorrelation coefficients for the prediction error')
    else
      title(['Autocorrelation coefficients for the prediction error (output # ',num2str(i) ')'])
    end
    grid
    subplot(111)
    drawnow


    % Histogram over errors
    subplot(212)
    hist(E(i,:),20)
    title('Histogram over prediction errors')
    subplot(111)
    drawnow
  end
end
