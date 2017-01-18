function [W1,W2]=wrescale(nntype,W1,W2,Uscale,Yscale,NN)
% WRESCALE
% --------
%         Rescale weights of a trained network so that the network
%         can work on unscaled data.
%
% Calling [W1,W2]=wrescale(method,W1,W2,Uscale,Yscale,NN) rescales the weights
% for networks with LINEAR OUTPUT UNITS. The function works for
% ordinary feed-forward networks as well as for input-output models of dynamic
% systems (i.e. NNAR(X), NNARMA(X) and NNOE type models). Notice that when
% the function is used on a pruned network it will re-introduce the biases removed
% during the pruning session.
%
% NN is left out for neural networks that are not models of dynamic systems.
%
% INPUTS:
%      method - The function applied for generating the model. For example
%               method='nnarx' or method='nnoe'. Use method='inverse' for 
%               inverse models.
%      W1     - Input-to-hidden of network trained on scaled data
%      W2     - Hidden-to-output weights
%      Uscale - Matrix containing sample mean and standard deviation for
%               each input. For time series an empty matrix, [], is passed. 
%      Yscale - Matrix containing mean and std's for each output.
%      NN     - Vector containing lag spaces (see nnarx, nnarmax, nnoe ..)
%
% OUTPUTS:
%      W1, W2 - "Unscaled" weight matrices.
%
% See the function DSCALE on how to scale the data before training.

% Written by Magnus Norgaard, IAU/IMM, Technical University of Denmark
% LastEditDate: Jan. 15, 2000

if strcmpi(nntype,'nniol')
   if ~iscell(W1) | ~iscell(W2)
      error('arguments W1 and W2 must cell structures when method is "nniol"');
   else
      W1f = W1{1}; W2f = W2{1};
      W1g = W1{2}; W2g = W2{2};
      [hiddenf,netinputs] = size(W1f);
      [hiddeng,netinputs] = size(W1g);
   end
else
   [hidden,netinputs] = size(W1);
end
[inputs,dummy]  = size(Uscale);
[outputs,dummy] = size(Yscale);

% ----- Rescale input-to-hidden weights -----
% -- Ordinary networks --
if strcmpi(nntype,'marq') | strcmpi(nntype,'marqlm') | strcmpi(nntype,'rpe') | ...
      strcmpi(nntype,'batbp') | strcmpi(nntype,'incbp') | strcmpi(nntype,'igls') 
   meanvec = Uscale(:,1);
   stdvec  = Uscale(:,2);
   
% -- Inverse models --
elseif strcmpi(nntype,'inverse') | strcmpi(nntype,'general')
  meanvec = repmat(Yscale(1),NN(1)+1,1);
  stdvec  = repmat(Yscale(2),NN(1)+1,1);
  meanvec = [meanvec;repmat(Uscale(1),NN(2)-1,1)];
  stdvec  = [stdvec;repmat(Uscale(2),NN(2)-1,1)];
  outputs = 1;
  
% -- Dynamic models --
else
  lagtot =  sum(sum(NN(1:size(NN,2)-size(Uscale,1))))+1;
  if (~strcmpi(nntype,'nnarmax1') & ~strcmpi(nntype,'nniol') & lagtot ~= netinputs) | ...
         (strcmpi(nntype,'nniol')  & lagtot ~= netinputs+1) ,
    error('You have made an error in the call of "wrescale"')
  end
  meanvec = Yscale(1)*ones(NN(1),1);
  stdvec  = Yscale(2)*ones(NN(1),1);
  nc = 0;
  
  % NARMAX2 model
  if strcmpi(nntype,'nnarmax2') | strcmpi(nntype,'nnramx2')
     nab=sum([NN(1:inputs) NN(inputs+1:inputs+outputs)]);
     nc = NN(inputs+outputs+1);
     
  % NARMAX model, linear noise filter
  elseif strcmpi(nntype,'nnarmax1') | strcmpi(nntype,'nnrarmx1')
     nab=sum([NN(1:inputs) NN(inputs+1:inputs+outputs)]);
     nc = NN(inputs+outputs+1);
     if nab==netinputs-1,
       nc = 0;
     end
  end
  
  % For all dynamic models
  for k=1:inputs,
    meanvec = [meanvec;Uscale(k,1)*ones(NN(k+1),1)];
    stdvec  = [stdvec;Uscale(k,2)*ones(NN(k+1),1)];
  end
  if nc>0,         % NARMAX models
    meanvec = [meanvec;zeros(nc,1)];
    stdvec  = [stdvec;Yscale(2)*ones(nc,1)];
  end
end

% ----- Rescale input-to-hidden weights -----
if strcmpi(nntype,'nniol')
   meanvec = meanvec(1:end-1);
   stdvec  = stdvec(1:end-1);
   for k=1:netinputs-1,
     W1f(:,k) = W1f(:,k)/stdvec(k);
     W1f(:,netinputs) = W1f(:,netinputs) - W1f(:,k)*meanvec(k);
     W1g(:,k) = W1g(:,k)/stdvec(k);
     W1g(:,netinputs) = W1g(:,netinputs) - W1g(:,k)*meanvec(k);
  end
else
  for k=1:netinputs-1,
     W1(:,k) = W1(:,k)/stdvec(k);
     W1(:,netinputs) = W1(:,netinputs) - W1(:,k)*meanvec(k);
  end
end

% ----- Rescale hidden-to-output weights -----
if strcmpi(nntype,'inverse') | strcmpi(nntype,'general')
   meanvec = Uscale(:,1);
   stdvec  = Uscale(:,2);
else
   meanvec = Yscale(:,1);
   stdvec  = Yscale(:,2);
end   
if strcmpi(nntype,'nniol')
      k=1;
      W2f(k,:) = W2f(k,:)*stdvec(k);
      W2f(k,hiddenf+1) = W2f(k,hiddenf+1)+meanvec(k);
      W2g(k,:) = W2g(k,:)*stdvec(k)/Uscale(2);
      
      W1f = [W1f;W1g];
      W2f = [W2f(1:hiddenf) -W2g(1:hiddeng)*Uscale(1) W2f(hiddenf+1)-W2g(hiddeng+1)*Uscale(1)];
    %  NetDeff = [repmat('H',1,hiddenf+hiddeng);'L' repmat('-',1,hiddenf+hiddeng-1)];
      W1 = {W1f,W1g};
      W2 = {W2f,W2g};
else
  for k=1:outputs,
    W2(k,:) = W2(k,:)*stdvec(k);
    W2(k,hidden+1) = W2(k,hidden+1)+meanvec(k);
  end
end
