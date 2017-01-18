% DEMONSTRATION PROGRAM FOR ILLUSTRATING THE EFFECT OF REGULARIZATION
%
% Written by Magnus Norgaard, IAU/IMM, Technical Univ. of Denmark
% LastEditDate: Jan. 15, 2000.

close all
StopDemo=0;
figure
guihand=gcf;
for k=1:1, %dummy loop

  % >>>>>>>>>>>>>>>>  BUILD GUI INTERFACE  <<<<<<<<<<<<<<<<<
  [guihand,edmulti,contbut,quitbut]=pmnshow;
  set(guihand,'Name','Demonstration of regularization');

  % >>>>>>>>>>>>>>>>  SCREEN 1  <<<<<<<<<<<<<<<<<
  s0='1';
  s1='In this demo it is shown how regularization by simple';
  s2='weight decay might be of help when dealing with';
  s3='overparametrization in neural networks.';
  s4=[];
  s5='The problem, which will be the subject of our';
  s6='investigation, is to use a neural network for fitting';
  s7='the underlying sine wave from the points marked';
  s8='''training data.''';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6,s7,s8);

  % ---------- Generate Data ----------
  randn('seed',1);
  PHI=2*pi*rand(1,600);
  Y=sin(PHI)+0.2*randn(1,length(PHI));
  PHI1 = PHI(:,1:300);
  Y1   = Y(:,1:1:300);
  PHI2 = PHI(:,301:600);
  Y2   = Y(:,301:600);

  %----- Plot Data -----
  sub1=subplot('position',[0.1 0.55 0.38 0.38]);
  plot(PHI1,Y1,'+');
  set(gca,'Xlim',[min(PHI1) max(PHI1)]);
  title('Training data');
  sub2=subplot('position',[0.57 0.55 0.38 0.38]);
  plot(PHI2,Y2,'m+')
  set(gca,'Xlim',[min(PHI2) max(PHI2)]);
  title('Test data');
  drawnow
  
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 2  <<<<<<<<<<<<<<<<<
  s0='2';
  s1='Let''s begin by training a network with 15';
  s2='hidden ''tanh'' units and one linear output unit';
  s3='without using regularization.';
  smat=str2mat(s0,s1,s2,s3);
  NetDef = ['HHHHHHHHHHHHHHH'       
            'L--------------'];
  W1  = rand(15,2);             % Weights to hidden layer 
  W2  = rand(1,16);             % Weights to output
  delete(sub1);
  delete(sub2);
  sub1=subplot('position',[0.1 0.55 0.45 0.38]);
  drawnet(W1,W2,eps,{'x'},{'y'});
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 3  <<<<<<<<<<<<<<<<<
  % ----- Train network -----
  s0='3';
  s1=[];
  s2='    >> Training process in action!! <<';
  s3=[];
  s4=[];
  s5='We run up to 500 iterations so you may have to';
  s6='wait for a while!';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6);
  set(edmulti,'String',smat);
  drawnow
  trparms = settrain;
  trparms = settrain(trparms,'maxiter',500);
  [W1,W2,NSSEvec,iter,lambda2]=marq(NetDef,W1,W2,PHI1,Y1,trparms);
  delete(gca);
  subplot('position',[0.1 0.55 0.45 0.38]);
  semilogy(NSSEvec);
  xlabel('Iteration');
  ylabel('Criterion of fit');
  grid
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 4  <<<<<<<<<<<<<<<<<
  s0='4';
  s1='Next we compute training error, test error,';
  s2='FPE estimate, and LOO estimate to get an idea';
  s3='of how well the network fits the curve.';
  smat=str2mat(s0,s1,s2,s3);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end


  % >>>>>>>>>>>>>>>>  SCREEN 5  <<<<<<<<<<<<<<<<<
  [Yhat,E,NSSE_tr] = nneval(NetDef,W1,W2,PHI1,Y1,1);
  [Yhat,E,NSSE_te] = nneval(NetDef,W1,W2,PHI2,Y2,1);
  FPE = fpe(NetDef,W1,W2,PHI1,Y1,trparms);
  trparms2 = settrain(trparms,'maxiter',0);
  ELOO= loo(NetDef,W1,W2,PHI1,Y1,trparms2);
  s0='5';
  s1=['Training error:   ' num2str(NSSE_tr)];
  s2=['Test Error:         ' num2str(NSSE_te)];
  s3=['FPE estimate:   ' num2str(FPE)];
  s4=['LOO estimate:   ' num2str(ELOO)];
  smat=str2mat(s0,s1,s2,s3,s4);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 6  <<<<<<<<<<<<<<<<<
  s0='6';
  s1='This result is typical for networks having too many';
  s2='weights. The superflous weights will capture';
  s3='some of the noise on the training set, leading to';
  s4='a poor generalization ability. This phenomenon is';
  s5='usually referred to as ''overfitting''';
  s6=[];
  s7='Let''s try to train with a small weight decay (0.02)';
  s8='and see what happens.';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6,s7,s8);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end 
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 7  <<<<<<<<<<<<<<<<<
  % ----- Train network -----
  s0='7';
  s1=[];
  s2='    >> Training process in action!! <<';
  s3=[];
  s4=[];
  s5='We run up to 500 iterations so you may have to';
  s6='wait for a while!';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6);
  set(edmulti,'String',smat);
  drawnow
  trparms = settrain(trparms,'D',0.02);
  [W1,W2,NSSEvec,iter,lambda2]=marq(NetDef,W1,W2,PHI1,Y1,trparms);
  delete(gca);
  subplot('position',[0.1 0.55 0.45 0.38]);
  semilogy(NSSEvec);
  xlabel('Iteration');
  ylabel('Criterion of fit');
  grid
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 8  <<<<<<<<<<<<<<<<<
  [Yhat,E,NSSE_tr2] = nneval(NetDef,W1,W2,PHI1,Y1,1);
  [Yhat,E,NSSE_te2] = nneval(NetDef,W1,W2,PHI2,Y2,1);
  FPE2 = fpe(NetDef,W1,W2,PHI,Y,trparms);
  trparms2 = settrain(trparms,'maxiter',0);
  ELOO2= loo(NetDef,W1,W2,PHI,Y,trparms2);
  s0='8';
  s1='                     No regularization     Regularization';
  s2=['Training error:   ' num2str(NSSE_tr) '               ' num2str(NSSE_tr2)  ];
  s3=['Test Error:         ' num2str(NSSE_te) '              ' num2str(NSSE_te2)];
  s4=['FPE estimate:   ' num2str(FPE) '                ' num2str(FPE2)];
  s5=['LOO estimate:   ' num2str(ELOO) '                ' num2str(ELOO2)];
  s6=[];
  s7='It looks as if the weight decay actually improved';
  s8='generalization.';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6,s7,s8);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end 
  

  % >>>>>>>>>>>>>>>>  SCREEN 9  <<<<<<<<<<<<<<<<<
  s0='9';
  s1='To really proof the effect of regularization, we';
  s2='redo the experiment for 50 different values of the';
  s3='weight decay parameter. Also we train the network';
  s4='7 times for each weight decay, using different initial';
  s5='weights in order to reduce the influence from local';
  s6='minima.';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end 


  % >>>>>>>>>>>>>>>>  SCREEN 10  <<<<<<<<<<<<<<<<<
  s0='10';
  s1='Well OK I think we cheat by simply loading the';
  s2='results from a data file.';
  smat=str2mat(s0,s1,s2);
  load test6mat
  semilogx(D_vec,data1,'x',D_vec,data2,'o')
  hold on;plot([0.0339 0.0339],[0.015 0.045],'r--');hold off
  xlabel('Weight decay parameter')
  ylabel('Normalized SSE')
  title('x = training data,       o = test data')
  set(gca,'Ylim',[0.015 0.0339]);
  set(gca,'Xlim',[1e-6 1]);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 11  <<<<<<<<<<<<<<<<<
  s0='11';
  s1='Two important remarks can be made from this';
  s2='experiment:';
  s25=[];
  s3='1) When training an overparametrized network on';
  s4='    noisy data, regularization improves generalization.';
  s5='2) Regularization has a smoothing effect on the';
  s6='    criterion. This significantly reduces the number of';
  s7='    local minima.';
  smat=str2mat(s0,s1,s2,s25,s3,s4,s5,s6,s7);
  set(edmulti,'String',smat);
  drawnow
end
