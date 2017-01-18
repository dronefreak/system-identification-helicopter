% DEMONSTRATION PROGRAM FOR TESTING NNARMAX2
%
% Programmed by Magnus Norgaard, IAU/IMM/EI, Technical Univ. of Denmark
% LastEditDate: Jan 15, 2000

close all
StopDemo=0;
figure
guihand=gcf;
for k=1:1, %dummy loop

  % >>>>>>>>>>>>>>>>  BUILD GUI INTERFACE  <<<<<<<<<<<<<<<<<
  [guihand,edmulti,contbut,quitbut]=pmnshow;
  set(guihand,'Name','NNARMAX demonstration');

  % >>>>>>>>>>>>>>>>  SCREEN 1  <<<<<<<<<<<<<<<<<
  s0='1';
  s1='The purpose of this demo is to show how the NNARMAX2';
  s2='function can be used for modelling a second order';
  s3='nonlinear dynamic system.';
  s4=[];
  s5='To generate the model we will use the noise corrupted';
  s6='data set shown above.';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6);

  load spmdata
  subplot(411)
  plot(u1);
  title('Input and output sequence')
  subplot(412)
  plot(y1);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  % >>>>>>>>>>>>>>>>  SCREEN 2  <<<<<<<<<<<<<<<<<
  s0='2';
  s1='First we have to select a model structure inside';
  s2='which we wish to search for a good model.';
  s3=[];
  s4='The model structure selection consists of two';
  s5='subproblems: Choosing a regressor structure and';
  s6='choosing a network architecture.';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end


  % >>>>>>>>>>>>>>>>  SCREEN 3  <<<<<<<<<<<<<<<<<
  subplot(411);delete(gca);subplot(412);delete(gca)
  subplot('position',[0.1 0.5 0.6 0.5])
  W1=rand(5,7);
  W2=rand(1,6);
  drawnet(W1,W2,eps,{'y(t-1)' 'y(t-2)' 'u(t-1)' 'u(t-2)' 'e(t-1)' 'e(t-2)'},{'yhat(t)'});
  s0='3';
  s1='As regressors we will use two past inputs, two';
  s2='past outputs, and two past residuals.';
  s3=[];
  s4='Furthermore we will choose a network architecture';
  s5='with 5 hidden ''tanh'' units and one linear output unit.';
  smat=str2mat(s0,s1,s2,s3,s4,s5);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  

  % >>>>>>>>>>>>>>>>  SCREEN 4  <<<<<<<<<<<<<<<<<
  s0='4';
  s1='Now that a model structure has been selected, we are';
  s2='ready to begin training.';
  s3=[];
  s4='Let''s run the NNARMAX2 function, which uses a Levenberg-';
  s5='Marquardt algorithm for generating a NNARMAX-model.';
  smat=str2mat(s0,s1,s2,s3,s4,s5);
  
  % ----- Choose training parameters -----
  trparms = settrain;
  trparms = settrain(trparms,'maxiter',100,'D',1e-3,'skip',10);
  NetDef = ['HHHHH'
            'L----']; 
  NN = [2 2 2 1];
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end


  % >>>>>>>>>>>>>>>>  SCREEN 5  <<<<<<<<<<<<<<<<<
  % ----- Train network -----
  s0='5';
  s1=[];
  s2='    >> Training process in action!! <<';
  s3=[];
  s4=[];
  s5='We run up to 100 iterations so you may have to';
  s6='wait for a while.';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6);
  set(edmulti,'String',smat);
  drawnow
  [W1,W2,NSSEvec,iter,lambda]=nnarmax2(NetDef,NN,[],[],trparms,y1,u1);
  delete(gca);
  subplot('position',[0.1 0.55 0.45 0.38]);
  semilogy(NSSEvec);
  xlabel('Iteration');
  ylabel('Criterion of fit');
  grid
  

  % >>>>>>>>>>>>>>>>  SCREEN 6  <<<<<<<<<<<<<<<<<
  % ----- Evaluate network -----
  s0='6';
  s1='The network has now been trained and we conclude';
  s2='the session by validating the model on a fresh';
  s3='data set.';
  smat=str2mat(s0,s1,s2,s3);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end 
  [Yhat,NSSE]=nnvalid('nnarmax2',NetDef,NN,W1,W2,y2,u2);
  figure(1)
  
  % >>>>>>>>>>>>>>>>  SCREEN 7  <<<<<<<<<<<<<<<<<
  s1=[];
  s2='                  >> THE END <<';
  smat=str2mat(s1,s1,s1,s1,s2);
  set(edmulti,'String',smat);
  drawnow
end
