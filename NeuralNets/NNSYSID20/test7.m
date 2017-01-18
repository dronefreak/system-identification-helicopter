% DEMONSTRATION PROGRAM FOR ILLUSTRATING PRUNING BY OBS
%
% Programmed by Magnus Norgaard, IAU/IMM/EI, Technical Univ. of Denmark
% LastEditDate: Aug 21, 1995.

close all
StopDemo=0;
figure
guihand=gcf;
for k=1:1, %dummy loop

  % >>>>>>>>>>>>>>>>  BUILD GUI INTERFACE  <<<<<<<<<<<<<<<<<
  [guihand,edmulti,contbut,quitbut]=pmnshow;
  set(guihand,'Name','Demonstration of pruning');


  % >>>>>>>>>>>>>>>>  SCREEN 1  <<<<<<<<<<<<<<<<<
  s0='1';
  s1='The figure shows the famous sunspot benchmark data,';
  s2='well-known from the time-series analysis community.';
  s3='The data set will be used here to illustrate how';
  s4='pruning can enhance the network performance.';
  smat=str2mat(s0,s1,s2,s3,s4);

  % ---------- Load Data ----------
  load solplet.asc
  subplot('position',[0.1 0.55 0.45 0.38]);
  plot(solplet(:,1),solplet(:,2));
  xlabel('Year');
  ylabel('Sunspot activity');
  title('Sunspot benchmark data');
  set(gca,'Xlim',[min(solplet(:,1)) max(solplet(:,1))]);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 2  <<<<<<<<<<<<<<<<<
  s0='2';
  s1='Let''s begin by training a very large netork as';
  s2='the on-year-ahead predictor of the sunspot activity.';
  s3='The activity in a particular year is predicted';
  s4='from the activity in the 12 previous years.';
  smat=str2mat(s0,s1,s2,s3,s4);
  W1         = rand(8,13)-0.5;  % Weights to hidden layer 
  W2         = rand(1,9)-0.5;   % Weights to the output layer
  drawnet(W1,W2,eps,{'y(t-1) ' 'y(t-2) ' 'y(t-3) ' 'y(t-4) ' ...
                   'y(t-5) ' 'y(t-6) ' 'y(t-7) ' 'y(t-8) ' ...
                   'y(t-9) ' 'y(t-10)' 'y(t-11)' 'y(t-12)'},{'yhat(t)'});
  title('INITIAL FULLY CONNECTED NETWORK')
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 3  <<<<<<<<<<<<<<<<<
  s0='3';
  s1='Due to common practice, the data from 1700 to 1920';
  s2='is used for training while the data from 1921 to 1979';
  s3='is used for validation.';
  smat=str2mat(s0,s1,s2,s3);
  plot(solplet(:,1),solplet(:,2));
  xlabel('Year');
  ylabel('Sunspot activity');
  title('Sunspot benchmark data');
  set(gca,'Xlim',[min(solplet(:,1)) max(solplet(:,1))]);
  hold on
  plot([1921 1921],[0 1],'m--');
  hold off
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end 
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 4  <<<<<<<<<<<<<<<<<
  % ----- Train network -----
  s0='4';
  s1=[];
  s2='    >> Training process in action!! <<';
  s3=[];
  s4=[];
  s5='We run up to 800 iterations so you may have to';
  s6='wait for a while!';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6);
  set(edmulti,'String',smat);
  drawnow
  NetDef = ['HHHHHHHH';'L-------'];
  trparms = settrain;
  trparms=settrain(trparms,'maxiter',800,'D',[0.01 0.02]);
  NN = [12];
  y1=solplet(1:221,2)';
  y2=solplet(209:280,2)';
  [W1,W2,NSSEvec,iter]=nnarx(NetDef,NN,[],[],trparms,y1);
  delete(gca);
  subplot('position',[0.1 0.55 0.45 0.38]);
  semilogy(NSSEvec);
  set(gca,'Xlim',[0 iter]);
  xlabel('Iteration');
  title('Value of criterion');
  ylabel('Normalized SSE');
  grid
  

  % >>>>>>>>>>>>>>>>  SCREEN 5  <<<<<<<<<<<<<<<<<
  s0='5';
  s1='Now that the fully connected network has been';
  s2='trained to a minimum, we are ready to begin the';
  s3='pruning session.';
  s4=[];
  s5='We will eliminate the weights one by one, followed';
  s6='by a retraining.';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6);
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 6  <<<<<<<<<<<<<<<<<
  s0='6';
  s1='Pruning by optimal brain surgery is VERY time';
  s2='consuming and you may therefore choose to';
  s3='simply load the results.';
  smat=str2mat(s0,s1,s2,s3);
  set(edmulti,'String',smat);
  drawnow
  p=menu('Prune or Load?','Go prune','Just load the results','Quit');
  if p==3, close all, break; end
  if p==1,
   [thd,NSSEvec,FPEvec,NSSEtestvec,def,pv]=...
          nnprune('nnarx',NetDef,W1,W2,[],y1,NN,trparms,[50 0],[],y2);
   close(2)                        
  elseif p==2,
    load test7mat
    plot(pv,NSSEvec(pv),'x',pv,NSSEtestvec(pv),'o',pv,FPEvec(pv),'+')
    set(gca,'Xlim',[0 113]);
    xlabel('Parameters')
    title('x = training error,   + = FPE estimate,   o = test error')
    drawnow
  end
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 7  <<<<<<<<<<<<<<<<<
  s0='7';
  s1='The results are normalized and plotted again.';
  s2='The FPE estimate does not follow the test error';
  s3='very close, which has to do with the very abnormal';
  s4='progression of the sunspot activity towards the most';
  s5='recent years. Furthermore it looks as if not much is gained';
  s6='by pruning. The reason for this is, however, that the';
  s7='network has been trained using regularization.';
  smat=str2mat(s0,s1,s2,s3,s4,s5,s6,s7);
  minindex=27;
  vartotal=cov(solplet(:,2));
  plot(pv,2*FPEvec(pv)/vartotal,':',pv,2*NSSEvec(pv)/vartotal,'-',...
                                           pv,2*NSSEtestvec(pv)/vartotal,'--')
  hold on
  plot([minindex minindex],[0.07 0.3],'--m')
  hold off
  set(gca,'xlim',[0 120]);
  set(gca,'ylim',[0.07 0.3]);
  xlabel('Parameters')
  title('training (solid), test (dashed) and FPE (dotted)')
  pmnshow(smat,guihand,edmulti,contbut,quitbut);
  if StopDemo==1, close all, break; end
  
  
  % >>>>>>>>>>>>>>>>  SCREEN 8  <<<<<<<<<<<<<<<<<
  s0='8';
  s1='Let''s choose the architecture with 27 weigts as';
  s2='our final network.';
  [W1,W2]=netstruc(NetDef,thd,minindex);
  drawnet(W1,W2,eps,{'y(t-1) ' 'y(t-2) ' 'y(t-3) ' 'y(t-4) ' ...
                   'y(t-5) ' 'y(t-6) ' 'y(t-7) ' 'y(t-8) ' ...
                   'y(t-9) ' 'y(t-10)' 'y(t-11)' 'y(t-12)'},{'yhat(t)'});
  title('FINAL NETWORK ARCHITECTURE')
  smat=str2mat(s0,s1,s2);
  set(edmulti,'String',smat);
  drawnow
end
