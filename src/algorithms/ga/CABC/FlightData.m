inr(:,1)=RC(:,2);                          %Lateral input
inr(:,2)=RC(:,3);                          %Longitudinal input
inr(:,3)=RC(:,5);                          %Pedal input
inr(:,4)=RC(:,7);                          %Collective


outr(:,1)=LPOS(:,8);
outr(:,2)=LPOS(:,7);

atmax=size(ATT);
atmin=size(LPOS);
outr=resample(outr,atmax(1),atmin(1));
outr(:,6)=ATT(:,9);    
outr(:,5)=ATT(:,10);
outr(:,4)=ATT(:,6);
outr(:,3)=ATT(:,7);
outr(:,10)=ATT(:,11);
atmax=size(IMU1);
atmin=size(ATT);
outr=resample(outr,atmax(1),atmin(1));
outr(:,9)=IMU1(:,4);

max=size(outr);
min=size(inr);

inr=resample(inr,max(1),min(1));

outr(:,7)=-6.9*inr(:,2);
outr(:,8)=-6.9*inr(:,1);

% com=ATTC;
% pil=inr(:,3);
% maxfb=size(com);
% minfb=size(pil);
% diff=[];
% pil=resample(pil,maxfb,minfb);
% diff=com(:,4)-pil(:,3);
% 
% maxfb=size(outr);
% minfb=size(diff);
% 
% if maxfb>minfb
%     diff=resample(diff,maxfb,minfb);
% else
%     outr=resample(outr,minfb,maxfb);
% end
% 
% outr(:,11)=diff;

 clearvars -except inr outr population ATTC