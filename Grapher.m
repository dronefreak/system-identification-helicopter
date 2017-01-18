
figure;
%  plot(t,outr(:,6),'black',t,iwo(:,6),'--b',t,ga(:,6),'--r','LineWidth',2);
% axis([0 Inf -0.5 0.5]);
%  title('Pitch Angle');
%  y1=smooth(t,outr(:,3),0.05,'loess');
% y2=smooth(t,iwo(:,3),0.05,'loess');
% y3=smooth(t,ga(:,3),0.05,'loess');
% 
%  plot(t,y1,'black',t,y2,'--b',t,y3,'--r','LineWidth',2);
plot(t,outr(:,3));

 axis([0 inf -4 4]);
 title('Roll Rate Response');


%  plot(t,inr(:,3),'LineWidth',2);
% 
% axis([0 inf -0.25 0.25]);
% title('Pedal Input');