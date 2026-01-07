 function [hd,pop] = Initialize(ss1,ss2,ss3,ss4,ss5,u,t,io)
population=[ss1,ss2,ss3,ss4,ss5];
pop=[1,2,3,4,5];
y=[];
  for i=1:5
    y(:,i)=lsim(population(i),u,t);
  end
 hd=[];
 
 for j=1:5
     hd(j)=ModHausdorffDist(y(:,j),io);
 end
 
 

