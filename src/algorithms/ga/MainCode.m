t=0:0.1:10;
u=zeros(size(t));
hd=[];
population=[ss1,ss2,ss3,ss4,ss5];
for k=1:100
    
hd=Initialize(population(1),population(2),population(3),population(4),population(5),u,t,io);
prob=FitnessFun(hd);
sel=MateSelection(population,prob);
[ChromA,ChromB]=Chrominator(sel);
[modchrom1,modchrom2]=Chiasma(ChromA,ChromB);
[modsys1,modsys2]=Dechrominator(modchrom1,modchrom2);
y1=lsim(modsys1,u,t);
y2=lsim(modsys2,u,t);
hd1=ModHausdorffDist(y1,io);
hd2=ModHausdorffDist(y2,io);

[~,i]=max(hd);
hd(i)=hd1;
population(i)=modsys1;
[~,j]=max(hd);
hd(j)=hd2;
population(j)=modsys2;

[m,l]=min(hd);
ymin=lsim(population(l),u,t);
plot(ymin,io,t);
title('Generation',k);
end

