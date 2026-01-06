j=1;
Xu=population(1,j);
Yv=population(2,j);
Lu=population(3,j);
Lv=population(4,j);
Lb=population(5,j);
La=population(6,j);
Lq=population(7,j);
Mu=population(8,j);
Mv=population(9,j);
Ma=population(10,j);
Mb=population(11,j);
Mq=population(12,j);
Zw=population(13,j);
Zr=population(14,j);
Nw=population(15,j);
Nr=population(16,j);
Nrfb=population(17,j);
Kr=population(18,j);
Krfb=population(19,j);
Ts=population(20,j);
k1=population(21,j);


Alat=population(22,j);
Blon=population(23,j);
Zt=population(24,j);
Zm=population(25,j);
Nm=population(26,j);
Nt=population(27,j);
Xlat=population(28,j);
Ylon=population(29,j);
Mlat=population(30,j);
Mlon=population(31,j);
Llat=population(32,j);
Llon=population(33,j);

 A=[Xu 0 -9.81 0 0 0 -9.81 0 0 0 0
     0 Yv 0 9.81 0 0 0 9.81 0 0 0
     0 0 0 0 1 0 0 0 0 0 0 
     0 0 0 0 0 1 0 0 0 0 0
     Mu Mv 0 0 Mq 0 Ma Mb 0 0 0
     Lu Lv 0 0 0 Lq La Lb 0 0 0
     0 0 0 0 -1 k1 -1/Ts 0 0 0 0
     0 0 0 0 -k1 -1 0 -1/Ts 0 0 0
     0 0 0 0 0 0 0 0 Zw Zr 0
     0 0 0 0 0 0 0 0 Nw Nr Nrfb
     0 0 0 0 0 0 0 0 0 Kr -Krfb];

 B=[Xlat,0,0,0;
    0,Ylon,0,0;
    0,0,0,0;
    0,0,0,0;
    Mlat,Mlon,0,0;
    Llat,Llon,0,0;
    Alat,0,0,0;
    0,Blon,0,0;
    0,0,Zm,Zt;
    0,0,Nm,Nt;
    0,0,0,0];
C=eye(11);
D=0;
mod=ss(A,B,C,D);


y=lsim(mod,inr,t);
figure
subplot(4,3,7)
plot(t,y(:,3),'-.',t,outr(:,3))
axis([0 inf -1 1])
title('Pitch angle')

subplot(4,3,8)
plot(t,y(:,4),'-.',t,outr(:,4))
axis([0 inf -2 2])
title('Roll angle')

subplot(4,3,5)
plot(t,y(:,5),'-.',t,outr(:,5))
title('Pitch rate')

subplot(4,3,4)
plot(t,y(:,6),'-.',t,outr(:,6))
axis([0 inf -2 2])
title('Roll Rate')


subplot(4,3,6)
plot(t,y(:,10),'-.',t,outr(:,10))
axis([0 inf -5 5])
title('Yaw Rate')

subplot(4,3,10)
plot(t,y(:,1),'-.',t,outr(:,1))
axis([0 inf -20 20])
title('u')

subplot(4,3,11)
plot(t,y(:,2),'-.',t,outr(:,2))
axis([0 inf -5 5])
title('v')

subplot(4,3,2)
plot(t,inr(:,1))
title('Roll input')

subplot(4,3,1)
plot(t,inr(:,2))
title('Pitch input')

subplot(4,3,3)
plot(t,inr(:,3))
title('Yaw input')


subplot(4,3,12)
plot(t,y(:,9),'-.',t,outr(:,9))
axis([0 inf -10 10])
title('w')

subplot(4,3,9)
plot(t,inr(:,4))
title('Collective input')