j=1;
Xu=popul(1,j);
Xa=popul(2,j);
Yv=popul(3,j);
Yb=popul(4,j);
Lu=popul(5,j);
Lv=popul(6,j);
Lb=popul(7,j);
Lw=popul(8,j);
Mu=popul(9,j);
Mv=popul(10,j);
Ma=popul(11,j);
Mw=popul(12,j);
Tf=popul(13,j);
Ab=popul(14,j);
Ac=popul(15,j);
Ba=popul(16,j);
Bd=popul(17,j);
Za=popul(18,j);
Zb=popul(19,j);
Zw=popul(20,j);
Zr=popul(21,j);
Zrfb=0;
Nv=popul(22,j);
Np=popul(23,j);
Nw=popul(24,j);
Nr=popul(25,j);
Nrfb=popul(26,j);
Kr=popul(27,j);
Krfb=popul(28,j);
Ts=popul(29,j);

Yped=popul(30,j);
Mcol=popul(31,j);
Alat=popul(32,j);
Alon=popul(33,j);
Blat=popul(34,j);
Blon=popul(35,j);
Zcol=popul(36,j);
Nped=popul(37,j);
Ncol=popul(38,j);
Clon=popul(39,j);
Dlat=popul(40,j);


 A=[Xu,0,0,0,0,-9.81,Xa,0,0,0,0,0,0;
    0,Yv,0,0,9.81,0,0,Yb,0,0,0,0,0;
    Lu,Lv,0,0,0,0,0,Lb,Lw,0,0,0,0;
    Mu,Mv,0,0,0,0,Ma,0,Mw,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0,0,0,0;
    0,0,0,-Tf,0,0,-1,Ab,0,0,0,Ac,0;
    0,0,-Tf,0,0,0,Ba,-1,0,0,0,0,Bd;
    0,0,0,0,0,0,Za,Zb,Zw,Zr,Zrfb,0,0;
    0,Nv,Np,0,0,0,0,0,Nw,Nr,Nrfb,0,0;
    0,0,0,0,0,0,0,0,0,Kr,Krfb,0,0;
    0,0,0,-Ts,0,0,0,0,0,0,0,-1,0;
    0,0,-Ts,0,0,0,0,0,0,0,0,0,-1];

 B=[0,0,0,0;
    0,0,Yped,0;
    0,0,0,0;
    0,0,0,Mcol;
    0,0,0,0;
    0,0,0,0;
    Alat,Alon,0,0;
    Blat,Blon,0,0;
    0,0,0,Zcol;
    0,0,Nped,Ncol;
    0,0,0,0;
    0,Clon,0,0;
    Dlat,0,0,0];
C=[cpop(1,1),0,0,0,0,0,0,0,0,0,0,0,0;
    0,cpop(2,1),0,0,0,0,0,0,0,0,0,0,0;
    0,0,cpop(3,1),0,0,0,0,0,0,0,0,0,0; 
    0,0,0,cpop(4,1),0,0,0,0,0,0,0,0,0;
    0,0,0,0,cpop(5,1),0,0,0,0,0,0,0,0;
    0,0,0,0,0,cpop(6,1),0,0,0,0,0,0,0;
    0,0,0,0,0,0,cpop(7,1),0,0,0,0,0,0;
    0,0,0,0,0,0,0,cpop(8,1),0,0,0,0,0;
    0,0,0,0,0,0,0,0,cpop(9,1),0,0,0,0;
    0,0,0,0,0,0,0,0,0,cpop(10,1),0,0,0;
    0,0,0,0,0,0,0,0,0,0,cpop(11,1),0,0;
    0,0,0,0,0,0,0,0,0,0,0,cpop(12,1),0;
    0,0,0,0,0,0,0,0,0,0,0,0,cpop(13,1)];


D=0;
mod=ss(A,B,C,D);

y=lsim(mod,inr,t);
figure
subplot(4,3,7)
plot(t,y(:,6),'-.',t,outr(:,6))
axis([0 inf -1 1])
title('Pitch angle')

subplot(4,3,8)
plot(t,y(:,5),'-.',t,outr(:,5))
axis([0 inf -2 2])
title('Roll angle')

subplot(4,3,5)
plot(t,y(:,4),'-.',t,outr(:,4))
title('Pitch rate')

subplot(4,3,4)
plot(t,y(:,3),'-.',t,outr(:,3))
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
axis([0 inf -10 10])
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

