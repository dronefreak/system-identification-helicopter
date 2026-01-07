j=1;
Xu=population(1,j);
Xa=population(2,j);
Yv=population(3,j);
Yb=population(4,j);
Lu=population(5,j);
Lv=population(6,j);
Lb=population(7,j);
Lw=population(8,j);
Mu=population(9,j);
Mv=population(10,j);
Ma=population(11,j);
Mw=population(12,j);
Tf=population(13,j);
Ab=population(14,j);
Ac=population(15,j);
Ba=population(16,j);
Bd=population(17,j);
Za=population(18,j);
Zb=population(19,j);
Zw=population(20,j);
Zr=population(21,j);
Zrfb=0;
Nv=population(22,j);
Np=population(23,j);
Nw=population(24,j);
Nr=population(25,j);
Nrfb=population(26,j);
Kr=population(27,j);
Krfb=population(28,j);
Ts=population(29,j);

Yped=population(30,j);
Mcol=population(31,j);
Alat=population(32,j);
Alon=population(33,j);
Blat=population(34,j);
Blon=population(35,j);
Zcol=population(36,j);
Nped=population(37,j);
Ncol=population(38,j);
Clon=population(39,j);
Dlat=population(40,j);


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

