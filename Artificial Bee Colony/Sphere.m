%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA114
% Project Title: Implementation of Artificial Bee Colony in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function z = Sphere(x,in,ou,t)
x=x';
Xu=x(1,1);
Xa=x(2,1);
Yv=x(3,1);
Yb=x(4,1);
Lu=x(5,1);
Lv=x(6,1);
Lb=x(7,1);
Lw=x(8,1);
Mu=x(9,1);
Mv=x(10,1);
Ma=x(11,1);
Mw=x(12,1);
Tf=x(13,1);
Ab=x(14,1);
Ac=x(15,1);
Ba=x(16,1);
Bd=x(17,1);
Za=x(18,1);
Zb=x(19,1);
Zw=x(20,1);
Zr=x(21,1);

Nv=x(22,1);
Np=x(23,1);
Nw=x(24,1);
Nr=x(25,1);
Nrfb=x(26,1);
Kr=x(27,1);
Krfb=x(28,1);
Ts=x(29,1);

Yped=x(30,1);
Mcol=x(31,1);
Alat=x(32,1);
Alon=x(33,1);
Blat=x(34,1);
Blon=x(35,1);
Zcol=x(36,1);
Nped=x(37,1);
Ncol=x(38,1);
Clon=x(39,1);
Dlat=x(40,1);

 A=[Xu,0,0,0,0,-9.81,Xa,0,0,0,0,0,0;
    0,Yv,0,0,9.81,0,0,Yb,0,0,0,0,0;
    Lu,Lv,0,0,0,0,0,Lb,Lw,0,0,0,0;
    Mu,Mv,0,0,0,0,Ma,0,Mw,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0,0,0,0;
    0,0,0,-Tf,0,0,-1,Ab,0,0,0,Ac,0;
    0,0,-Tf,0,0,0,Ba,-1,0,0,0,0,Bd;
    0,0,0,0,0,0,Za,Zb,Zw,Zr,0,0,0;
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
C=eye(13);
D=0;
mod=ss(A,B,C,D);


 y(:,:,1)=lsim(mod,in,t);
 cor1=sum(abs(y(:,1,1)-ou(:,1)));
 cor2=sum(abs(y(:,2,1)-ou(:,2)));
 cor3=sum(abs(y(:,3,1)-ou(:,3)));
 cor4=sum(abs(y(:,4,1)-ou(:,4)));
 cor5=sum(abs(y(:,5,1)-ou(:,5)));
 cor6=sum(abs(y(:,6,1)-ou(:,6)));
 cor7=sum(abs(y(:,7,1)-ou(:,7)));
 cor8=sum(abs(y(:,8,1)-ou(:,8)));
 cor9=sum(abs(y(:,9,1)-ou(:,9)));
 cor10=sum(abs(y(:,10,1)-ou(:,10)));
 
 d1=cor1;
 d2=cor2;
 d3=cor3;
 d4=cor4;
 d5=cor5;
 d6=cor6;
 d7=sqrt(sqrt(cor7));
 d8=sqrt(sqrt(cor8));
 d9=(cor9);
 d10=cor10;
 
 fit=abs((d1+d2+d3+d4+d5+d6+d7+d8+d9+d10));
    if isnan(fit)
        fit=1000000000;
    end
 cor1=corrcoef(y(:,1,1),ou(:,1));
 cor2=corrcoef(y(:,2,1),ou(:,2));
 cor3=corrcoef(y(:,3,1),ou(:,3));
 cor4=corrcoef(y(:,4,1),ou(:,4));
 cor5=corrcoef(y(:,5,1),ou(:,5));
 cor6=corrcoef(y(:,6,1),ou(:,6));
 cor7=corrcoef(y(:,7,1),ou(:,7));
 cor8=corrcoef(y(:,8,1),ou(:,8));
 cor9=corrcoef(y(:,9,1),ou(:,9));
 cor10=corrcoef(y(:,10,1),ou(:,10));
 
 d1=cor1(1,2);
 d2=cor2(1,2);
 d3=cor3(1,2);
 d4=cor4(1,2);
 d5=cor5(1,2);
 d6=cor6(1,2);
 d7=(cor7(1,2))^4;
 d8=(cor8(1,2))^4;
 d9=sqrt(sqrt(cor9(1,2)));
 d10=cor10(1,2);
 
cor=(d1+d2+d3+d4+d5+d6+d7+d8+d9+d10);
 if cor>1
 z=fit/cor^2;
 else
    z=fit/sqrt(cor); 
 end


end