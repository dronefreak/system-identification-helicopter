
psize =20;
pop=1:psize;
t=0:0.01:85.89;
maxIter=10000;

for i=1:maxIter

% P=[Xu,Xa,Yv,Yb,Lu,Lv,Lb,Lw,Mu,Mv,Ma,Mw,Tf,Ab,Ac,Ba,Bd,Za,Zb,Zw,Zr,Nv,Np,Nw,Nr,Nrfb,Kr,Krfb,Ts,Yped,Mcol,Alat,Alon,Blat,Blon,Zcol,Nped,Ncol,Clon,Dlat];
for j=1:psize
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


 y(:,:,j)=lsim(mod,inr,t);
 cor1=sum(abs(y(:,1,j)-ou(:,1)));
 cor2=sum(abs(y(:,2,j)-ou(:,2)));
 cor3=sum(abs(y(:,3,j)-ou(:,3)));
 cor4=sum(abs(y(:,4,j)-ou(:,4)));
 cor5=sum(abs(y(:,5,j)-ou(:,5)));
 cor6=sum(abs(y(:,6,j)-ou(:,6)));
 cor7=sum(abs(y(:,7,j)-ou(:,7)));
 cor8=sum(abs(y(:,8,j)-ou(:,8)));
 cor9=sum(abs(y(:,9,j)-ou(:,9)));
 cor10=sum(abs(y(:,10,j)-ou(:,10)));
 
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
 cor1=corrcoef(y(:,1,j),ou(:,1));
 cor2=corrcoef(y(:,2,j),ou(:,2));
 cor3=corrcoef(y(:,3,j),ou(:,3));
 cor4=corrcoef(y(:,4,j),ou(:,4));
 cor5=corrcoef(y(:,5,j),ou(:,5));
 cor6=corrcoef(y(:,6,j),ou(:,6));
 cor7=corrcoef(y(:,7,j),ou(:,7));
 cor8=corrcoef(y(:,8,j),ou(:,8));
 cor9=corrcoef(y(:,9,j),ou(:,9));
 cor10=corrcoef(y(:,10,j),ou(:,10));
 
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
 fitness(j)=fit/cor^2;
 else
    fitness(j)=fit/sqrt(cor); 
 end
 
end

   for k=1:psize
       prob(k)=(fitness(k)/sum(fitness));
       if prob(k)<0
           prob(k)=0;
       end
   end
   
   sel=randsample(pop,2,true,prob);
   
   parA=population(:,sel(1));
   parB=population(:,sel(2));
   
   minA=abs(min(parA));
   minB=abs(min(parB));
   
   parA=floor((parA+minA)*1000);
   parB=floor((parB+minB)*1000);
   
   binparA=decimalToBinaryVector(parA);
   binparB=decimalToBinaryVector(parB);
   
   sizA=size(binparA);
   sizB=size(binparB);
   
   
   chiasma=[randsample(40,1),randsample(sizA(2),1)];  %For 40 parameters
   
   if rem(i,5)==0
   mutpoint=[randsample(40,1),randsample(sizA(2),1)]
   if binparA(mutpoint(1),mutpoint(2))==0
       binparA(mutpoint(1),mutpoint(2))=1;
   else
       binparA(mutpoint(1),mutpoint(2))=0;
   end
   end
   for l=1:chiasma(1)
        for m=1:chiasma(2)
            tempA(l,m)=binparA(l,m);
        end
    end
    
    for l=1:chiasma(1)
        for m=1:chiasma(2)
            binparA(l,m)=binparB(l,m);
        end
    end
    
    for l=1:chiasma(1)
        for m=1:chiasma(2)
            binparB(l,m)=tempA(l,m);
        end
    end
    
    parA=binaryVectorToDecimal(binparA);
    parB=binaryVectorToDecimal(binparB);
    
    parA= (parA/1000)-minA;
    parB= (parB/1000)-minB;
    
  %ParA fitness check
  Xu=parA(1);
  Xa=parA(2);
  Yv=parA(3);
Yb=parA(4);
Lu=parA(5);
Lv=parA(6);
Lb=parA(7);
Lw=parA(8);
Mu=parA(9);
Mv=parA(10);
Ma=parA(11);
Mw=parA(12);
Tf=parA(13);
Ab=parA(14);
Ac=parA(15);
Ba=parA(16);
Bd=parA(17);
Za=parA(18);
Zb=parA(19);
Zw=parA(20);
Zr=parA(21);
Zrfb=0;
Nv=parA(22);
Np=parA(23);
Nw=parA(24);
Nr=parA(25);
Nrfb=parA(26);
Kr=parA(27);
Krfb=parA(28);
Ts=parA(29);

Yped=parA(30);
Mcol=parA(31);
Alat=parA(32);
Alon=parA(33);
Blat=parA(34);
Blon=parA(35);
Zcol=parA(36);
Nped=parA(37);
Ncol=parA(38);
Clon=parA(39);
Dlat=parA(40);

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
C=eye(13);
D=0;
mod=ss(A,B,C,D);


 y(:,:,1)=lsim(mod,inr,t);
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
 fit1=fit/cor^2;
 else
    fit1=fit/sqrt(cor); 
 end
 
 [minfit,n]=min(fitness);
 
 if minfit<fit1
     population(:,n)=parA;
     dude=mod;
 end
   
 parA=parB;
 
  Xu=parA(1);
  Xa=parA(2);
  Yv=parA(3);
Yb=parA(4);
Lu=parA(5);
Lv=parA(6);
Lb=parA(7);
Lw=parA(8);
Mu=parA(9);
Mv=parA(10);
Ma=parA(11);
Mw=parA(12);
Tf=parA(13);
Ab=parA(14);
Ac=parA(15);
Ba=parA(16);
Bd=parA(17);
Za=parA(18);
Zb=parA(19);
Zw=parA(20);
Zr=parA(21);
Zrfb=0;
Nv=parA(22);
Np=parA(23);
Nw=parA(24);
Nr=parA(25);
Nrfb=parA(26);
Kr=parA(27);
Krfb=parA(28);
Ts=parA(29);

Yped=parA(30);
Mcol=parA(31);
Alat=parA(32);
Alon=parA(33);
Blat=parA(34);
Blon=parA(35);
Zcol=parA(36);
Nped=parA(37);
Ncol=parA(38);
Clon=parA(39);
Dlat=parA(40);

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
C=eye(13);
D=0;
mod=ss(A,B,C,D);


 y(:,:,1)=lsim(mod,inr,t);
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
 fit1=fit/cor^2;
 else
    fit1=fit/sqrt(cor); 
 end
 
 [minfit,n]=min(fitness);
 
 if minfit<fit1
     population(:,n)=parA;
     dude=mod;
 end
 
 display(i,'Generation');
 display(fitness,'Fitness');
 his(i)=minfit;
end
semilogy(his,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

population=population';