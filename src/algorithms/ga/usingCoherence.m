%Using COherence

psize = 5;
pop=1:psize;
t=0:0.01:8.18;
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

outsize=size(outr);

 y(:,:,j)=lsim(mod,inr,t);
 cor1=mscohere(y(:,1,j),outr(:,1));
 cor2=mscohere(y(:,2,j),outr(:,2));
 cor3=mscohere(y(:,3,j),outr(:,3));
 cor4=mscohere(y(:,4,j),outr(:,4));
 cor5=mscohere(y(:,5,j),outr(:,5));
 cor6=mscohere(y(:,6,j),outr(:,6));
 cor7=mscohere(y(:,9,j),outr(:,9));
 cor8=mscohere(y(:,10,j),outr(:,10));
 d1=sum(cor1)/outsize(1);
 d2=sum(cor2)/outsize(1);
 d3=sum(cor3)/outsize(1);
 d4=sum(cor4)/outsize(1);
 d5=sum(cor5)/outsize(1);
 d6=sum(cor6)/outsize(1);
 d7=sum(cor7)/outsize(1);
 d8=sum(cor8)/outsize(1);
 
 
 fitness(j)=(d1+d2+d3+d4+d5+d6+d7+d8)*12.5;
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


 y(:,:,1)=lsim(mod,inr,t);
   cor1=mscohere(y(:,1,1),outr(:,1));
 cor2=mscohere(y(:,2,1),outr(:,2));
 cor3=mscohere(y(:,3,1),outr(:,3));
 cor4=mscohere(y(:,4,1),outr(:,4));
 cor5=mscohere(y(:,5,1),outr(:,5));
 cor6=mscohere(y(:,6,1),outr(:,6));
 cor7=mscohere(y(:,9,1),outr(:,9));
 cor8=mscohere(y(:,10,1),outr(:,10));
 d1=sum(cor1)/outsize(1);
 d2=sum(cor2)/outsize(1);
 d3=sum(cor3)/outsize(1);
 d4=sum(cor4)/outsize(1);
 d5=sum(cor5)/outsize(1);
 d6=sum(cor6)/outsize(1);
 d7=sum(cor7)/outsize(1);
 d8=sum(cor8)/outsize(1);
 
 
 
 
 
 fit1=((d1+d2+d3+d4+d5+d6+d7+d8)*12.5);
 
 [minfit,n]=min(fitness);
 
 if minfit<fit1
     population(:,n)=parA;
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


 y(:,:,1)=lsim(mod,inr,t);
  cor1=mscohere(y(:,1,1),outr(:,1));
 cor2=mscohere(y(:,2,1),outr(:,2));
 cor3=mscohere(y(:,3,1),outr(:,3));
 cor4=mscohere(y(:,4,1),outr(:,4));
 cor5=mscohere(y(:,5,1),outr(:,5));
 cor6=mscohere(y(:,6,1),outr(:,6));
 cor7=mscohere(y(:,9,1),outr(:,9));
 cor8=mscohere(y(:,10,1),outr(:,10));
 d1=sum(cor1)/outsize(1);
 d2=sum(cor2)/outsize(1);
 d3=sum(cor3)/outsize(1);
 d4=sum(cor4)/outsize(1);
 d5=sum(cor5)/outsize(1);
 d6=sum(cor6)/outsize(1);
 d7=sum(cor7)/outsize(1);
 d8=sum(cor8)/outsize(1);
 fit1=((d1+d2+d3+d4+d5+d6+d7+d8)*12.5);
 [minfit,n]=min(fitness);
 
 if minfit<fit1
     population(:,n)=parA;
 end
 
 display(i,'Generation');
 display(fitness,'Fitness');
 
end