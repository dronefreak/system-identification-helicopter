psize = 5;
pop=1:psize;
t=0:0.01:136.21;
maxIter=10000;
for i=1:maxIter

% P=[Xu,Xa,Yv,Yb,Lu,Lv,Lb,Lw,Mu,Mv,Ma,Mw,Tf,Ab,Ac,Ba,Bd,Za,Zb,Zw,Zr,Nv,Np,Nw,Nr,Nrfb,Kr,Krfb,Ts,Yped,Mcol,Alat,Alon,Blat,Blon,Zcol,Nped,Ncol,Clon,Dlat];
for j=1:psize
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


 y(:,:,j)=lsim(mod,inr,t);
 cor1=corrcoef(y(:,1,j),outr(:,1));
 cor2=corrcoef(y(:,2,j),outr(:,2));
 cor3=corrcoef(y(:,3,j),outr(:,3));
 cor4=corrcoef(y(:,4,j),outr(:,4));
 cor5=corrcoef(y(:,5,j),outr(:,5));
 cor6=corrcoef(y(:,6,j),outr(:,6));
 cor7=corrcoef(y(:,9,j),outr(:,9));
 cor8=corrcoef(y(:,10,j),outr(:,10));
 cor9=corrcoef(y(:,7,j),outr(:,7));
 cor10=corrcoef(y(:,8,j),outr(:,8));
 cor11=corrcoef(y(:,11,j),outr(:,11));
 d1=cor1(1,2);
 d2=cor2(1,2);
 d3=cor3(1,2);
 d4=cor4(1,2);
 d5=cor5(1,2);
 d6=cor6(1,2);
 d7=cor7(1,2);
 d8=cor8(1,2);
 d9=cor9(1,2);
 d10=cor10(1,2);
 d11=cor11(1,2);
 
 fitness(j)=(d1+d2+d3+d4+d5+d6+d7+d8+d9+d10+d11)*9.09;
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
   
   sizpop=size(population);
   chiasma=[randsample(sizpop(1),1),randsample(sizA(2),1)];  
   
   if rem(i,5)==0
   mutpoint=[randsample(sizpop(1),1),randsample(sizA(2),1)]
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
Yv=parA(2);
Lu=parA(3);
Lv=parA(4);
Lb=parA(5);
La=parA(6);
Lq=parA(7);
Mu=parA(8);
Mv=parA(9);
Ma=parA(10);
Mb=parA(11);
Mq=parA(12);
Zw=parA(13);
Zr=parA(14);
Nw=parA(15);
Nr=parA(16);
Nrfb=parA(17);
Kr=parA(18);
Krfb=parA(19);
Ts=parA(20);
k1=parA(21);


Alat=parA(22);
Blon=parA(23);
Zt=parA(24);
Zm=parA(25);
Nm=parA(26);
Nt=parA(27);
Xlat=parA(28);
Ylon=parA(29);
Mlat=parA(30);
Mlon=parA(31);
Llat=parA(32);
Llon=parA(33);

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


 y(:,:,1)=lsim(mod,inr,t);
  cor1=corrcoef(y(:,1,1),outr(:,1));
 cor2=corrcoef(y(:,2,1),outr(:,2));
 cor3=corrcoef(y(:,3,1),outr(:,3));
 cor4=corrcoef(y(:,4,1),outr(:,4));
 cor5=corrcoef(y(:,5,1),outr(:,5));
 cor6=corrcoef(y(:,6,1),outr(:,6));
 cor7=corrcoef(y(:,9,1),outr(:,9));
 cor8=corrcoef(y(:,10,1),outr(:,10));
 cor9=corrcoef(y(:,7,1),outr(:,7));
 cor10=corrcoef(y(:,8,1),outr(:,8));
 cor11=corrcoef(y(:,11,1),outr(:,11));
 d1=cor1(1,2);
 d2=cor2(1,2);
 d3=cor3(1,2);
 d4=cor4(1,2);
 d5=cor5(1,2);
 d6=cor6(1,2);
 d7=cor7(1,2);
 d8=cor8(1,2);
 d9=cor9(1,2);
 d10=cor10(1,2);
 d11=cor11(1,2);
 
 fit1=(d1+d2+d3+d4+d5+d6+d7+d8+d9+d10+d11)*9.09;
 
 [minfit,n]=min(fitness);
 
 if minfit<fit1
     population(:,n)=parA;
     dude=mod;
 end
   
 parA=parB;
 
Xu=parA(1);
Yv=parA(2);
Lu=parA(3);
Lv=parA(4);
Lb=parA(5);
La=parA(6);
Lq=parA(7);
Mu=parA(8);
Mv=parA(9);
Ma=parA(10);
Mb=parA(11);
Mq=parA(12);
Zw=parA(13);
Zr=parA(14);
Nw=parA(15);
Nr=parA(16);
Nrfb=parA(17);
Kr=parA(18);
Krfb=parA(19);
Ts=parA(20);
k1=parA(21);


Alat=parA(22);
Blon=parA(23);
Zt=parA(24);
Zm=parA(25);
Nm=parA(26);
Nt=parA(27);
Xlat=parA(28);
Ylon=parA(29);
Mlat=parA(30);
Mlon=parA(31);
Llat=parA(32);
Llon=parA(33);


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


 y(:,:,1)=lsim(mod,inr,t);
 cor1=corrcoef(y(:,1,1),outr(:,1));
 cor2=corrcoef(y(:,2,1),outr(:,2));
 cor3=corrcoef(y(:,3,1),outr(:,3));
 cor4=corrcoef(y(:,4,1),outr(:,4));
 cor5=corrcoef(y(:,5,1),outr(:,5));
 cor6=corrcoef(y(:,6,1),outr(:,6));
 cor7=corrcoef(y(:,9,1),outr(:,9));
 cor8=corrcoef(y(:,10,1),outr(:,10));
 cor9=corrcoef(y(:,7,1),outr(:,7));
 cor10=corrcoef(y(:,8,1),outr(:,8));
 cor11=corrcoef(y(:,11,1),outr(:,11));
 d1=cor1(1,2);
 d2=cor2(1,2);
 d3=cor3(1,2);
 d4=cor4(1,2);
 d5=cor5(1,2);
 d6=cor6(1,2);
 d7=cor7(1,2);
 d8=cor8(1,2);
 d9=cor9(1,2);
 d10=cor10(1,2);
 d11=cor11(1,2);
 
 fit1=(d1+d2+d3+d4+d5+d6+d7+d8+d9+d10+d11)*9.09;
 [minfit,n]=min(fitness);
 
 if minfit<fit1
     population(:,n)=parA;
     dude=mod;
 end
 
 display(i,'Generation');
 display(fitness,'Fitness');
 
end


[A,B,x,D]=ssdata(dude);
csize=5;
citer=0
for o=0:citer 
  for j=1:csize
      
C=[cpop(1,j),0,0,0,0,0,0,0,0,0,0,0,0;
    0,cpop(2,j),0,0,0,0,0,0,0,0,0,0,0;
    0,0,cpop(3,j),0,0,0,0,0,0,0,0,0,0;
    0,0,0,cpop(4,j),0,0,0,0,0,0,0,0,0;
    0,0,0,0,cpop(5,j),0,0,0,0,0,0,0,0;
    0,0,0,0,0,cpop(6,j),0,0,0,0,0,0,0;
    0,0,0,0,0,0,cpop(7,j),0,0,0,0,0,0;
    0,0,0,0,0,0,0,cpop(8,j),0,0,0,0,0;
    0,0,0,0,0,0,0,0,cpop(9,j),0,0,0,0;
    0,0,0,0,0,0,0,0,0,cpop(10,j),0,0,0;
    0,0,0,0,0,0,0,0,0,0,cpop(11,j),0,0;
    0,0,0,0,0,0,0,0,0,0,0,cpop(12,j),0;
    0,0,0,0,0,0,0,0,0,0,0,0,cpop(13,j)];

mod=ss(A,B,C,D);


 y(:,:,j)=lsim(mod,inr,t);

 
 bd1=(norm(y(:,1,j)-outr(:,1)));
 bd2=(norm(y(:,2,j)-outr(:,2)));
 bd3=(norm(y(:,3,j)-outr(:,3)));
 bd4=(norm(y(:,4,j)-outr(:,4)));
 bd5=(norm(y(:,5,j)-outr(:,5)));
 bd6=(norm(y(:,6,j)-outr(:,6)));
 bd7=(norm(y(:,9,j)-outr(:,9)));
 bd8=(norm(y(:,10,j)-outr(:,10)));
 
d=(bd1+bd2+bd3+bd4+bd5+bd6+bd7+bd8);
 
 
 
 fitness(j)=1/d;
  end

   for k=1:csize
       prob(k)=(fitness(k)/sum(fitness));
       if prob(k)<0
           prob(k)=0;
       end
   end
   
   sel=randsample(pop,2,true,prob);
   
   parA=cpop(:,sel(1));
   parB=cpop(:,sel(2));
   
   minA=abs(min(parA));
   minB=abs(min(parB));
   
   parA=floor((parA+minA)*1000);
   parB=floor((parB+minB)*1000);
   
   binparA=decimalToBinaryVector(parA,32);
   binparB=decimalToBinaryVector(parB,32);
   
   sizA=size(binparA);
   sizB=size(binparB);
   
   
   chiasma=[randsample(13,1),randsample(sizA(2),1)];  %For 13 parameters
   
   if rem(i,5)==0
   mutpoint=[randsample(13,1),randsample(sizA(2),1)]
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
  C=[parA(1),0,0,0,0,0,0,0,0,0,0,0,0;
    0,parA(2),0,0,0,0,0,0,0,0,0,0,0;
    0,0,parA(3),0,0,0,0,0,0,0,0,0,0;
    0,0,0,parA(4),0,0,0,0,0,0,0,0,0;
    0,0,0,0,parA(5),0,0,0,0,0,0,0,0;
    0,0,0,0,0,parA(6),0,0,0,0,0,0,0;
    0,0,0,0,0,0,parA(7),0,0,0,0,0,0;
    0,0,0,0,0,0,0,parA(8),0,0,0,0,0;
    0,0,0,0,0,0,0,0,parA(9),0,0,0,0;
    0,0,0,0,0,0,0,0,0,parA(10),0,0,0;
    0,0,0,0,0,0,0,0,0,0,parA(11),0,0;
    0,0,0,0,0,0,0,0,0,0,0,parA(12),0;
    0,0,0,0,0,0,0,0,0,0,0,0,parA(13)];
  
    mod=ss(A,B,C,D);


 y(:,:,1)=lsim(mod,inr,t);
  
 
 bd1=(norm(y(:,1,1)-outr(:,1)));
 bd2=(norm(y(:,2,1)-outr(:,2)));
 bd3=(norm(y(:,3,1)-outr(:,3)));
 bd4=(norm(y(:,4,1)-outr(:,4)));
 bd5=(norm(y(:,5,1)-outr(:,5)));
 bd6=(norm(y(:,6,1)-outr(:,6)));
 bd7=(norm(y(:,9,1)-outr(:,9)));
 bd8=(norm(y(:,10,1)-outr(:,10)));
 
 
d=(bd1+bd2+bd3+bd4+bd5+bd6+bd7+bd8);
 
 
 fit1=1/d;
 
 [minfit,n]=min(fitness);
 
 if minfit<fit1
     cpop(:,n)=parA;
 end
   
 parA=parB;
 
  C=[parA(1),0,0,0,0,0,0,0,0,0,0,0,0;
    0,parA(2),0,0,0,0,0,0,0,0,0,0,0;
    0,0,parA(3),0,0,0,0,0,0,0,0,0,0;
    0,0,0,parA(4),0,0,0,0,0,0,0,0,0;
    0,0,0,0,parA(5),0,0,0,0,0,0,0,0;
    0,0,0,0,0,parA(6),0,0,0,0,0,0,0;
    0,0,0,0,0,0,parA(7),0,0,0,0,0,0;
    0,0,0,0,0,0,0,parA(8),0,0,0,0,0;
    0,0,0,0,0,0,0,0,parA(9),0,0,0,0;
    0,0,0,0,0,0,0,0,0,parA(10),0,0,0;
    0,0,0,0,0,0,0,0,0,0,parA(11),0,0;
    0,0,0,0,0,0,0,0,0,0,0,parA(12),0;
    0,0,0,0,0,0,0,0,0,0,0,0,parA(13)];
  
    mod=ss(A,B,C,D);


 y(:,:,1)=lsim(mod,inr,t);
  
 
 bd1=(norm(y(:,1,1)-outr(:,1)));
 bd2=(norm(y(:,2,1)-outr(:,2)));
 bd3=(norm(y(:,3,1)-outr(:,3)));
 bd4=(norm(y(:,4,1)-outr(:,4)));
 bd5=(norm(y(:,5,1)-outr(:,5)));
 bd6=(norm(y(:,6,1)-outr(:,6)));
 bd7=(norm(y(:,9,1)-outr(:,9)));
 bd8=(norm(y(:,10,1)-outr(:,10)));
 
d=(bd1+bd2+bd3+bd4+bd5+bd6+bd7+bd8);
 
 
 fit1=1/d;
 
 if minfit<fit1
     cpop(:,n)=parA;
 end
 
 display(o,'Generation');
 display(fitness,'Fitness');
end
