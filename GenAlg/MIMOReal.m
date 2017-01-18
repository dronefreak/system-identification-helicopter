%% Real Coded Genetic Algorithm for System Identification of Small Scaled Unmanned Helicopter
%Navaneethkrishnan B
%%Initialisation
population={ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9,ss10,ss11,ss12,ss13,ss14,ss15,ss16,ss17,ss18,ss19,ss20};                                           %Initial Population
psize=20;
pop=1:psize;                                                               %Population index
maxIter=50000;                                                                 %Maximum number of iteration
t=0:6696;
t=transpose(t);
%Loop
for ll=1:maxIter
   
   
    for i=1:psize
    y(:,:,i)=lsim(population{1,i},inr,t);
    end
    for j=1:psize
      corr1(:,:,j)=corrcoef(y(:,1,j),outr(:,1));
      corr2(:,:,j)=corrcoef(y(:,2,j),outr(:,2));
      corr3(:,:,j)=corrcoef(y(:,3,j),outr(:,3));
      corr4(:,:,j)=corrcoef(y(:,4,j),outr(:,4));
      corr5(:,:,j)=corrcoef(y(:,5,j),outr(:,5));
      corr6(:,:,j)=corrcoef(y(:,6,j),outr(:,6));
      
    end
    
    for j=1:psize
      d1(j)=corr1(1,2,j);
      d2(j)=corr2(1,2,j);
      d3(j)=corr3(1,2,j);
      d4(j)=corr4(1,2,j);
      d5(j)=corr5(1,2,j);
      d6(j)=corr6(1,2,j);
      hd(j)= d1(j)+d2(j)+d3(j)+d4(j)+d5(j)+d6(j);
    end
    for j=1:psize
        if(hd(j)<0)
            hd(j)=0;
        end
        if(isnan(hd(j)))
            hd(j)=0;
        end
    end
    %% Genetic Algorithm
    for k=1:psize
       prob(k)=(hd(k)/sum(hd));
    end
    sel=[];
    sel=randsample(pop,2,true,prob);
 
     A1=[];
     B1=[];
     C1=[];
     
    [A1,B1,C1,D1] = ssdata(population{1,sel(1)});
    
    [A2,B2,C2,D2] = ssdata(population{1,sel(2)});
    %display(size(A1),'Size of A');
    minA1=abs(min(min(A1)));
    minA2=abs(min(min(A2)));
    minB1=abs(min(min(B1)));
    minB2=abs(min(min(B2)));
    minC1=abs(min(min(C1)));
    minC2=abs(min(min(C2)));
    minD1=abs(min(min(D1)));
    minD2=abs(min(min(D2)));
    %Makeup for A
    A1= floor((A1+minA1)*1000);
    A2=floor((A2+minA2)*1000);
    %Makeup for B
    
    B1= floor((B1+minB1)*1000);
    B2=floor((B2+minB2)*1000);
    
    
    %Makeup for C
    C1= floor((C1+minC1)*1000);
    C2=floor((C2+minC2)*1000);
    
 
    %{
    D1= floor((D1+minD1)*1000);
    D2=floor((D2+minD2)*1000);
    
    %ChromA=[A1,B1,C1,D1];
    %ChromB=[A2,B2,C2,D2];
    %}
    sizA=[];
    sizB=[];
    sizC=[];
    
    
    D=zeros(6,4);
    
    %Crossover for A
   % A1=ChromA(1);
    %A2=ChromB(1);
     sizA=size(A1);
     maxA=max(max(A1));
     %%Mutation for A every 5 generations
    %if rem(l,5)==0
        %Mutation
        if rem(ll,5)==0
           
               mutpoint1=[randsample(sizA(1),1),randsample(sizA(2),1)];
            mutpoint2=[randsample(sizA(1),1),randsample(sizA(2),1)];
              temp= A1(mutpoint1(1),mutpoint1(2));
              A1(mutpoint1(1),mutpoint1(2))= A1(mutpoint2(1),mutpoint2(2));
              A1(mutpoint2(1),mutpoint2(2))=temp;
           
           % end
        end
    
   
    chisma=[randsample(sizA(1),1),randsample(sizA(2),1)];
    tempA=[];
    for i=1:chisma(1)
        for j=1:chisma(2)
            tempA(i,j)=A1(i,j);
        end
    end
    
    for i=1:chisma(1)
        for j=1:chisma(2)
            A1(i,j)=A2(i,j);
        end
    end
    
    for i=1:chisma(1)
        for j=1:chisma(2)
            A2(i,j)=tempA(i,j);
        end
    end
    % Crossover for B
      
   % B1=ChromA(2);
    %B2=ChromB(2);
    
    sizB=size(B1);
    maxB=max(max(B1));
    %Mutation
    if rem(ll,5)==2
        mutpoint1=[randsample(sizB(1),1),randsample(sizB(2),1)];
            mutpoint2=[randsample(sizB(1),1),randsample(sizB(2),1)];
              temp= B1(mutpoint1(1),mutpoint1(2));
              B1(mutpoint1(1),mutpoint1(2))= B1(mutpoint2(1),mutpoint2(2));
              B1(mutpoint2(1),mutpoint2(2))=temp;
    end
    
    
    chisma=[randsample(sizB(1),1),randsample(sizB(2),1)];
    tempB=[];
    for i=1:chisma(1)
        for j=1:chisma(2)
            tempB(i,j)=B1(i,j);
        end
    end
    
    for i=1:chisma(1)
        for j=1:chisma(2)
            B1(i,j)=B2(i,j);
        end
    end
    
    for i=1:chisma(1)
        for j=1:chisma(2)
            B2(i,j)=tempB(i,j);
        end
    end
    
    %Crossover of C
    %  C1=ChromA(3);
    %C2=ChromB(3);
    
    sizC=size(C1);
    maxC=max(max(C1));
    %Mutation
    if rem(ll,5)==4
      mutpoint1=[randsample(sizC(1),1),randsample(sizC(2),1)];
            mutpoint2=[randsample(sizC(1),1),randsample(sizC(2),1)];
              temp= C1(mutpoint1(1),mutpoint1(2));
              C1(mutpoint1(1),mutpoint1(2))= C1(mutpoint2(1),mutpoint2(2));
              C1(mutpoint2(1),mutpoint2(2))=temp;
    end
    
    chisma=[randsample(sizC(1),1),randsample(sizC(2),1)];
    tempC=[];
    for i=1:chisma(1)
        for j=1:chisma(2)
            tempC(i,j)=C1(i,j);
        end
    end
    
    for i=1:chisma(1)
        for j=1:chisma(2)
            C1(i,j)=C2(i,j);
        end
    end
    
    for i=1:chisma(1)
        for j=1:chisma(2)
            C2(i,j)=tempC(i,j);
        end
    end
    
    
   %{ 
     sizD=size(D1);
     maxD=max(max(D1));
     %%Mutation for D every 5 generations
    
        if rem(ll,10)==0
            mutpoint1=[randsample(sizD(1),1),randsample(sizD(2),1)];
            mutpoint2=[randsample(sizD(1),1),randsample(sizD(2),1)];
              temp= D1(mutpoint1(1),mutpoint1(2));
              D1(mutpoint1(1),mutpoint1(2))= D1(mutpoint2(1),mutpoint2(2));
             D1(mutpoint2(1),mutpoint2(2))=temp;
           
        end
    
   
    chisma=[randsample(sizD(1),1),randsample(sizD(2),1)];
    tempD=[];
    for i=1:chisma(1)
        for j=1:chisma(2)
            tempD(i,j)=D1(i,j);
        end
    end
    
    for i=1:chisma(1)
        for j=1:chisma(2)
            D1(i,j)=D2(i,j);
        end
    end
    
    for i=1:chisma(1)
        for j=1:chisma(2)
            D2(i,j)=tempD(i,j);
        end
    end
    %}
   
A1=(A1/1000)-minA1;
A2=(A2/1000)-minA2;
B1=(B1/1000)-minB1;
B2=(B2/1000)-minB2;
C1=(C1/1000)-minC1;
C2=(C2/1000)-minC2;
D1=0;
D2=0;

modsys1=ss(A1,B1,C1,D1);
modsys2=ss(A2,B2,C2,D2);


    y(:,:,1)=lsim(modsys1,inr,t);
   
      corr1(:,:,1)=corrcoef(y(:,1,1),outr(:,1));
      corr2(:,:,1)=corrcoef(y(:,2,1),outr(:,2));
      corr3(:,:,1)=corrcoef(y(:,3,1),outr(:,3));
      corr4(:,:,1)=corrcoef(y(:,4,1),outr(:,4));
      corr5(:,:,1)=corrcoef(y(:,5,1),outr(:,5));
      corr6(:,:,1)=corrcoef(y(:,6,1),outr(:,6));
      
      d1=corr1(1,2,1);
      d2=corr2(1,2,1);
      d3=corr3(1,2,1);
      d4=corr4(1,2,1);
      d5=corr5(1,2,1);
      d6=corr6(1,2,1);
      hd1= d1+d2+d3+d4+d5+d6;
    
    y(:,:,2)=lsim(modsys2,inr,t);
      corr1(:,:,2)=corrcoef(y(:,1,1),outr(:,1));
      corr2(:,:,2)=corrcoef(y(:,2,1),outr(:,2));
      corr3(:,:,2)=corrcoef(y(:,3,1),outr(:,3));
      corr4(:,:,2)=corrcoef(y(:,4,1),outr(:,4));
      corr5(:,:,2)=corrcoef(y(:,5,1),outr(:,5));
      corr6(:,:,2)=corrcoef(y(:,6,1),outr(:,6));
      
      d1=corr1(1,2,2);
      d2=corr2(1,2,2);
      d3=corr3(1,2,2);
      d4=corr4(1,2,2);
      d5=corr5(1,2,2);
      d6=corr6(1,2,2);
      hd2= d1+d2+d3+d4+d5+d6;
    
    
    [~,i]=min(hd);
if hd(i)<hd1
hd(i)=hd1;
population{1,i}=modsys1;
end
[q,j]=min(hd);
if hd(j)<hd2;
hd(j)=hd2;
population{1,j}=modsys2;
end
q=max(hd);
clc;
display(ll,'Generation');
display(q*20,'Fitness');
worstdistance(ll)=q;



end
close all;
semilogy(worstdistance,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
