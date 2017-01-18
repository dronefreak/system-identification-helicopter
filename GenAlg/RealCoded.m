   %%Consolidated
%population={ss1,ss2,ss3,ss4,ss5};
pop=[1,2,3,4,5];
for ll=1:100
y=[];
  for i=1:5
    y(:,i)=lsim(population{1,i},u,t);
  end
 hd=[];
 
 for j=1:5
     hd(j)=ModHausdorffDist(y(:,j),io);
 end

  prob=[];
    for k=1:5
       prob(k)= 1-(hd(k)/(hd(1)+hd(2)+hd(3)+hd(4)+hd(5)));
    end
    sel=[];
   sel=randsample(pop,2,true,prob);
 
     A1=[];
     B1=[];
     C1=[];
     
    [A1,B1,C1,D1] = ssdata(population{1,sel(1)});
    
    [A2,B2,C2,D2] = ssdata(population{1,sel(2)});
    display(size(A1),'Size of A');
    minA1=abs(min(min(A1)));
    minA2=abs(min(min(A2)));
    minB1=abs(min(min(B1)));
    minB2=abs(min(min(B2)));
    minC1=abs(min(min(C1)));
    minC2=abs(min(min(C2)));
    %Makeup for A
    A1= floor((A1+minA1)*1000);
    A2=floor((A2+minA2)*1000);
    %Makeup for B
    
    B1= floor((B1+minB1)*1000);
    B2=floor((B2+minB2)*1000);
    
    
    %Makeup for C
    C1= floor((C1+minC1)*1000);
    C2=floor((C2+minC2)*1000);
    
 
    D1=0;
    D2=0;     
    
    %ChromA=[A1,B1,C1,D1];
    %ChromB=[A2,B2,C2,D2];
    
    sizA=[];
    sizB=[];
    sizC=[];
    
    %Crossover for A
   % A1=ChromA(1);
    %A2=ChromB(1);
     sizA=size(A1);
     maxA=max(max(A1));
     %%Mutation for A every 5 generations
    %if rem(l,5)==0
        %Mutation
        if rem(ll,3)==0
            mutpoint=[randsample(sizA(1),1),randsample(sizA(2),1)];
               A1(mutpoint(1),mutpoint(2))=randsample(maxA,1);
           
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
    if rem(ll,3)==0
            mutpoint=[randsample(sizB(1),1),randsample(sizB(2),1)];
               B1(mutpoint(1),mutpoint(2))=randsample(maxB,1);
           
    end
    
    
    chisma=[randsample(sizB(1),1),randsample(sizB(2),1)];
    tempB=[];
    for i=1:chisma(1)
        for j=1:chisma(2)
            tempB(i,j)=A1(i,j);
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
    if rem(ll,3)==0
            mutpoint=[randsample(sizC(1),1),randsample(sizC(2),1)];
               C1(mutpoint(1),mutpoint(2))=randsample(maxC,1);
           
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
    
   
A1=(A1/1000)-minA1;
A2=(A2/1000)-minA2;
B1=(B1/1000)-minB1;
B2=(B2/1000)-minB2;
C1=(C1/1000)-minC1;
C2=(C2/1000)-minC2;

modsys1=ss(A1,B1,C1,D1);
modsys2=ss(A2,B2,C2,D2);

y1=lsim(modsys1,u,t);
y2=lsim(modsys2,u,t);
hd1=ModHausdorffDist(y1,io);
hd2=ModHausdorffDist(y2,io);

[~,i]=max(hd);
if hd(i)>hd1
hd(i)=hd1;
population{1,i}=modsys1;
end
[~,j]=max(hd);
if hd(j)>hd2;
hd(j)=hd2;
population{1,j}=modsys2;
end
[m,l]=min(hd);
ymin=lsim(population{1,l},u,t);
display(ll,'Generation');
display(m,'Minimum Distance')
display(hd,'Distances')
%display(size(A1),'Size of A')
end
