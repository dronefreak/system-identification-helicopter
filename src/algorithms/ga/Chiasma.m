function [modchrom1,modchrom2]=Chiasma(ChromA,ChromB)
    sizA=[];
    sizB=[];
    sizC=[];
    
    %Crossover for A
    A1=ChromA(1);
    A2=ChromB(1);
    
    sizA=size(A1);
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
      
    B1=ChromA(2);
    B2=ChromB(2);
    
    sizB=size(B1);
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
      C1=ChromA(3);
    C2=ChromB(3);
    
    sizC=size(C1);
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
    
    
    modchrom1=[A1,B1,C1,0];
    modchrom2=[A2,B2,C2,0];
    
 
    