function [ChromA,ChromB]=Chrominator(sel)
    
    [A1,B1,C1,D1] = ssdata(population(sel(1)));
    [A2,B2,C2,D2] = ssdata(population(sel(2)));
    %Makeup for A
    A1= A1*1000+500;
    A2=A2*1000+500;
    
    A1=de2bi(A1);
    A2=de2bi(A2);
    
    %Makeup for B
    
    B1= B1*1000+500;
    B2=B2*1000+500;
    
    B1=de2bi(B1);
    B2=de2bi(B2);
    
    %Makeup for C
    C1= C1*1000+500;
    C2=C2*1000+500;
    
    C1=de2bi(C1);
    C2=de2bi(C2);
    
    D1=0;
    D2=0;
    
    ChromA=[A1,B1,C1,D1];
    ChromB=[A2,B2,C2,D2];
    