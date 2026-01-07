function [modsys1,modsys2]=Dechrominator(modchrom1,modchrom2)
A1=bi2de(modchrom1(1));
B1=bi2de(modchrom1(2));
C1=bi2de(modchrom1(3));
D1=0;

A2=bi2de(modchrom2(1));
B2=bi2de(modchrom2(2));
C2=bi2de(modchrom2(3));
D2=0;

A1=transpose(reshape(A1,2,2));
A2=transpose(reshape(A2,2,2));



C1=transpose(C1);
C2=transpose(C2);

modsys1=ss(A1,B1,C1,D1);
modsys2=ss(A2,B2,C2,D2);
