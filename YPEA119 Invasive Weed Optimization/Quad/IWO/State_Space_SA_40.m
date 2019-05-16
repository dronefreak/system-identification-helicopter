function sys = State_Space_SA_40(BestPosition_Para, BestPosition_Mat)
g=9.8;
%------------Declaring Parameters(F Matrix)--------------------------------
x=BestPosition_Para';
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
%------------Declaring Parameters(G Matrix)--------------------------------
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
%------------Declaring Parameters(C Matrix)--------------------------------
z = BestPosition_Mat';
C1 = z(1,1);
C2 = z(2,1);
C3 = z(3,1);
C4 = z(4,1);
C5 = z(5,1);
C6 = z(6,1);
C7 = z(7,1);
C8 = z(8,1);
%------------Declaring Matrices--------------------------------------------
F=[Xu 0   0   0 0 -g  Xa     0     0   0  0     0      0
   0  Yv  0   0 g  0  0      Yb    0   0  0     0      0
   Lu Lv  0   0 0  0  0      Lb    Lw  0  0     0      0
   Mu Mv  0   0 0  0  Ma     0     Mw  0  0     0      0
   0  0   1   0 0  0  0      0     0   0  0     0      0
   0  0   0   1 0  0  0      0     0   0  0     0      0
   0  0   0  -1 0  0 -1/Tf   Ab/Tf 0   0  0     Ac/Tf  0
   0  0  -1   0 0  0  Ba/Tf -1/Tf  0   0  0     0      Bd/Tf 
   0  0   0   0 0  0  Za     Zb    Zw  Zr 0     0      0
   0  Nv  Np  0 0  0  0      0     Nw  Nr Nrfb  0      0
   0  0   0   0 0  0  0      0     0   Kr Krfb  0      0
   0  0   0  -1 0  0  0      0     0   0  0    -1/Ts   0
   0  0  -1   0 0  0  0      0     0   0  0     0     -1/Ts];
G=[0       0       0       0
   0       0       0       0
   0       0       Yped    0
   0       0       0       Mcol
   0       0       0       0
   0       0       0       0
   Alat/Tf Alon/Tf 0       0
   Blat/Tf Blon/Tf 0       0
   0       0       0       Zcol
   0       0       Nped    Ncol
   0       0       0       0
   0       Clon/Ts 0       0
   Dlat/Ts 0       0       0];
C=[C1 0  0  0  0  0  0  0  0  0  0  0  0   %u
   0  C2 0  0  0  0  0  0  0  0  0  0  0   %v
   0  0  C3 0  0  0  0  0  0  0  0  0  0   %p
   0  0  0  C4 0  0  0  0  0  0  0  0  0   %q
   0  0  0  0  C5 0  0  0  0  0  0  0  0   %phy
   0  0  0  0  0  C6 0  0  0  0  0  0  0   %theta
   0  0  0  0  0  0  0  0  0  0  0  0  0   %a
   0  0  0  0  0  0  0  0  0  0  0  0  0   %b
   0  0  0  0  0  0  0  0  C7 0  0  0  0   %w
   0  0  0  0  0  0  0  0  0  C8 0  0  0   %r
   0  0  0  0  0  0  0  0  0  0  0  0  0   %r_fb
   0  0  0  0  0  0  0  0  0  0  0  0  0   %c 
   0  0  0  0  0  0  0  0  0  0  0  0  0]; %d
D=[0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0
   0 0 0 0];
%---------------------State Space Representation---------------------------
sys = ss(F,G,C,D);
end