function [distance]=Bhatta(A,B)
% varA=var(A);
%varB=var(B);
%distance=0.25*log(0.25*((varA/varB)+(varB/varA)+2))+0.25*((mean(A)-mean(B))^2)/varA+varB;
siz1=size(A);
min1=abs(min(min(A)));
min2=abs(min(min(B)));
if min1>min2
A=A+min1;
B=B+min1;
else
A=A+min2;
B=B+min2;
end
bc=0;
for i=1:siz1(1)
    bc=bc+sqrt(A(i)*B(i));
end
     distance=log(bc);
       if isnan(distance)
       distance= 1000;
       end
       if isinf(distance)
       distance= 1000;
       end

end


