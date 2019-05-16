

function z = Sphere(x,in,ou,t)
%x=BestSol_Para;
x=x';
%------------Declaring Parameters(F Matrix)--------------------------------
Lu=x(1,1);
Lp=x(2,1);
Lv=x(3,1);
Lw=x(4,1);
Lq=x(5,1);
Lr=x(6,1);

%------------Declaring Matrices--------------------------------------------
A=[-Lu 0 0 0 1 0  0 Lp     
   0 -Lv 0 1 0 0 Lp  0    
   0 0 -Lw 0 0 0 0 0    
   0 0 0 -Lp 0 0 0 0     
   0 0 0 0 -Lq 0 0 0     
   0 0 0 0 0 -Lr 0 0     
   0 0 0 1 0 0 0 0 
   0 0 0 0 1 0 0 0  ];

B=[0       0       0       0
   0       0       0       0
   Lw      0       0       0
   0       Lp      0       0
   0       0       Lq      0
   0       0       0       Lr
   0       0       0       0
   0       0       0       0
   ];

C=eye(8);
D=0;
mod=ss(A,B,C,D);
 y(:,:,1)=lsim(mod,in,t);
%  out_data=fft(ou);
% test_out=fft(y);
%  z=0;
%  for i=1:10
%  z=z+ goodnessOfFit(test_out(:,i),out_data(:,i),'MSE');
%      %z=z + sum(abs(out_data(:,i)-test_out(:,i)));
%  end
%  
%  if isreal(z)
%  else
%      z=inf;
%  end
%  
     
 
 
 
 
 
 
%  
%  cor1=sum(abs(y(:,1,1)-ou(:,1)));
%  cor2=sum(abs(y(:,2,1)-ou(:,2)));
%  cor3=sum(abs(y(:,3,1)-ou(:,3)));
%  cor4=sum(abs(y(:,4,1)-ou(:,4)));
%  cor5=sum(abs(y(:,5,1)-ou(:,5)));
%  cor6=sum(abs(y(:,6,1)-ou(:,6)));
%  cor7=sum(abs(y(:,7,1)-ou(:,7)));
%  cor8=sum(abs(y(:,8,1)-ou(:,8)));
%  cor9=sum(abs(y(:,9,1)-ou(:,9)));
%  cor10=sum(abs(y(:,10,1)-ou(:,10)));
%  
%  d1=cor1;
%  d2=cor2;
%  d3=cor3;
%  d4=cor4;
%  d5=cor5;
%  d6=cor6;
%  d7=sqrt(sqrt(cor7));
%  d8=sqrt(sqrt(cor8));
%  d9=(cor9);
%  d10=cor10;
%  
%  fit=abs((d1+d2+d3+d4+d5+d6+d7+d8+d9+d10));
%     if isnan(fit)
%         fit=1000000000;
%     end
%     
 cor1=corrcoef(y(:,1,1),ou(:,1));
 cor2=corrcoef(y(:,2,1),ou(:,2));
 cor3=corrcoef(y(:,3,1),ou(:,3));
 cor4=corrcoef(y(:,4,1),ou(:,4));
 cor5=corrcoef(y(:,5,1),ou(:,5));
 cor6=corrcoef(y(:,6,1),ou(:,6));
 cor7=corrcoef(y(:,7,1),ou(:,7));
 cor8=corrcoef(y(:,8,1),ou(:,8));
 
 d1=abs(cor1(1,2));
 d2=abs(cor2(1,2));
 d3=abs(cor3(1,2));
 d4=abs(cor4(1,2));
 d5=abs(cor5(1,2));
 d6=abs(cor6(1,2));
 d7=abs(cor7(1,2));
 d8=abs(cor8(1,2));
 
cor= d1+d2+d3+d4+d5+d6+d7+d8;
if isreal(cor)
  z=8-cor;
end
    if isnan(cor)
        z=Inf;
    end

 
  


end