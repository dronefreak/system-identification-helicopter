function [prob]=FitnessFun(hd)
    prob=[];
    for i=1:5
       prob(i)= 1-(hd(i)/(hd(1)+hd(2)+hd(3)+hd(4)+hd(5)));
    end
    