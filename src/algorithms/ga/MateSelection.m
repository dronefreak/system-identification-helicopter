function [sel]=MateSelection(pop,prob)
   sel=[];
   sel=randsample(pop,2,true,prob);
  
   