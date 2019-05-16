% clear all;
% load('freq1.mat','-mat');
% clearvars inr outr;
% load('2.mat','-mat');
% %load('BestSol_Para','-mat');
% Flightdatacreator;
% % bbo;
% % population(1,:)=BestSol.Position;
% % disp('Starting BBO');
for xn=1:10
iwo;
BC(xn,:)=BestSol.Cost;
trial(xn,:)=BestSol.Position;
end