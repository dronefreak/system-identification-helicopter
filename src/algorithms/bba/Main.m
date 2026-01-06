%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BBA source codes version 1.1                                     %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, X. Yang              %
%               Binary Bat Algorithm, Neural Computing and          %
%               Application, in press,                              %
%               DOI: 10.1007/s00521-013-1525-5                      %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all
clc

CostFunction=@(x) MyCost(x); % Modify or replace Mycost.m according to your cost funciton

Max_iteration=500; % Maximum number of iterations
noP=30; % Number of particles
noV=100;

A=.25;      % Loudness  (constant or decreasing)
r=.1;      % Pulse rate (constant or decreasing)

%BPSO with s-shaped family of transfer functions
[gBest, gBestScore ,ConvergenceCurve]=BBA(noP, A, r, noV, Max_iteration, CostFunction);

plot(ConvergenceCurve,'DisplayName','BBA','Color', 'r');
hold on


title(['\fontsize{12}\bf Convergence curve']);
xlabel('\fontsize{12}\bf Iteration');ylabel('\fontsize{12}\bf Average Best-so-far');
legend('\fontsize{10}\bf BBA',1);
grid on
axis tight

save resuls

