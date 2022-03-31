%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
% First Iteration Project Build: Experiment 1                        %
%                                                                    %
% This experiment involves performing parameter extraction on        %
% the client’s data set using the state space identification method  %
% described as the build of this project report.                     %
%                                                                    %
% The first iteration of the project build will be used to perform   %
% this experiment.                                                   % 
%                                                                    %
% The data for this experiment is confidential and is therefore not  %
% included in the associated experiment directory.                   %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. obtain client data file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = 'MDCRvT_PANASONICNCRBB_101315_High_110315.xlsx';
channel = 'Channel_1-037';
[~,~,~,channel] = get_battery_data(file, channel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. obtain client current input data u and output voltage data y 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 1; % sampling frequency
Ts = 1/fs; % sampling period
x = 1:800; % the data range x defines the length of the data vector that is pertinent for parameter estimation
t = channel.Test_Time_s(x); % test times
u = channel.Current_A(x); % current (input) data
y = channel.Voltage_V(x); % voltage (output) data

% resample data to ensure that pulses are evenly spaced in time
[u, tu] = resample(u, t, Ts, 1, 1);
[y, ty] = resample(y, t, Ts, 1, 1);

% make sure the time vectors are the same
assert(all(tu == ty));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. obtain discrete state space matrix approximations using project build
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 40;
order = 2;
fs = 1; Ts = 1/fs;
[r0p, r1p, r2p, c1p, c2p, sysc_build, mse] = build_iter1(u,y,p,order,Ts);



