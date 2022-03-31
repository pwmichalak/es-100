%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
% Third Iteration Project Build: Experiment 1                        %
%                                                                    %
% This experiment involves performing parameter extraction on        %
% the client’s data set using the state space identification method  %
% described as the build of this project report.                     %
%                                                                    %
% The third iteration of the project build will be used to perform   %
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
[r0p, r1p, r2p, c1p, c2p, sysc_build, mse] = build_iter3(u,y,p,order,Ts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. get client RC parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate = "rate2";
soc = 0.2;
temp = 20;

r1 = getRCparams(rate, "r1", soc, num2str(temp));
r2 = getRCparams(rate, "r2", soc, num2str(temp));
c1 = getRCparams(rate, "c1", soc, num2str(temp));
c2 = getRCparams(rate, "c2", soc, num2str(temp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. obtain parameter estimation errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1err = abs(max(r1p,r2p) - max(r1,r2)) / max(r1,r2);
r2err = abs(min(r1p,r2p) - min(r1,r2)) / min(r1,r2);
c1err = abs(max(c1p,c2p) - max(c1,c2)) / max(c1,c2);
c2err = abs(min(c1p,c2p) - min(c1,c2)) / min(c1,c2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. obtain client MSE for better comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. define state space matrices
Ac = [-1/(r1 * c1) 0; 0 -1/(r2 * c2)];
Bc = [1/c1; 1/c2];
Cc = [1 1];
Dc = 0;

sysc = ss(Ac,Bc,Cc,Dc); % define continuous time state space
sysd = c2d(sysc,Ts); % convert state space to discrete time
yr = dlsim(sysd.A,sysd.B,sysd.C,sysd.D,u); % simulate output data
ym = y - y(1);
mse_client = immse(ym,yr);