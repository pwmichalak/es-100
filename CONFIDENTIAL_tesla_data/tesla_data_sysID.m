%% 0. clear and close everything; add OKID and ERA paths to this session
clear;
close all;
addpath('/Users/winstonmichalak/Desktop/harvard_academic_material/senior_year/thesis/es100-scripts/');

%% 1. obtain input current data u and output voltage data y 
[Channel_1_037, columns_037] = xlsread('MDCRvT_PANASONICNCRBB_101315_High_110315.xlsx','Channel_1-037');
C = num2cell(Channel_1_037);

% clean up the column names
columns_037 = strrep(columns_037, '(', '_');
columns_037 = strrep(columns_037, ')', '');
columns_037 = strrep(columns_037, ' ', '');
columns_037 = strrep(columns_037, '/', '_');
Channel_1_037 = cell2table(C, 'VariableNames', columns_037);
clear C columns_037;

%% 2. obtain data
Ts = 1; % sampling rate
x = 1:200; % data range
t = Channel_1_037.Test_Time_s(x); % test times
u = Channel_1_037.Current_A(x); % current (input) data
y = Channel_1_037.Voltage_V(x); % voltage (output) data

[u, tu] = resample(u, t, Ts, 1, 1);
[y, ty] = resample(y, t, Ts, 1, 1);

% make sure the time vectors are the same
assert(all(tu == ty));

%% 3. obtain state space matrix approximations using OKID / ERA
order = 2;
p = 40;
markovParams = okid(u,y,p);
[Ar, Br, Cr, Dr] = era(markovParams, order);

%% 4. obtain RC values from estimated state space matrices
sys_discrete = ss(Ar,Br,Cr,Dr,Ts);
%     sys = d2c(sys_discrete, 'foh');
sys = sys_discrete;
figure();
bode(sys);
title(strcat('OKID/ERA-Reconstructed System Bode Plot; p=', num2str(p)));
legend(strcat('Data Range: 1-', num2str(length(x))));
sysc = d2c(sys_discrete, 'foh');

% check that reconstucted is the same as original
yr = dlsim(Ar,Br,Cr,Dr,u);
figure()
hold on
plot(y)
plot(yr+y(1))

% extract RC parameters
[R0, R1, R2, Cap1, Cap2, T] = extract_2rc_params(Ar,Br,Cr,Dr);

%% 5. MOESP
% [ss_moesp,ssfun] = moesp(y,u,order);
% [Am,Bm,Cm,Dm] = ssfun(order); 
% sys_moesp_discrete = ss(Am,Bm,Cm,Dm,Ts);
% sys_moesp = d2c(sys_moesp_discrete);
% 
% figure()
% bode(sys_moesp);
% title('MOESP-Reconstructed System Bode Plot');