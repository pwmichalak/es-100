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
x = 1:800; % data range
t = Channel_1_037.Test_Time_s(x); % test times
u = Channel_1_037.Current_A(x); % current (input) data
y = Channel_1_037.Voltage_V(x); % voltage (output) data

[u, tu] = resample(u, t, Ts, 1, 1);
[y, ty] = resample(y, t, Ts, 1, 1);

% make sure the time vectors are the same
assert(all(tu == ty));

%% 3. obtain state space matrix approximations using OKID / ERA
% pvals = 10:10:50;
% i = 1;
% RC2_params = zeros(5,length(pvals));
% for p = pvals
p = 40;
order = 2;
markovParams = okid(u,y,p);
[Ar, Br, Cr, Dr] = era(markovParams, order);

%% 4. obtain RC values from estimated state space matrices
sys_discrete = ss(Ar,Br,Cr,Dr,Ts);
% sys = d2c(sys_discrete, 'foh');
sys = sys_discrete;
figure();
bode(sys);
title(strcat('OKID/ERA-Reconstructed System Bode Plot; p=', num2str(p)));
legend(strcat('Data Range: 1-', num2str(length(x))));

% check that reconstucted is the same as original
figure()
hold on
plot(y)
yr = dlsim(Ar,Br,Cr,Dr,u);
plot(yr+y(1),'DisplayName',strcat('p=',num2str(p)))
sys_discrete = ss(Ar,Br,Cr,Dr,Ts);

%     sysc = d2c(sys_discrete, 'foh');
%     sysc = sys_discrete;
%     Ac = sysc.A;
%     Bc = sysc.B;
%     Cc = sysc.C;
%     Dc = sysc.D;
%     [R0, R1, R2, Cap1, Cap2, T] = extract_2rc_params(Ac,Bc,Cc,Dc);
%     RC2_params(:,i) = [R0 R1 R2 Cap1 Cap2]';
%     i = i+1;
% end

xlabel('Time(s)');
ylabel('Voltage(V)');
title('Reconstruction of time domain current-voltage data');
legend();
% legend('Original data',strcat('Reconstructed data, o=',num2str(order)));

% figure()
% hold on
% plot(pvals, RC2_params(1:3,:))
% xlabel('p value');
% ylabel('Resistance value');
% legend('R0','R1','R2');
% 
% figure()
% hold on
% plot(pvals, RC2_params(4:5,:))
% xlabel('p value');
% ylabel('Capacitance value');
% legend('Cap1','Cap2');

% extract RC parameters 
% NEED TO CONVERT BACK TO CONTINUOUS TIME STATE SPACE FIRST
% Ac = sysc.A;
% Bc = sysc.B;
% Cc = sysc.C;
% Dc = sysc.D;
% [R0, R1, R2, R3, Cap1, Cap2, Cap3, T] = extract_3rc_params(Ac,Bc,Cc,Dc);

%% 5. MOESP
% [ss_moesp,ssfun] = moesp(y,u,order);
% [Am,Bm,Cm,Dm] = ssfun(order); 
% sys_moesp_discrete = ss(Am,Bm,Cm,Dm,Ts);
% sys_moesp = d2c(sys_moesp_discrete);
% 
% figure()
% bode(sys_moesp);
% title('MOESP-Reconstructed System Bode Plot');