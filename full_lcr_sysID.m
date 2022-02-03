%% 0. clear and close everything
clear;
close all;

%% 1. define circuit elements
% RLC Circuit
L = 1e1 * rand(1,1); % define inductance
R = 1e3 * rand(1,1); % define resistance
Cap = 1e1 * rand(1,1); % define capacitance

% 2RC Circuit
R0 = 100; 
R1 = 1000; 
R2 = 100;
R3 = 150;
Cap1 = 0.1;
Cap2 = 0.15;
Cap3 = 0.2;

%% 2. define state space parameters

% The state variable x: [capacitor voltage ; inductor current]
% Input variable u: voltage before the inductor
% Output variable y: capacitor voltage
% 
% Ac: continuous state matrix A
% Bc: continuous input matrix B
% Cc: continuous output matrix C
% Dc: continuous feedthrough matrix D

% % THESE SYSTEM MATRICES ARE ALREADY DISCRETE
% Ac = [0.603 0.603 0 0;-0.603 0.603 0 0;0 0 -0.603 -0.603;0 0 0.603 -0.603];
% Bc = [1.1650,-0.6965;0.6268 1.6961;0.0751,0.0591;0.3516 1.7971];
% Cc = [0.2641,-1.4462,1.2460,0.5774;0.8717,-0.7012,-0.6390,-0.3600];
% Dc = [-0.1356,-1.2704;-1.3493,0.9846];
% 
% % THESE SYSTEM MATRICES ARE CONTINUOUS
% Ac = [-1/(R * Cap) 1/Cap; -1/L 0];
% Bc = [0; 1/L];
% Cc = [1 0];
% Dc = 0;
%
% % THESE SYSTEM MATRICES ARE CONTINUOUS
Ac = [-1/(R1 * Cap1) 0 0; 0 -1/(R2 * Cap2) 0; 0 0 -1/(R3 * Cap3)];
Bc = [1/Cap1; 1/Cap2; 1/Cap3];
Cc = [1 1 1];
Dc = R0;

%% 3. obtain random input data voltage u and output data voltage y
q = size(Dc,1); % output size
m = size(Dc,2); % input size
order = size(Ac,2); % order of the model
n = 1000; % number of samples
Ts = 1; rand(1,1); % sampling time

% simulate data
u = randn(n,m);  
sys_orig = ss(Ac,Bc,Cc,Dc); % define continuous time state space
sysd = c2d(sys_orig,Ts); % convert state space to discrete time
y = dlsim(sysd.A,sysd.B,sysd.C,sysd.D,u); % simulate output data
% y = dlsim(Ac,Bc,Cc,Dc,u); % IF THE SYSTEM IS ALREADY DISCRETE
% sys_orig = ss(Ac,Bc,Cc,Dc,Ts); % IF THE SYSTEM IS ALREADY DISCRETE

%% 4. obtain pulse response using OKID 
p = 7;
markovParams = okid(u,y,p); % these markov parameters are in discrete time. 
[bruntonMarkovParams, ~] = bruntonOKID(y',u',order);

%% 5. analyze pulse response using ERA
[Ar, Br, Cr, Dr] = era(markovParams, order);
mco = floor((length(bruntonMarkovParams)-1)/2);
[Arb, Brb, Crb, Drb, ~] = bruntonERA(bruntonMarkovParams, mco, mco, m, q, order);

%% 6. MOESP
% [ss_moesp,ssfun] = moesp(y,u,order);
% [Am,Bm,Cm,Dm] = ssfun(order); 
% sys_moesp = ss(Am,Bm,Cm,Dm,Ts);
% 
% % plot
% figure();
% bode(sys_moesp,'b', sys_orig , 'r');
% legend('MOESP','original');

%% 7. N4SID 
% data = iddata(y,u,Ts);
% sys_n4sid = n4sid(data,order, 'Ts',Ts, 'Feedthrough',1, 'DisturbanceModel','None');
% 
% % plot
% figure();
% bode(sys_n4sid,'b', sys_orig , 'r');
% legend('N4SID','original');

%% 8. obtain OKID / ERA reconstructed state space parameters and compare

% convert from discrete to continuous time reconstructed OKID state space
sys_okid_discrete = ss(Ar,Br,Cr,Dr,Ts);
sys_okid = d2c(sys_okid_discrete);
% sys_okid = ss(Ar,Br,Cr,Dr,Ts);

sys_brunton_discrete = ss(Arb, Brb, Crb, Drb,Ts);
sys_brunton = d2c(sys_brunton_discrete);
% sys_brunton = ss(Arb, Brb, Crb, Drb,Ts);

% [R0r, R1r, R2r, Cap1r, Cap2r, Tr] = extract_2rc_params(sys_okid.A,sys_okid.B,sys_okid.C,sys_okid.D);
% [R0r, R1r, R2r, R3r, Cap1r, Cap2r, Cap3r, Tr] = extract_3rc_params(sys_okid.A,sys_okid.B,sys_okid.C,sys_okid.D);

% compare reconstructed state spaces
figure()
bode(sys_okid, 'b', sys_orig , 'r');
legend('OKID / ERA','Original');

figure()
bode(sys_orig , 'r', sys_brunton, 'black');
legend('Original','Brunton');

% figure();
% bode(sys_moesp,'b', sys_orig , 'r', sys_okid, 'black', sys_n4sid, 'green');
% legend('MOESP','original', 'OKID/ERA', 'N4SID');