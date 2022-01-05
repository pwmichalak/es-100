%% 0. clear and close everything
clear;
close all;

%% 1. define circuit elements

% define inductance
L = 1e1 * rand(1,1); 

% define resistance
R = 1e3 * rand(1,1);

% define capacitance
Cap = 1e1 * rand(1,1);

%% 2. define state space parameters

% The state variable x: [capacitor voltage ; inductor current]
% Input variable u: voltage before the inductor
% Output variable y: capacitor voltage
% 
% Ac: continuous state matrix A
% Bc: continuous input matrix B
% Cc: continuous output matrix C
% Dc: continuous feedthrough matrix D
Ac = [-1/(R * Cap) 1/Cap 1 1; -1/L 0 1 1; 1 1 1 1; 1 1 1 1];
Bc = [0; 1/L; 1; 1];
Cc = [1 0 0 0; 0 1 0 0];
Dc = [0;0];

Ac = [0.603 0.603 0 0;-0.603 0.603 0 0;0 0 -0.603 -0.603;0 0 0.603 -0.603];
Bc = [1.1650,-0.6965;0.6268 1.6961;0.0751,0.0591;0.3516 1.7971];
Cc = [0.2641,-1.4462,1.2460,0.5774;0.8717,-0.7012,-0.6390,-0.3600];
Dc = [-0.1356,-1.2704;-1.3493,0.9846];

Ac = [-1/(R * Cap) 1/Cap; -1/L 0];
Bc = [0 1; 1/L 1];
Cc = [1 0; 0 1];
Dc = [0 1;0 1];
% 
% Ac = [-1/(R * Cap) 1/Cap 1; -1/L 0 1; 1 1 1];
% Bc = [0; 1/L; 1];
% Cc = [1 0 0];
% Dc = [0];
% 
Ac = [1 1 1; 1 1 0; 1 0 1];
Bc = [1; 1; 1];
Cc = [1 0 0];
Dc = 0;

Ac = [-1/(R * Cap) 1/Cap; -1/L 0];
Bc = [0; 1/L];
Cc = [1 0; 0 1];
Dc = [0;0];

Ac = [-1/(R * Cap) 1/Cap; -1/L 0];
Bc = [0; 1/L];
Cc = [1 0];
Dc = 0;

% Ac = [1 1; 1 0];
% Bc = [0; 1];
% Cc = [1 0];
% Dc = 0;

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

% for debugging
% u = [0.1403 -1.8421 0.5369 -1.1585 -1.1526 0.0467 0.0213 -0.6087 0.4874 0.9303]';
% y = [0 0.0112 -0.1054 -0.4168 -0.9845 -2.1745 -4.3752 -8.3128 -15.4987 -28.5991]';

%% 4. obtain pulse response using OKID 
p = 7;
markovParams = okid(u,y,p); % these markov parameters are in discrete time. 
tmp_markovParams = tmp_okid(u,y,p);
[bruntonMarkovParams, ~] = bruntonOKID(y',u',order);

% markovParams variable for sanity checking later on
trueMP = {[Dc] [Cc * Bc]};
for i=1:p-1
    trueMP{i+2} = Cc * Ac ^ i * Bc;
end

%% 5. analyze pulse response using ERA
[Ar, Br, Cr, Dr] = era(markovParams, order);
[Art, Brt, Crt, Drt] = tmp_era(tmp_markovParams, order);
mco = floor((length(bruntonMarkovParams)-1)/2);
[Arb, Brb, Crb, Drb, ~] = bruntonERA(bruntonMarkovParams, mco, mco, m, q, order);

%% 6. obtain OKID / ERA reconstructed state space parameters and compare

% % n4sid
% data = iddata(y,u,Ts);
% sys_n4sid = n4sid(data,order,'DisturbanceModel','None');

% moesp
[ss_moesp,ssfun] = moesp(y,u,order);
[Am,Bm,Cm,Dm] = ssfun(order); 
sys_moesp = ss(Am,Bm,Cm,Dm);

% covert from discrete to continuous time reconstructed OKID state space
sys_okid_discrete = ss((Ar),(Br),(Cr),Dr,Ts);
sys_okid = d2c(sys_okid_discrete);

tmp_sys_okid_discrete = ss((Art),(Brt),(Crt),Drt,Ts);
tmp_sys_okid = d2c(tmp_sys_okid_discrete);

sys_brunton_discrete = ss(Arb, Brb, Crb, Drb,Ts);
sys_brunton = d2c(sys_brunton_discrete);
tfb = tf(sys_brunton);

% compare OKID and reconstructed state spaces
figure()
% bode(sys_okid,'b', sys_n4sid,'r', sys_moesp, 'g', sys_orig , 'y');
% legend('OKID / ERA','N4SID', 'MOESP','Original');
bode(sys_okid, 'b', sys_orig , 'r', sys_brunton, 'black');
legend('OKID / ERA','Original','Brunton');

figure()
bode(tmp_sys_okid, 'b', sys_orig , 'r');
legend('11-17-21 OKID / ERA','Original');

figure()
bode(sys_moesp,'b', sys_orig , 'r');
legend('MOESP','original')