%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
% Third Iteration Project Build: Experiment 4                        %
%                                                                    %
% This experiment involves performing parameter extraction on        %
% an open source real-world driving data set using the state space   % 
% identification method described as the build of this project       % 
% report.                                                            %
%                                                                    %
% The third iteration of the project build will be used to perform   %
% this experiment.                                                   %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
close all;

% the data from Wang et al. 2017 (included in folder) was preprocessed by hand in excel.
udds_channel = get_data('data/udds-Wang-Liu-Pan-Chen-2017.xlsx','Sheet1');
udds_channel.SampleTime = [];

% define input / output data
y_udds = udds_channel.voltage_V;
u_udds = udds_channel.current_A;
t_udds = udds_channel.timestep;

% perform subspace identification
x = 1:1000;
y_udds_train = y_udds(x); 
u_udds_train = u_udds(x); 
t_udds_train = t_udds(x);
p_udds = 56;
order = 2;
fs = 1; Ts = 1/fs;
[r0p, r1p, r2p, c1p, c2p, sysc_build, sysd_build, mse] = build_iter3(u_udds_train,...
                                                         y_udds_train,...
                                                         p_udds,...
                                                         order,...
                                                         Ts);

Audds = sysd_build.A;
Budds = sysd_build.B;
Cudds = sysd_build.C;
Dudds = sysd_build.D;  

yr_udds_train = dlsim(Audds, Budds, Cudds, Dudds, u_udds_train);
OCV_udds_train = y_udds_train(1);
                                                     
% plot training results
font = 20; % fontsize
figure(); 
axes('FontSize', font, 'NextPlot', 'add');
plot(y_udds_train,'DisplayName','Original Voltage','LineWidth',1); 
hold on;
plot(OCV_udds_train - yr_udds_train,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Training Reconstruction of Sampled Voltage From UDDS Data','FontSize', font); 

figure(); 
axes('FontSize', font, 'NextPlot', 'add');
yyaxis left; plot(y_udds_train,'DisplayName','Original Voltage','LineWidth',1); 
hold on; yyaxis right; 
plot(OCV_udds_train - yr_udds_train,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Training Reconstruction of Sampled Voltage From UDDS Data (Different Axes)','FontSize', font); 
    
                                                     
% obtain values for testing technical specification                                                    
xtest = 1001:4000;
y_dst_test = y_udds(xtest); 
u_dst_test = u_udds(xtest); 
t_dst_test = t_udds(xtest);

yr_dst_test = dlsim(Audds, Budds, Cudds, Dudds, u_dst_test);
OCV_dst_test = y_dst_test(1);
mse_test = immse(y_dst_test,OCV_dst_test - yr_dst_test);

% plot testing results
figure(); 
axes('FontSize', font, 'NextPlot', 'add');
plot(y_dst_test,'DisplayName','Original Voltage','LineWidth',1); 
hold on;
plot(OCV_dst_test - yr_dst_test,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Testing Reconstruction of Sampled Voltage From UDDS Data','FontSize', font); 

figure(); 
axes('FontSize', font, 'NextPlot', 'add');
yyaxis left; plot(y_dst_test,'DisplayName','Original Voltage','LineWidth',1); 
hold on; yyaxis right; 
plot(OCV_dst_test - yr_dst_test,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Testing Reconstruction of Sampled Voltage From UDDS Data (Different Axes)','FontSize', font); 
    
