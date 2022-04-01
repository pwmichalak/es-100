%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
% First Iteration Project Build: Experiment 3                        %
%                                                                    %
% This experiment involves performing parameter extraction on        %
% an open source laboratory data set using the state space           % 
% identification method described as the build of this project       % 
% report.                                                            %
%                                                                    %
% The first iteration of the project build will be used to perform   %
% this experiment.                                                   %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
close all;

% the data from Wang et al. 2017 (included in folder) was preprocessed by hand in excel.
dst_channel = get_data('data/dst-Wang-Liu-Pan-Chen-2017.xlsx','Sheet1');
dst_channel.SampleTime = [];

% define input / output data
y_dst = dst_channel.voltage_V;
u_dst = dst_channel.current_A;
t_dst = dst_channel.timestep;

% perform subspace identification
x = 1:1000;
y_dst_train = y_dst(x); 
u_dst_train = u_dst(x); 
t_dst_train = t_dst(x);
p_dst = 15;
order = 2;
fs = 1; Ts = 1/fs;
[r0p, r1p, r2p, c1p, c2p, sysc_build, sysd_build, mse] = build_iter1(u_dst_train,...
                                                         y_dst_train,...
                                                         p_dst,...
                                                         order,...
                                                         Ts);

Adst = sysd_build.A;
Bdst = sysd_build.B;
Cdst = sysd_build.C;
Ddst = sysd_build.D;  

yr_dst_train = dlsim(Adst, Bdst, Cdst, Ddst, u_dst_train);
OCV_dst_train = y_dst_train(1);
                                                     
% plot training results
font = 20; % fontsize
figure(); 
axes('FontSize', font, 'NextPlot', 'add');
plot(y_dst_train,'DisplayName','Original Voltage','LineWidth',1); 
hold on;
plot(OCV_dst_train - yr_dst_train,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Training Reconstruction of Sampled Voltage From DST Data','FontSize', font); 

figure(); 
axes('FontSize', font, 'NextPlot', 'add');
yyaxis left; plot(y_dst_train,'DisplayName','Original Voltage','LineWidth',1); 
hold on; yyaxis right; 
plot(OCV_dst_train - yr_dst_train,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Training Reconstruction of Sampled Voltage From DST Data (Different Axes)','FontSize', font); 
    
                                                     
% obtain values for testing technical specification                                                    
xtest = 1001:4000;
y_dst_test = y_dst(xtest); 
u_dst_test = u_dst(xtest); 
t_dst_test = t_dst(xtest);

yr_dst_test = dlsim(Adst, Bdst, Cdst, Ddst, u_dst_test);
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
title('Testing Reconstruction of Sampled Voltage From DST Data','FontSize', font); 

figure(); 
axes('FontSize', font, 'NextPlot', 'add');
yyaxis left; plot(y_dst_test,'DisplayName','Original Voltage','LineWidth',1); 
hold on; yyaxis right; 
plot(OCV_dst_test - yr_dst_test,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Testing Reconstruction of Sampled Voltage From DST Data (Different Axes)','FontSize', font); 
    
