%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
% Second Iteration Project Build: Experiment 4                       %
%                                                                    %
% This experiment involves performing parameter extraction on        %
% an open source real-world driving data set using the state space   % 
% identification method described as the build of this project       % 
% report.                                                            %
%                                                                    %
% The second iteration of the project build will be used to perform  %
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
OCV_udds_train = y_udds_train(1);

xtest = 1001:4000;
y_udds_test = y_udds(xtest); 
u_udds_test = u_udds(xtest); 
t_udds_test = t_udds(xtest);
OCV_udds_test = y_udds_test(1);

p_udds = 56;
order = 2;
fs = 1; Ts = 1/fs;

% perform several trials of identification
ntrials = 10000;
mses = zeros(1,ntrials);
mses_test = zeros(1,ntrials);
min_mse = Inf; min_mse_test = Inf; 
yr_udds_train_best = y_udds_train;
yr_udds_test_best = y_udds_test;
for i = 1:ntrials
    disp(i);
    [r0p, r1p, r2p, c1p, c2p, sysc_build, sysd_build, mse] = build_iter2(u_udds_train,...
                                                             y_udds_train,...
                                                             p_udds,...
                                                             order,...
                                                             Ts);

    Audds = sysd_build.A;
    Budds = sysd_build.B;
    Cudds = sysd_build.C;
    Dudds = sysd_build.D;  

    yr_udds_train = dlsim(Audds, Budds, Cudds, Dudds, u_udds_train);

    % obtain values for testing technical specification                                                    
    yr_udds_test = dlsim(Audds, Budds, Cudds, Dudds, u_udds_test);
    mse_test = immse(y_udds_test,OCV_udds_test - yr_udds_test);
    
    % get mse values
    mses(i) = mse;
    mses_test(i) = mse_test;
    if mse < min_mse
        min_mse = mse;
        yr_udds_train_best = yr_udds_train;
    end
    if mse_test < min_mse_test
        min_mse_test = mse_test;
        yr_udds_test_best = yr_udds_test;
    end
end; clear i;

% replace nan and inf values in the trials to obtain statistics
mses(isinf(mses) | isnan(mses)) = realmax; 
mses_test(isinf(mses_test) | isnan(mses_test)) = realmax; 
                              
% plot training results
font = 20; % fontsize
figure(); 
axes('FontSize', font, 'NextPlot', 'add');
plot(y_udds_train,'DisplayName','Original Voltage','LineWidth',1); 
hold on;
plot(OCV_udds_train - yr_udds_train_best,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Training Reconstruction of Sampled Voltage From UDDS Data','FontSize', font); 

figure(); 
axes('FontSize', font, 'NextPlot', 'add');
yyaxis left; plot(y_udds_train,'DisplayName','Original Voltage','LineWidth',1); 
hold on; yyaxis right; 
plot(OCV_udds_train - yr_udds_train_best,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Training Reconstruction of Sampled Voltage From UDDS Data (Different Axes)','FontSize', font); 

% plot testing results
figure(); 
axes('FontSize', font, 'NextPlot', 'add');
plot(y_udds_test,'DisplayName','Original Voltage','LineWidth',1); 
hold on;
plot(OCV_udds_test - yr_udds_test_best,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Testing Reconstruction of Sampled Voltage From UDDS Data','FontSize', font); 

figure(); 
axes('FontSize', font, 'NextPlot', 'add');
yyaxis left; plot(y_udds_test,'DisplayName','Original Voltage','LineWidth',1); 
hold on; yyaxis right; 
plot(OCV_udds_test - yr_udds_test_best,'DisplayName','Estimated Voltage','LineWidth',1); 
legend('FontSize', font);
xlabel('Time (s)','FontSize', font); 
ylabel('Voltage (V)','FontSize', font);
title('Testing Reconstruction of Sampled Voltage From UDDS Data (Different Axes)','FontSize', font); 
    
