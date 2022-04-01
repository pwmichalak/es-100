%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
% Third Iteration Project Build: Experiment 2                        %
%                                                                    %
% This experiment consists of a sensitivity analysis of the          %
% parameter estimation method described in the project build. The    %
% purpose of this sensitivity analysis is to support the validation  %
% of the project build and present a description of the              %
% considerations that implementation of the project build would have %
% to consider.                                                       %
%                                                                    %
% The third iteration of the project build will be used to perform   %
% this experiment.                                                   %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Define 2RC Circuit elements and state space matrices
% 
% The state variable x: [capacitor 1 voltage ; capacitor 2 voltage]
% Input variable u: battery current
% Output variable y: battery voltage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r0 = rand(1,1); 
r1 = rand(1,1); 
r2 = rand(1,1); 
c1 = 1e3 * rand(1,1);
c2 = 1e3 * rand(1,1);

A = [-1/(r1 * c1) 0; 0 -1/(r2 * c2)];
B = [1/c1; 1/c2];
C = [1 1];
D = r0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Sensitivity 1: adding noise to the input and output data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n");
fprintf("Sensitivity 1: adding noise to the input and output data\n");
fprintf("--------------------------------------------------------------\n");
n = 1000;
fs = 1; Ts = 1/fs;
rms_valsy = 0:10e-6:80e-6;
rms_valsu = 0:1e-3:5e-3;
errs1 = zeros(5,length(rms_valsu),length(rms_valsy));
preds1 = zeros(5,length(rms_valsu),length(rms_valsy));
mses1 = zeros(length(rms_valsu),length(rms_valsy));
        
% simulate input and output data
[u,y,order] = simulate(A,B,C,D,Ts,n);

i = 1; j = 1;
for rmsy = rms_valsy
    for rmsu = rms_valsu
        fprintf(strcat("i=",num2str(i)," j=",num2str(j)),"\n");

        % add noise to data
        snr_y = 6 * rms(y) / rmsy;
        y_y = awgn(y, snr_y, 'measured');
        
        snr_u = 6 * rms(u) / rmsu;
        u_u = awgn(u, snr_u, 'measured');

        p = randi([4 5],1) * order; % memory factor; should be 4-5 times the order
        [r0p, r1p, r2p, c1p, c2p, sysc_build, sysd_build, mse] = build_iter3(u_u,y_y,p,order,Ts);

        % obtain percent error of parameter predictions 
        r0err = abs(r0p - r0) / r0;
        r1err = abs(max(r1p,r2p) - max(r1,r2)) / max(r1,r2);
        r2err = abs(min(r1p,r2p) - min(r1,r2)) / min(r1,r2);
        c1err = abs(max(c1p,c2p) - max(c1,c2)) / max(c1,c2);
        c2err = abs(min(c1p,c2p) - min(c1,c2)) / min(c1,c2);

        errs1(:,j,i) = [r0err;r1err;r2err;c1err;c2err];
        preds1(:,j,i) = [r0p;r1p;r2p;c1p;c2p];
        mses1(j,i) = mse;
        j = j+1;
    end
    i = i+1;
end; clear i; clear j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Sensitivity 2: changing the sampling rate of the data sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Sensitivity 2: changing the sampling rate of the data sets\n");
fprintf("--------------------------------------------------------------\n");

% the time constant is on the order of 1e2; we must not subceed beneath
% the cutoff frequency of the system, which is 1/2time constant = 5e-3
fsvals = [5e-2 1e-1 1 1e1];
n = 1000;
errs2 = zeros(5,length(fsvals));
preds2 = zeros(5,length(fsvals));
mses2 = zeros(1,length(fsvals));
i = 1;
for fs = fsvals
    % obtain sampling period 
    Ts = 1/fs;
    fprintf(strcat("fs=",num2str(fs),"; Ts=",num2str(Ts),"\n"));
    
    % generate input current data and output voltage data
    [u,y,order] = simulate(A,B,C,D,Ts,n);
    
    % obtain predicted parameter values
    p = randi([4 5],1) * order; % memory factor; should be 4-5 times the order
    [r0p, r1p, r2p, c1p, c2p, sysc_build, sysd_build, mse] = build_iter3(u,y,p,order,Ts);
    
    % obtain percent error of parameter predictions 
    r0err = abs(r0p - r0) / r0;
    r1err = abs(max(r1p,r2p) - max(r1,r2)) / max(r1,r2);
    r2err = abs(min(r1p,r2p) - min(r1,r2)) / min(r1,r2);
    c1err = abs(max(c1p,c2p) - max(c1,c2)) / max(c1,c2);
    c2err = abs(min(c1p,c2p) - min(c1,c2)) / min(c1,c2);
    
    errs2(:,i) = [r0err;r1err;r2err;c1err;c2err];
    preds2(:,i) = [r0p;r1p;r2p;c1p;c2p];
    mses2(i) = mse;
    i = i+1;
end; clear i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Sensitivity 3: changing the predicted system order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Sensitivity 3: changing the predicted system order\n");
fprintf("--------------------------------------------------------------\n");
n = 1000;
fs = 1; Ts = 1/fs;
orders = 2:10;
mses3 = zeros(1,length(orders));
        
% simulate input and output data
[u,y,~] = simulate(A,B,C,D,Ts,n);

i = 1;
for order = orders
    fprintf(strcat("i=",num2str(order)),"\n");

    p = randi([4 5],1) * order; % memory factor; should be 4-5 times the order
    [r0p, r1p, r2p, c1p, c2p, sysc_build, sysd_build, mse] = build_iter3(u,y,p,order,Ts);
    
    mses3(i) = mse;
    i = i+1;
end; clear i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Sensitivity 4: changing the memory factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Sensitivity 4: changing the memory factor\n");
fprintf("--------------------------------------------------------------\n");

plower = 10;
pupper = 500;
pvals = plower:pupper;

errs4 = zeros(5,length(pvals));
preds4 = zeros(5,length(pvals));
mses4 = zeros(1,length(pvals));

n = 5000;
fs = 1; Ts = 1/fs;

% generate input current data and output voltage data
[u,y,order] = simulate(A,B,C,D,Ts,n);

i = 1;
for p = pvals
    disp(i);
    
    % obtain predicted parameter values
    [r0p, r1p, r2p, c1p, c2p, sysc_build, sysd_build, mse] = build_iter3(u,y,p,order,Ts);
    
    % obtain percent error of parameter predictions 
    r0err = abs(r0p - r0) / r0;
    r1err = abs(max(r1p,r2p) - max(r1,r2)) / max(r1,r2);
    r2err = abs(min(r1p,r2p) - min(r1,r2)) / min(r1,r2);
    c1err = abs(max(c1p,c2p) - max(c1,c2)) / max(c1,c2);
    c2err = abs(min(c1p,c2p) - min(c1,c2)) / min(c1,c2);
    
    errs4(:,i) = [r0err;r1err;r2err;c1err;c2err];
    preds4(:,i) = [r0p;r1p;r2p;c1p;c2p];
    mses4(i) = mse;
    i = i+1;
end; clear i;

% process data so that extremely large values are removed from analysis
mses4(mses4 > 10) = NaN;

% plot mses curve
font = 20; % fontsize
figure(); axes('FontSize', font, 'NextPlot', 'add');
plot(pvals, mses4, 'LineWidth', 3);
hold on; xline(pvals(mses4==min(mses4)), 'LineWidth', 3, 'Color', 'r');
xlabel('Memory Factor', 'FontSize', font); 
ylabel('Mean-squared-error (MSE)', 'FontSize', font); 
title('MSE as a Function of Memory Factor, Iteration 3', 'FontSize', font);

% plot prediction errors curves
names = ["r0 "; "r1 "; "r2 "; "c1 "; "c2 "];
colors = ["red";"blue";"magenta";"#FD5602";"black"];
for i = 1:5
    figure();
    axes('FontSize', font, 'NextPlot', 'add');
    plot(pvals, errs4(i,:), 'LineWidth', 2, 'DisplayName', strcat(names(i),"Error"), 'Color', colors(i))
    xlabel('Memory Factor', 'FontSize', font); 
    ylabel(strcat(names(i),"Parameter Prediction Error (%)"), 'FontSize', font); 
    title(strcat(names(i),"Parameter Error as a Function of Memory Factor, Iteration 3"), 'FontSize', font);
    legend('FontSize',font);
end; clear i;