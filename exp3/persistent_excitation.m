% check for persistent excitation (at every order) of client data
clear; close all;

% obtain client data
channel = get_data('data/dst-Wang-Liu-Pan-Chen-2017.xlsx','Sheet1');
channel.SampleTime = [];

fs = 1; % sampling frequency
Ts = 1/fs; % sampling period
y = channel.voltage_V;
u = channel.current_A;
t = channel.timestep;

% resample data to ensure that pulses are evenly spaced in time
[u, tu] = resample(u, t, Ts, 1, 1);
[y, ty] = resample(y, t, Ts, 1, 1);

% obtain correlation matrices
N = length(u);
stops = 1:200;
dets = zeros(1,length(stops));
j = 1;
for stop = stops
    Rs = zeros(1,stop);
    for i = 1:stop
        Rs(i) = (1/N) * dot(u(i:end), u(1:end-i+1));
    end
    Rsf = fi(Rs);
    R = str2num(toeplitz(Rsf).Value); %#ok<ST2NM>

    % check eignevalues of R matrix
    dets(j) = det(R);
    j = j + 1;
end
font = 20; figure(); axes('FontSize', font, 'NextPlot', 'add');
plot(dets, 'LineWidth',2); 
xlabel('Correlation Matrix Size', 'FontSize', font);
ylabel('Determinant', 'FontSize', font); 
title('Order of Persistent Excitation for DST data', 'FontSize', font);
