function [u,y,order] = simulate(A,B,C,D,Ts,n)
    % obtain random input data current u and output data voltage y
    % 
    % INPUTS
    % A,B,C,D: system state space matrices
    % Ts: sampling period
    % n: number of data points
    % 
    % OUTPUTS
    % u: generated input current data vector
    % y: simulated output voltage data vector
    % order: predicted order of the model
    
    q = size(D,1); % output size
    m = size(D,2); % input size
    order = size(A,2); % order of the model

    % make sure the matrices are consistent in their dimensions
    assert(m == size(B,2)); 
    assert(q == size(C,1));
    assert(order == size(C,2));
    assert(order == size(B,1));

    % simulate data
    nu = n / Ts;
    u = randn(nu ,m); %rms(u) ~ 1
    sysc = ss(A,B,C,D); % define continuous time state space
    sysd = c2d(sysc,Ts); % convert state space to discrete time

    % output; we have to scale u so that rms(y) ~ [0.005, 0.01]
    % this is consistent with the literature for what the voltage values
    % of a battery should be.
    y = dlsim(sysd.A,sysd.B,sysd.C,sysd.D,u);
    rms_scale = ((randn(1) + 1) / 2) * 0.005 + 0.005; % random value between 0.005 and 0.01
    scale = rms_scale / rms(y); 
    u = u * scale;
    y= dlsim(sysd.A,sysd.B,sysd.C,sysd.D,u); 
end
