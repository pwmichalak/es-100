function [r0p, r1p, r2p, c1p, c2p, sysc_build] = build_iter3(u,y,p,order,Ts)
    % First Iteration Build
    % 
    % INPUTS
    % u: vector of current inputs to state space
    % y: vector of voltage outputs from state space
    % p: memory factor
    % order: order of the system being reconstructed (2 for a 2RC ECM)
    % Ts: sampling period
    %
    % OUTPUTS:
    % r0p: prediction for resistance r0 from the battery 2RC ECM
    % r1p: prediction for resistance r1 from the battery 2RC ECM
    % r2p: prediction for resistance r2 from the battery 2RC ECM
    % c1p: prediction for capacitance c1 from the battery 2RC ECM
    % c2p: prediction for capacitance c2 from the battery 2RC ECM
    
    % obtain estimated system state space matrices
    markovParams = okid(u,y,p);
    [Ar,Br,Cr,Dr] = era(markovParams, order);
    
    % convert reconstructed OKID state space matrices from discrete to continuous time 
    sysd_build = ss(Ar,Br,Cr,Dr,Ts);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                    %
    % THIRD ITERATION BUILD                                              %
    %                                                                    %
    % The main innovation of this iteration of the third iteration of    %
    % the project build is the utilization of the Tustin or bilinear     %
    % method of conversion from the discrete time to the continuous time.%
    %                                                                    %
    % Where the ZOH method is arguably the most often used method and    %
    % achieves a discretization that perfectly matches the system in the % 
    % time-domain (i.e. the true system at different points in time),    %
    % the Tustin method achieves a discretization that best approximates %
    % the system in the frequency domain (i.e. the frequency behavior of %
    % the system is best approximated by the Tustin method).             %
    %                                                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sysc_build = d2c(sysd_build,'tustin');
    
    Aco = sysc_build.A;
    Bco = sysc_build.B;
    Cco = sysc_build.C;
    Dco = sysc_build.D;
    
    % obtain RC parameter estimates
    [r0p, r1p, r2p, c1p, c2p, ~] = extract_2rc_params(Aco, Bco, Cco, Dco);
end


%% ERA; same as first iteration build (see first iteration build for comments)
function [Ar, Br, Cr, Dr] = era(markovParams, sz)
    Dr = markovParams{1};
    m = size(Dr,2);
    q = size(Dr,1);
    p = size(markovParams,2) - 1;
    [H0, H1] = hankelMatrices(markovParams(2:end),p,sz);
    [U0,S0,V0] = svd(H0);
%     system_order_strength = diag(S0)
%     p = char(inputdlg('Plot system singular value magnitudes (y/n)?')); 
%     if p == 'y'
%         figure();
%         semilogy(system_order_strength);
%         title('System Order Strength Plot');
%         xlabel('Order (Singular Value Index)');
%         ylabel('Magnitude of Order Strength');
%         grid on;
%     end
    S0 = S0(1:sz,1:sz);
    n = size(S0,1);
    U0 = U0(:,1:n);
    V0 = V0(:,1:n);
    tmp_B = S0^(1/2) * transpose(V0);
    Br = tmp_B(:,1:m);
    tmp_C = U0 * S0^(1/2);
    Cr = tmp_C(1:q,:);
    Ar = S0^(-1/2) * transpose(U0) * H1 * V0 * S0^(-1/2);
end


%% Create Hankel Matrices; same as first iteration build (see first iteration build for comments)
function [H0,H1] = hankelMatrices(M,p,sz)
    b = sz;
    a = p - sz;
    H0 = [];
    for i=1:a
        row = [];
        for j=1:b
            row = [row M{i + (j-1)}];
        end
        H0 = [H0; row];
    end
    H1 = [];
    for i=2:(a+1)
        row = [];
        for j=2:(b+1)
            row = [row M{(i-1) + (j-1)}];
        end
        H1 = [H1; row];
    end
end


%% OKID; same as first iteration build (see first iteration build for comments
function markovParams = okid(u,y,p)
    m = size(u,2);
    k = size(u,1);
    q = size(y,2); 
    V = makeAugmentedDataV(u,y,p,m,k,q);
    yt = y';
    observerMarkovParams = yt * V' / (V * V');
    markovParams = recoverMarkovParams(observerMarkovParams, m, q, p);
end


%% augmented data matrix; same as first iteration build (see first iteration build for comments)
function V = makeAugmentedDataV(u,y,p,m,n,q)
    V = zeros((m+q) * p + m,n);
    V(1:m,:) = u';
    for i = 2:p+1
        v = [u(1:end - i + 1,:)' ; y(1:end - i + 1,:)'];
        V((m+q)*(i-2) + m+1 : (m+q)*(i-1) + m, i : n) = v;
    end
end


%% recover system Markov Parameters; same as first iteration build (see first iteration build for comments)
function markovParams = recoverMarkovParams(observerParams, m, q, p)
    markovParams = cell(1,p+1);
    D = observerParams(:,1:m);
    markovParams{1} = D;
    for k=0:(p-1)
        Ybar_k_1 = observerParams(:, m + (q+m)*k + 1: m + (q+m)*k + m);
        Ybar_k_2 = observerParams(:, m + (q+m)*k + m + 1: m + (q+m)*k + (q+m));
        conv = 0;
        for i=0:(k-1)
            Ybar_i_2 = observerParams(:,m + (q+m)*i + m + 1:m + (q+m)*i + (q+m));
            conv = conv + Ybar_i_2 * markovParams{k - i + 1};
        end
        markovParams{k+2} = Ybar_k_1 + Ybar_k_2 * D + conv;
    end
end


%% extract RC parameters from state space matrices; same as first iteration build (see first iteration build for comments)
function [R0, R1, R2, Cap1, Cap2, T] = extract_2rc_params(A,B,C,D)
    assert(all(size(A)==[2 2]))
    assert(all(size(B)==[2 1]))
    assert(all(size(C)==[1 2]))
    assert(all(size(D)==[1 1]))

    R0 = D;

    tau_plus = (-trace(A) + sqrt(trace(A)^2 - 4 * det(A))) / (2 * det(A));
    tau_minus = (-trace(A) - sqrt(trace(A)^2 - 4 * det(A))) / (2 * det(A));
    tau1 = max(tau_plus, tau_minus);
    tau2 = min(tau_plus, tau_minus);

    T11 = 1 / C(1) - (C(2) * A(1,1) * tau1 + C(2)) / ((C(1) * C(2) * A(1,1) - C(1)^2 * A(1,2)) * tau1 + C(1) * C(2));
    T12 = 1 / C(1) - (C(2) * A(1,1) * tau2 + C(2)) / ((C(1) * C(2) * A(1,1) - C(1)^2 * A(1,2)) * tau2 + C(1) * C(2));
    T21 = (A(1,1) * tau1 + 1) / ((C(2) * A(1,1) - C(1) * A(1,2)) * tau1 + C(2));
    T22 = (A(1,1) * tau2 + 1) / ((C(2) * A(1,1) - C(1) * A(1,2)) * tau2 + C(2));
    T = [T11 T12; T21 T22];

    Cap1 = (T12 * T21 - T22 * T11) / (B(2) * T12 - B(1) * T22);
    Cap2 = (T12 * T21 - T22 * T11) / (B(1) * T21 - B(2) * T11);

    R1 = tau1 / Cap1;
    R2 = tau2 / Cap2;
end