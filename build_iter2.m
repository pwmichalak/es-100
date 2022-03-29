function build_iter2()  
    
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


%% OKID
function markovParams = okid_iter2(u,y,p,avg,stdev,skew,kurt)
    % OKID constants from https://people.duke.edu/~hpgavin/SystemID/References/Juang+Phan+etal-JGCD-1993.pdf
    % For description of new parametrs mean,std,skew,kurt, see below.
    
    % same as first iteration build
    m = size(u,2);
    k = size(u,1);
    q = size(y,2); 
    V = makeAugmentedDataV(u,y,p,m,k,q);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                    %
    % SECOND ITERATION BUILD                                             %
    %                                                                    %
    % The new aspect of the second iteration build is the use of         %
    % weighted least squares to extract the observer Markov parameters.  %
    % The supposed benefit of this iteration is that more information    %
    % specific to the batteries is used (i.e. the alleged noise in the   %
    % battery), so the extraction should be improved in theory.          %
    %                                                                    %
    % In particular, a Pearson system of random numbers will be          %
    % generated to simulate electrochemical battery noise, and the       % 
    % covariance of said data matrix will be used in the weighted least- %
    % squares procedure.                                                 %
    %                                                                    %
    % New parameters:                                                    %
    % - avg: mean value of the electrochemical noise distribution        %
    % - stdev: standard deviation of the electrochemical noise           %
    % - skew: skewness of the electrochemical noise distribution         %
    % - kurt: kurtosis of the electrochemical noise distribution         %
    %                                                                    %
    % NB: Some combinations of moments are not valid; in particular, the %
    % kurtosis must be greater than the square of the skewness plus 1.   %
    %                                                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yt = y'; sz = length(yt);
    P = pearsrnd(avg, stdev, skew, kurt, sz, sz);
    R = cov(P);
    Rinv = eye(length(R)) / R; % compute the inverse of the noise matrix
    observerMarkovParams = yt * Rinv * V' / (V * Rinv * V');
    
    % same as first iteration build
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


%% recover Markov Parameters; same as first iteration build (see first iteration build for comments)
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


%% Extract RC parameters from state space matrices; same as first iteration build (see first iteration build for comments)
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