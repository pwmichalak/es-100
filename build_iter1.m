function build_iter1()
    
end


%% ERA
function [Ar, Br, Cr, Dr] = era(markovParams, sz)
    % Implement the Eigensystem Realization Algorithm to obtain
    %
    % Arguemnts
    % markovParams: a vector containing the system Markov Parameters D, CB,
    %   CAB, CA^2B, etc.
    % sz: the order of the system we are approximating. 
    % 
    % Outputs
    % Ar: reduced order state matrix A
    % Br: reduced order input-to-state matrix B
    % Cr: reduced order state-to-output matrix C
    % D: feedthrough matrix D
    
    % extract D immediately
    Dr = markovParams{1};
    m = size(Dr,2);
    q = size(Dr,1);
    p = size(markovParams,2) - 1;
    
    % create hankel matrices 
    [H0, H1] = hankelMatrices(markovParams(2:end),p,sz);
    
    % obtain SVD of hankel matrices
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
    
    % reduce the dimensions down to obtain a square S0
    S0 = S0(1:sz,1:sz);
    n = size(S0,1);
    
    % obtain the first n columns of U and V
    U0 = U0(:,1:n);
    V0 = V0(:,1:n);
    
    % extract B, the first m columns of S0^(1/2) * V^T
    tmp_B = S0^(1/2) * transpose(V0);
    Br = tmp_B(:,1:m);
    
    % extract C, the first q rows of U0 *S0^(1/2)
    tmp_C = U0 * S0^(1/2);
    Cr = tmp_C(1:q,:);
    
    % extract A, basically just an unfolding of svd(H1)
    Ar = S0^(-1/2) * transpose(U0) * H1 * V0 * S0^(-1/2);
end


%% Create Hankel Matrices (ERA helper function)
function [H0,H1] = hankelMatrices(M,p,sz)
    % Create Hankel matrices H0 and H1, using notation from Juang and Phan,
    % "Identification and Control of Mechanical Systems." In particular, H0
    % shape is a x b where a and b are such that a x m > b x q, for m
    % number of inputs and q number of outputs. 
    %
    % Arguments
    % M: matrix of values to use to construct Hankel matrices
    % p: number of parameters in our Markov system
    % sz: the order of the system we are approximating. I think this is
    %   what b has to be equal to at least because it defines how many
    %   singular values there are, and we need at least one for each order
    %   of the system
    
    % define a and b; for convenience, a+b <= p 
    b = sz;
    a = p - sz;
    
    % construct Hankel matrix H0
    H0 = [];
    for i=1:a
        row = [];
        for j=1:b
            row = [row M{i + (j-1)}];
        end
        H0 = [H0; row];
    end
    
    % construct Hankel matrix H1
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
function markovParams = okid(u,y,p)
    % OKID constants from https://people.duke.edu/~hpgavin/SystemID/References/Juang+Phan+etal-JGCD-1993.pdf
    % 
    % p: number of time steps back that we go for prediction -- like the 
    %   forgetting factor alpha 
    
    % 1. extract input size, output size, and number of data points
    % m - number of inputs; u should be organized so that individual inputs
    %   are column vectors
    m = size(u,2);
    
    % k - number of data samples (ell); u should be organized so that 
    %   individual inputs are column vectors
    k = size(u,1);
    
    % q - number of outputs; y should be a column vector (each column is an
    %   output)
    q = size(y,2); 

    % 2. form the (zero initial condition) data matrix V of augmented data
    V = makeAugmentedDataV(u,y,p,m,k,q);
    
    % 3. compute the least-squares solution of the observer Markov
    % parameters (Ybar in Eq (7)).
    % Shape = q x ((m + q) * p + m);
    yt = y';
    observerMarkovParams = yt * V' / (V * V');
    
    % 4. Recover the combined system and Kalman filter Markov parameters
    markovParams = recoverMarkovParams(observerMarkovParams, m, q, p);
end


%% Obtain augmented (preprocessed) data matrix V (OKID helper function)
function V = makeAugmentedDataV(u,y,p,m,n,q)

    % u - a column vector of input data
    % y - a column vector of output data
    % p - the factor which defines how many steps back in time to go 
    % m - number of inputs = size(u,2)
    % n - number of data samples (ell) = size(u,1) = size(y,1)
    % q - number of outputs = size(y,2)
    
    % The augmented data matrix V contains the data used to predict the
    % output value at the current time step. The data in column j of the
    % augmented data matrix contains the input data at time step j-1 in
    % the first row, and the remaining rows contain the pairings of input
    % data and output data at time steps j-k for k = 2 to j. We only look p
    % time steps back to make our prediction (where the larger p gets, the
    % more the current estimate depends on previous values and the less
    % "Markovian" the system becomes (?)). 
    
    % initialize V to be a matrix of zeros
    V = zeros((m+q) * p + m,n);
    
    % populate the first rows of V with the current input data
    V(1:m,:) = u';
    
    % populate the remaining rows of V with the augmented past input-output
    % data vectors
    for i = 2:p+1
        
        % define augmented data vector v
        v = [u(1:end - i + 1,:)' ; y(1:end - i + 1,:)'];
        
        % put the augmented data vector in the augmented data matrix
        V((m+q)*(i-2) + m+1 : (m+q)*(i-1) + m, i : n) = v;
        
    end
end


%% Recover observer Markov Parameters (OKID helper function)
function markovParams = recoverMarkovParams(observerParams, m, q, p)
    
    % observerParams - a vector containing the observer Markov parameters
    %    identified from the least-squares solution of the observer problem
    % m - number of input parameters
    % q - number of output parameters
    % p - the factor which defines how many steps back in time to go 
    
    % Obtain the system and observer Markov parameters at the same time, as
    % given by eq. (29) in https://people.duke.edu/~hpgavin/SystemID/
    % References/Juang+Phan+etal-JGCD-1993.pdf
    
    % initialize the markovParameters matrix 
    markovParams = cell(1,p+1);
    
    % obtain the D parameter (first element of the observer Markov parameters)
    D = observerParams(:,1:m);
    markovParams{1} = D;
    
    % obtain the remaining Markov parameters iteratively
    for k=0:(p-1)
        % extract the coordinates from the observer Markov parameter k
        Ybar_k_1 = observerParams(:, m + (q+m)*k + 1: m + (q+m)*k + m);
        Ybar_k_2 = observerParams(:, m + (q+m)*k + m + 1: m + (q+m)*k + (q+m));
        
        % obtain the convolution of previous Markov parameters and observer
        % parameters
        conv = 0;
        for i=0:(k-1)
            Ybar_i_2 = observerParams(:,m + (q+m)*i + m + 1:m + (q+m)*i + (q+m));
            conv = conv + Ybar_i_2 * markovParams{k - i + 1};
        end
        
        % use the first coordinate from eq. (29) to obtain the Markov
        % parameters
        markovParams{k+2} = Ybar_k_1 + Ybar_k_2 * D + conv;
    end
end


%% Extract RC parameters from state space matrices
function [R0, R1, R2, Cap1, Cap2, T] = extract_2rc_params(A,B,C,D)
    % return the resistor and capacitor values in a 2RC equivalent circuit
    % battery model. Also return the transformation matrix T that converts
    % between the obtained state space and the desired 2RC state space

    % check that all matrices are of the correct size
    assert(all(size(A)==[2 2]))
    assert(all(size(B)==[2 1]))
    assert(all(size(C)==[1 2]))
    assert(all(size(D)==[1 1]))

    % perform parameter extraction
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