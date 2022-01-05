function [Ar, Br, Cr, Dr] = tmp_era(markovParams, sz)
    % Implement the Eigensystem Realization Algorithm to obtain
    %
    % Arguemnts
    % markovParams: a vector containing the system Markov Parameters D, CB,
    %   CAB, CA^2B, etc.
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
    [H0, H1] = hankelMatrices(markovParams(2:end),q,m,p);
    
    % obtain SVD of hankel matrices
    [U0,S0,V0] = svd(H0);
    
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



%% 
function [H0,H1] = hankelMatrices(M,q,m,p)
    % Create Hankel matrices H0 and H1, using notation from Juang and Phan,
    % "Identification and Control of Mechanical Systems." In particular, H0
    % shape is a x b where a and b are such that a x m > b x q, for m
    % number of inputs and q number of outputs. 
    %
    % Arguments
    % M: matrix of values to use to construct Hankel matrices
    % m: number of inputs
    % q: number of outputs
    % p: number of parameters in our Markov system
    
    % define a and b; for convenience, a+b <= p 
    a = round((q * p + 1) / (m + q)) + 1;
    b = p - a;
    
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