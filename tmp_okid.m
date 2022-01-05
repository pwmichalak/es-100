function markovParams = tmp_okid(u,y,p)
    % OKID constants from https://people.duke.edu/~hpgavin/SystemID/References/Juang+Phan+etal-JGCD-1993.pdf
    % 
    % ell: n (number of data samples)
    % q: 1 (number of outputs; 1 here for voltage output — note this is 
    %    different from the number of state outputs, which is 2 in this 
    %    case for the voltage and current)
    % m: 1 (number of inputs; 1 here for voltage input)
    % p: TBD (number of time steps back that we go for prediction -- like 
    %    the forgetting factor)

    %% 1. extract input size, output size, and number of data points
    
    % m - number of inputs; u should be organized so that individual inputs
    %   are column vectors
    m = size(u,2);
    
    % n - number of data samples (ell); u should be organized so that 
    %   individual inputs are column vectors
    k = size(u,1);
    
    % q - number of outputs; y should be a column vector
    q = size(y,2); 

    %% 2. form the (zero initial condition) data matrix V of augmented data
    V = makeAugmentedDataV(u,y,p,m,k,q);
    
    %% 3. compute the least-squares solution of the observer Markov
    % parameters (Ybar in Eq (7)).
    % Shape = q x ((m + q) * p + m); 
    observerMarkovParams = lsqr(V', y(:,1))';

    %% 4. Recover the combined system and Kalman filter Markov parameters
    markovParams = recoverMarkovParams(observerMarkovParams, m, q);
end




%%
function V = makeAugmentedDataV(u,y,p,m,n,q)

    % u - a column vector of input data
    % y - a column vector of output data
    % p - the factor which defines how many steps back in time we want to
    %   go 
    % m - number of inputs
    % n - number of data samples (ell)
    % q - number of outputs
    
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
    V(1:size(u,2),:) = u';
    
    % populate the remaining rows of V with the augmented past input-output
    % data vectors
    for i = 2:p
        
        % define augmented data vector v
        v = [u(1:end - i + 1,:)' ; y(1:end - i + 1,:)'];
        
        % put the augmented data vector in the augmented data matrix
        V((size(u,2) + size(y,2))*(i-1) : (size(u,2) + size(y,2))*(i)-1, i : n) = v;
        
    end
end



%%
function markovParams = recoverMarkovParams(observerParams, m, q)
    
    % observerParams - a vector containing the observer Markov parameters
    %    identified from the least-squares solution of the observer problem
    
    % Obtain the system and observer Markov parameters at the same time, as
    % given by eq. (29) in https://people.duke.edu/~hpgavin/SystemID/
    % References/Juang+Phan+etal-JGCD-1993.pdf
    
    % initialize the markovParameters matrix 
    p = ((size(observerParams,2) - m) / (m + q));
    markovParams = cell(1,p+1);
    
    % obtain the D parameter (first element of the observer parameters)
    D = observerParams(:,1);
    markovParams{1} = D;
    
    % obtain the remaining Markov parameters iteratively
    for k=0:(p-1)
        
        % extract the coordinates from the observer Markov parameter k
        Ybar_k_1 = observerParams(1+2*k+1);
        Ybar_k_2 = observerParams(1+2*k+2);
        
        % obtain the convolution of previous Markov parameters and observer
        % parameters
        conv = 0;
        for i=0:(k-1)
            Ybar_i_2 = observerParams(1+2*i+2);
            conv = conv + Ybar_i_2 * markovParams{k - i + 1};
        end
        
        % use the first coordinate from eq. (29) to obtain the Markov
        % parameters
        markovParams{k+2} = Ybar_k_1 + Ybar_k_2 * D + conv;
    end
end