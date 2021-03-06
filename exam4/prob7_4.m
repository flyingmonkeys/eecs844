% Problem 7, Exam 4

clear variables;

load p6.mat; % loads X1..X4

% Cycle through all four data sets
for a=1:4
    % Select data set to analyze
    switch a
        case 1
            A = X1;
        case 2
            A = X2;
        case 3
            A = X3;
        case 4
            A = X4;
    end

    [M,N] = size(A); % N = number of time samples
    p = M; % number of sensors
    n = N; % number of time samples

    % Form a reflection matrix and forward-backward estimate of R
    J = zeros(M,M);
    for idx=1:M
        J(M-idx+1,idx) = 1;
    end
    Rfb = (1/(2*N))*(A*ctranspose(A) + J*conj(A)*transpose(A)*J);

    % since R is square PSDH, it can be decomposed via eigen decomposition
    % since R is Hermitian, lambdas are real-valued
    % Get eigenvalues and eigenvectors of Rfb
    [V, lambda] = eig(Rfb);
    lambda = diag(lambda);
    l = real(flipud(sort(lambda)));
    
    % Calculate the BIC (Bayesian Information Criterion)
    BIC = zeros(p,1);
    for q=0:p-1
        sum1 = 0; % first summation term
        sum2 = 0; % second summation term
        for i=q+1:p
            sum1 = sum1 + log(l(i));
            sum2 = sum2 + (l(i)/(p-q));
        end
        log_lq = n * (sum1 - (p-q)*log(sum2));
        BIC(q+1) = (-2*log_lq) + (((q*((2*p)-q))+1)*log(n));
    end

    [min_BIC, min_index] = min(BIC);
    fprintf('Estimated number of signals = %d\n',min_index-1);
end

