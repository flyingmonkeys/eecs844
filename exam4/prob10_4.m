% Problem 10, Exam 4

clear variables;
close all;
hold off;

load p6.mat; % loads X1..X4
order = [4 1 4 1; ...
         8 2 4 2; ...
         4 4 4 11]; % array of BIC-derived orders for each problem

M_hat = 12;
     
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

    n = N;           % number of time samples
    L = N;           % number of time samples
    K = M-M_hat+1;   % number of sub-arrays to average
    p = M_hat;       % number of sub-array sensors

    % Calculate autocorrelation matrix R
    Rxx = (1/N)*A*ctranspose(A);

    % Form a reflection matrix and forward-backward estimate of R
    J = zeros(M,M);
    for idx=1:M
        J(M-idx+1,idx) = 1;
    end
    Rfb = (1/(2*N))*(A*ctranspose(A) + J*conj(A)*transpose(A)*J);

    % Perform spatial smoothing on the received signals
    Rss = zeros(M_hat,M_hat);
    for m=1:K
        % Create Xm matrix
        Xm = A(m:m+M_hat-1,:);

        % Calculate autocorrelation matrix Rss
        Rss = Rss + Xm*ctranspose(Xm);
    end
    Rss = Rss / (L*K);

    % Choose autocorrelation matrix to use for analysis
    R = Rfb;
    %M = M_hat;
    
    % since R is square PSDH, it can be decomposed via eigen decomposition
    % since R is Hermitian, lambdas are real-valued
    % Get eigenvalues and eigenvectors of R
    [V, lambda] = eig(R);
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
    
    % Declare number of detected signals
    p = min_index-1;
    
    % Compute pseudo-spectrum
    numsamps = 20*M;
    theta = linspace(-pi,pi,numsamps);
    degrees = linspace(-180,180,numsamps);
    
    PS = zeros(numsamps,1);
    for idx=1:numsamps
        U = 0.0;
        for k=p+1:M
            % Compute steering matrix (electrical angle)
            s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
            U = U + abs(ctranspose(s)*V(:,k))^2;
        end
        PS(idx) = 1/U;
    end

    % MVDR Power Spectrum (electrical angle)
    plot(degrees,20*log10(abs(PS)));
    hold on;
end

title('MUSIC Pseudo-Spectrum');
xlabel('Electrical Theta (deg)');
ylabel('Magnitude (dB)');
grid on;
legend({'X1','X2','X3','X4'});

figure(2)
pmusic(ctranspose(X1),4); % comparison MUSIC from MATLAB built-in function
