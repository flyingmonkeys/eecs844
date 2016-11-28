% Problem 8, Exam 4

clear variables;
close all;
hold off;

load p6.mat; % loads X1..X4
M_hat = 12;  % sub-array size
    
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

    % Perform spatial smoothing on the received signals
    Rss = zeros(M_hat,M_hat);
    for m=1:K
        % Create Xm matrix
        Xm = A(m:m+M_hat-1,:);

        % Calculate autocorrelation matrix R
        Rss = Rss + Xm*ctranspose(Xm);
    end
    Rss = Rss / (L*K);

    % since R is square PSDH, it can be decomposed via eigen decomposition
    % since R is Hermitian, lambdas are real-valued
    % Get eigenvalues and eigenvectors of Rss
    [V, lambda] = eig(Rss);
    lambda = diag(lambda);
    l = real(flipud(sort(lambda)));
    fprintf('Condition number = %d\n',cond(Rss));
    
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
    
    % Compute MVDR
    R_inv = inv(Rss);
    numsamps = 20*M_hat;
    theta = linspace(-pi,pi,numsamps);
    degrees = linspace(-180,180,numsamps);
    MVDR = zeros(numsamps,1);
    for idx=1:numsamps
        % Compute steering matrix (electrical angle)
        s = transpose(exp(-j*theta(idx)*linspace(0,M_hat-1,M_hat)));
        MVDR(idx) = 1.0 / ( ctranspose(s)*R_inv*s );
    end

    % MVDR Power Spectrum (electrical angle)
    plot(degrees,20*log10(abs(MVDR)));
    hold on;
end

title('MVDR Power Spectrum');
xlabel('Electrical Theta (deg)');
ylabel('Magnitude (dB)');
grid on;
legend({'X1','X2','X3','X4'});
