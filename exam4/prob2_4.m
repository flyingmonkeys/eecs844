% Exam 4 Problem 2

clear all;
close all;
hold off;

load p1.mat; % loads x,d

M = 60; % Wiener filter order
K = length(x);
N = K-M+1;

mu_hat_lms = 0.5;  % step size
mu_hat_KLT = 0.125/M;  % step size
delta      = 0.03;     % leakage factor

%% Wiener filter---------------------

% Create X snapshot matrix
X = zeros(M,K);
for k=1:N
    X(:,k) = flipud(x(k:k+M-1));
end

% Calculate autocorrelation matrix R
Rxx = (1/N)*X*ctranspose(X);

% Calculate P
P = zeros(M,1);
for k=M:K
    P = P + (flipud(x(k-M+1:k)) * conj(d(k)));
end
P = P / N;

% Calculate Wiener filter (wo)
wo = inv(Rxx)*P;

%% KLT transform LMS and Normalized LMS ---------------------
[Q,D] = eig(Rxx);
D = D + delta*eye(size(D));

% KLT initialization
wt              = zeros(M,1);   % init weights to zero
w_klt           = zeros(M,1);   % init weights to zero
e_squared_klt   = zeros(K,1);   % squared error
squared_dev_klt = zeros(K,1);   % squared deviation

% NLMS initialization
w_lms   = zeros(M,1);           % init NLMS weights to zero
e_squared_lms   = zeros(K,1);   % squared error
squared_dev_lms = zeros(K,1);   % squared deviation

for n=M:K
    % Do KLT transform LMS-------------
    y = ctranspose(w_klt)*flipud(x(n-M+1:n));       % compute KLT filter output
    e = d(n) - y;                                   % compute error
    z = ctranspose(Q)*flipud(x(n-M+1:n));           % compute transformed input vector
    wt = wt + (mu_hat_KLT * inv(D) * z * conj(e));  % compute transformed weight vector update
    w_klt = Q*wt;                                   % re-map weight vector back to "normal" domain
    
    e_squared_klt(n)   = conj(e)*e;                         % KLT squared error
    squared_dev_klt(n) = ctranspose(w_klt-wo)*(w_klt-wo);   % KLT squared deviation
    
    % Do NLMS------------------------
    y = ctranspose(w_lms)*flipud(x(n-M+1:n));   % compute LMS filter output
    e = d(n) - y;                               % compute error
    x_norm = ctranspose(x(n-M+1:n))*x(n-M+1:n); % get input energy for normalization
    mu = mu_hat_lms / (delta + x_norm);         % compute step size
    w_lms = w_lms + mu * conj(e) * flipud(x(n-M+1:n));  % LMS update
    
    e_squared_lms(n) = conj(e)*e;                           % NLMS squared error
    squared_dev_lms(n) = ctranspose(w_lms-wo)*(w_lms-wo);   % NLMS squared deviation
end

%% Plots----------------------------

figure(1);
plot(20*log10(e_squared_klt));
hold on;
plot(20*log10(e_squared_lms));
title('Squared Error (KLT vs. NLMS)');
xlabel('Iteration number');
ylabel('squared error (dB)');
grid on;
legend({'KLT','NLMS'});

figure(2);
plot(20*log10(squared_dev_klt));
hold on;
plot(20*log10(squared_dev_lms));
title('Squared Deviation (KLT vs. NLMS)');
xlabel('Iteration number');
ylabel('squared deviation (dB)');
grid on;
legend({'KLT','NLMS'});
