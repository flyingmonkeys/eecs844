% Exam 4 Problem 1

clear all;
close all;
hold off;

load p1.mat; % loads x,d

M = 60; % Wiener filter order
K = length(x);
N = K-M+1;

mu_hat = 0.5;  % step size
delta  = 0.03; % leakage factor

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

%% Normalized LMS
w   = zeros(M,1);           % init NLMS weights to zero
e_squared   = zeros(K,1);   % squared error
squared_dev = zeros(K,1);   % squared deviation

for n=M:K
    y_l = ctranspose(w)*flipud(x(n-M+1:n));     % compute LMS filter output
    e = d(n) - y_l;                             % compute error
    x_norm = ctranspose(x(n-M+1:n))*x(n-M+1:n); % get input energy for normalization
    mu = mu_hat / (delta + x_norm);             % compute step size
    w = w + mu * conj(e) * flipud(x(n-M+1:n));  % LMS update
    
    e_squared(n) = conj(e)*e;                   % squared error
    squared_dev(n) = ctranspose(w-wo)*(w-wo);   % squared deviation
end

%% Plots----------------------------

figure(1);
plot(20*log10(e_squared));
title('Squared Error (NLMS)');
xlabel('Iteration number');
ylabel('squared error (dB)');
grid on;

figure(2);
plot(20*log10(squared_dev));
title('Squared Deviation (NLMS)');
xlabel('Iteration number');
ylabel('squared deviation (dB)');
grid on;
