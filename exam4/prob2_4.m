% Exam 4 Problem 2

clear all;
close all;
hold off;

load p1.mat; % loads x,d

M = 60; % Wiener filter order
K = length(x);

mu_hat = 0.125/M; % step size
delta = 0.03;     % leakage factor

%% Wiener filter---------------------
X = zeros(M,K);
% Create X matrix
for k=1:K-M+1
    X(:,k) = flipud(x(k:k+M-1));
end

% Calculate autocorrelation matrix R
Rxx = (1/(K-M+1))*X*ctranspose(X);

% Calculate P
P = zeros(M,1);
for k=M:K
    P = P + (flipud(x(k-M+1:k)) * conj(d(k)));
end
P = P / (K-M+1);

% Calculate Wiener filter (ww)
wo = inv(Rxx)*P;

%% KLT transform LMS---------------------
[Q,D] = eig(Rxx);
D = D + delta*eye(size(D));

y = zeros(K,1);
wt = zeros(M,1);
w = zeros(M,1);
e_squared = zeros(K,1);
squared_dev = zeros(K,1);
for n=M+1:K
    y(n) = ctranspose(w)*flipud(x(n-M+1:n));
    e = d(n) - y(n);
    
    z = ctranspose(Q)*flipud(x(n-M+1:n));
    wt = wt + (mu_hat * inv(D) * z * conj(e));
    w = Q*wt;
    
    e_squared(n) = conj(e)*e;
    squared_dev(n) = ctranspose(w-wo)*(w-wo);
end

%% Plots----------------------------
figure(1);

plot(abs(d));
hold on
plot(abs(y),'color','red');

figure(2);
plot(20*log10(e_squared));
title(‘Squared Error’);
xlabel(‘Iteration number’);
ylabel(‘squared error (dB)’);

figure(3);
plot(20*log10(squared_dev));
title(‘Squared Deviation’);
xlabel(‘Iteration number’);
ylabel(‘squared deviation (dB)’);

