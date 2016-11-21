% Exam 4 Problem 1

clear all;
close all;
hold off;

load p1.mat; % loads x,d

M = 60; % Wiener filter order
K = length(x);

mu_hat = 0.5;     % step size
delta = 0.03; % leakage factor

%% Wiener filter---------------------
% Calculate R
rc = xcorr(x,M-1,'unbiased');
Rxx = toeplitz(conj(flipud(rc(1:M))));

% Calculate P
P = zeros(M,1);
for k=M:K
    P = P + (flipud(x(k-M+1:k)) * conj(d(k)));
end
P = P / (K-M+1);

% Calculate Wiener filter (ww)
wo = inv(Rxx)*P;

% Run x through Wiener filter
y_w = zeros(K,1);
for n=M:K
    y_w(n) = ctranspose(wo)*flipud(x(n-M+1:n));
end

%% Normalized LMS

y_l = zeros(K,1);
w = zeros(M,1);
e_squared = zeros(K,1);
squared_dev = zeros(K,1);
for n=M+1:K
    y_l(n) = ctranspose(w)*flipud(x(n-M+1:n));
    e = d(n) - y_l(n);
    x_norm = ctranspose(x(n-M+1:n))*x(n-M+1:n);
    mu = mu_hat / (delta + x_norm);
    w = w + mu * conj(e) * flipud(x(n-M+1:n));
    
    e_squared(n) = conj(e)*e;
    squared_dev(n) = ctranspose(w-wo)*(w-wo);
end

%% Plots----------------------------
figure(1);

plot(abs(d));
hold on
plot(abs(y_w),'color','red');
plot(abs(y_l),'color','green');

figure(2);
plot(20*log10(e_squared));

figure(3);
plot(20*log10(squared_dev));
