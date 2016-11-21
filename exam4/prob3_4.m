% Exam 4 Problem 3

clear all;
close all;
hold off;

load p1.mat; % loads x,d

M = 60; % Wiener filter order
K = length(x);

lambda_array = [0.99 0.998 0.9995];

%% Wiener filter---------------------
% Calculate R
rc = xcorr(x,M-1,'unbiased');
Rxx = toeplitz(conj(flipud(rc(1:M))));

% Calculate Pdx
Pdx = zeros(M,1);
for k=M:K
    Pdx = Pdx + (flipud(x(k-M+1:k)) * conj(d(k)));
end
Pdx = Pdx / (K-M+1);

% Calculate Wiener filter (ww)
wo = inv(Rxx)*Pdx;

%% RLS---------------------
e_squared = zeros(K,length(lambda_array));
squared_dev = zeros(K,length(lambda_array));
for m=1:length(lambda_array)
    P = eye(M);
    y = zeros(K,1);
    w = zeros(M,1);
    lambda = lambda_array(m);
    for n=M+1:K
        xn = flipud(x(n-M+1:n));

        pi = P*xn;
        k = pi / (lambda + ctranspose(xn)*pi);
        y(n) = ctranspose(w)*xn;
        e = d(n) - y(n);
        w = w + k*conj(e);
        P = (P/lambda) - (k*ctranspose(xn)*P)/lambda;

        e_squared(n,m) = conj(e)*e;
        squared_dev(n,m) = ctranspose(w-wo)*(w-wo);
    end
end

%% Plots----------------------------
figure(1);

plot(abs(d));
hold on
plot(abs(y),'color','red');

figure(2);
plot(20*log10(e_squared(:,1)));
hold on
plot(20*log10(e_squared(:,2)),'color','red');
plot(20*log10(e_squared(:,3)),'color','green');

figure(3);
plot(20*log10(squared_dev(:,1)));
hold on;
plot(20*log10(squared_dev(:,2)),'color','red');
plot(20*log10(squared_dev(:,3)),'color','green');

