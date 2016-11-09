% Exam 3 Problem 1

clear all;
close all;
hold off;

load p1.mat; % loads x,d

K = length(x);
M = 29;
rc = xcorr(d,M);
rc1 = xcorr(d,M,'coeff');
r = flipud(rc1(1:M));
r1 = flipud(rc1(2:M+1));
R = toeplitz(r1);

% Create X matrix
X = zeros(M,K-M+1);
for k=1:K-M+1
    X(:,k) = flipud(d(k:k+M-1));
end

% Calculate R
R = (1/(K-M+1))*X*ctranspose(X);

r(1:M-1) = conj(R(2:M,1));
r(M) = 0;

% Calculate LP filter
wf = inv(R)*r;

% Calculate P
P = zeros(M,1);
for k=M:length(x)
    P = P + (flipud(d(k-M+1:k)) * conj(x(k)));
end
P = P / (length(x)-M+1);

% Calculate Wiener filter
ww = inv(R)*P;

a = abs(conv(wf,d));

f = zeros(10000,1);
for n=M+1:10000
    sum = 0;
    for k=1:M
        sum = sum + conj(wf(k))*d(n-k);
    end
    f(n) = d(n) - sum;
end
f = conv(conj(wf),d);
fw = conv(conj(ww),d);

plot(abs(xcorr(f,30)));