% Exam 3 Problem 1

clear all;
close all;
hold off;

load p1.mat; % loads x,d

K = length(x);
M_lp = 29; % Forward linear prediction filter order
M_w  = 30; % Wiener filter order

%rc = xcorr(d,M_lp);
%rc1 = xcorr(d,M_lp,'coeff');
%r = flipud(rc1(1:M_lp));
%r1 = flipud(rc1(2:M_lp+1));
%Rxx = toeplitz(r1);

% Create D matrix
D = zeros(M_w,K-M_w+1);
for k=1:K-M_w+1
    D(:,k) = flipud(d(k:k+M_w-1));
end

% Calculate Rdd
Rdd = (1/(K-M_w+1))*D*ctranspose(D);

%r(1) = var(d);
%r = zeros(M_lp,1);
r = Rdd(2:M_lp+1,1);

% Calculate LP filter (w_lp)
w_lp = inv(Rdd(1:M_lp,1:M_lp))*r;
a = zeros(M_lp+1,1);
a(1) = 1;
a(2:M_lp+1) = -w_lp;

% Calculate P
P = zeros(M_w,1);
for k=M_w:length(x)
    P = P + (flipud(d(k-M_w+1:k)) * conj(x(k)));
end
P = P / (length(x)-M_w+1);

% Calculate Wiener filter (ww)
ww = inv(Rdd)*P;

%a = abs(conv(w_lp,d));

% Run d through LP filter
y_lp = zeros(10000,1);
for n=M_lp+1+1:10000
    sum = 0;
    for k=1:M_lp+1
%        sum = sum + conj(w_lp(k))*d(n-k);
        sum = sum + conj(a(k))*d(n-k+1);
    end
    y_lp(n) = sum;
end
%y_lp = conv(conj(w_lp),d);

% Run d through Wiener filter
y = zeros(10000,1);
for n=M_w+1:10000
    sum = 0;
    for k=1:M_w
        sum = sum + conj(ww(k))*d(n-k+1);
    end
    y(n) = sum;
end


figure(1);
plot(abs(y));
title('Wiener filter output, M=30');

figure(2);
plot(abs(y_lp));
title('Forward linear predictor output, M=29');

figure(3);
freqz(ww);
title('Frequency/Phase response, Wiener filter');

figure(4);
freqz(a);
title('Frequency/Phase response, forward linear prediction filter');

%for n=1:30
%    figure(n+10);
%    e = y(n:10000) - x(1:10000-n+1);
%    plot(abs(e));
%end
figure(5);
e = y - x;
plot(abs(e));
title('Wiener filter error');

figure(6);
e = y_lp - x;
plot(abs(e));
title('LP filter error');

figure(2);