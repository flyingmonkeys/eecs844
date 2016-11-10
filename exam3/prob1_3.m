% Exam 3 Problem 1

clear all;
close all;
hold off;

load p1.mat; % loads x,d

K = length(x);
M_lp = 29; % Forward linear prediction filter order
M_w  = 30; % Wiener filter order

rc = xcorr(d,M_w-1,'unbiased');
Rdd = toeplitz(conj(flipud(rc(1:M_w))));
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

% Run d through LP filter
y_lp = zeros(10000,1);
for n=M_lp+1+1:10000
    sum = 0;
    for k=1:M_lp+1
        sum = sum + conj(a(k))*d(n-k+1);
    end
    y_lp(n) = sum;
end

% Run d through Wiener filter
y = zeros(10000,1);
for n=M_w+1:10000
    sum = 0;
    for k=1:M_w
        sum = sum + conj(ww(k))*d(n-k+1);
    end
    y(n) = sum;
end

%% Plots--------------
figure(1);
freqz(ww);
title('Frequency/Phase response, Wiener filter');

figure(2);
freqz(a);
title('Frequency/Phase response, forward linear prediction filter');

figure(3);
plot(abs(ww),'color','blue');
grid on;
hold on;
plot(abs(a),'color','red');
title('Filter Magnitude');
legend({'Wiener filter','Forward LP'});

figure(4);
plot(angle(ww),'color','blue');
grid on;
hold on;
plot(angle(a),'color','red');
title('Filter Phase');
legend({'Wiener filter','Forward LP'});

figure(5);
e = y - x;
plot(abs(e));
title('Wiener filter error');

figure(6);
e = y_lp - x;
plot(abs(e));
title('LP filter error');

figure(7);
plot(abs(y));
title('Filter output');
hold on;
plot(abs(y_lp),'color','red');
legend({'Wiener','LP'});