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
[H1,W] = freqz(ww);
[H2,W] = freqz(a);
W = linspace(0,1,length(H1));
h0(1)=subplot(2,1,1);
plot(W,20*log10(abs(H1)));
grid on;
hold on;
plot(W,20*log10(abs(H2)),'color','red');
title('Frequency/Phase Response');
ylabel('Magnitude (dB)');
xlabel('Normalized Frequency (x pi rad/sample)');
legend({'Wiener filter','Forward LP filter'});
h0(2)=subplot(2,1,2);
plot(W,angle(H1)*180/pi);
grid on;
hold on;
plot(W,angle(H2)*180/pi,'color','red');
ylabel('Phase (deg)');
xlabel('Normalized Frequency (x pi rad/sample)');
linkaxes(h0,'x');
legend({'Wiener filter','Forward LP filter'});

figure(2);
plot(abs(ww),'color','blue');
grid on;
hold on;
plot(abs(a),'color','red');
title('Filter Magnitude');
ylabel('Magnitude');
xlabel('Tap number');
legend({'Wiener filter','Forward LP'});

figure(3);
plot(angle(ww)*180/pi,'color','blue');
grid on;
hold on;
plot(angle(a)*180/pi,'color','red');
title('Filter Phase');
ylabel('Phase (deg)');
xlabel('Tap number');
legend({'Wiener filter','Forward LP'});













%{
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
%}