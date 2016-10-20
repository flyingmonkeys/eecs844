% Exam 2 Problem 4

clear all;
close all;
hold off;

load p3.mat; % loads X

L = length(X);   % number of snapshots
M = 10;          % number of array elements
s = zeros(M,1); % steering vector matrix

% Compute correlation matrix based on L snapshots
R = (1/L)*X*ctranspose(X);
R_inv = inv(R);

%% Compute MVDR beamformers
s = transpose(exp(-j*0.0*linspace(0,M-1,M))); % boresight direction (0.0 deg)
w = (R_inv*s) / (ctranspose(s)*R_inv*s);

y_mvdr = ctranspose(w)*X; % MVDR
y_na = (1/M)*(ctranspose(s)*X); % non-adaptive

P_mvdr = sqrt((1/length(y_mvdr))*y_mvdr*ctranspose(y_mvdr))

%% Plots

% MVDR Power Spectrum (electrical angle)
figure(1);
plot(abs(X(1,:)));
title('Exam 2 Problem 4');
xlabel('Time');
ylabel('Magnitude');
hold on;
plot(abs(y_na(1,:)),'color','green');
plot(abs(y_mvdr(1,:)),'color','red');
legend({'Raw','Non-adaptive','Boresight MVDR'});
grid on;
