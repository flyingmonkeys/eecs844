% Exam 2 Problem 4

clear all;
close all;
hold off;

load p3.mat; % loads X

L = length(X);   % number of snapshots
M = 10;          % number of array elements
numsamps = 20*M; % oversampled angular increments

s = zeros(M,1); % steering vector matrix
MVDR = zeros(numsamps,1); % MVDR power spectrum
p = zeros(numsamps,1); % non-adaptive power spectrum

% Compute correlation matrix based on L snapshots
R = (1/L)*X*ctranspose(X);
R_inv = inv(R);


%% Compute MVDR beamformers
w = zeros(M,1);

s = transpose(exp(-j*0.0*linspace(0,M-1,M)));

w = (R_inv*s) / (ctranspose(s)*R_inv*s);
y_0 = ctranspose(w)*X;

y_na = (1/M)*(ctranspose(s)*X);

%% Plots

% MVDR Power Spectrum (electrical angle)
figure(1);
plot(abs(X(1,:)));
title('Exam 2 Problem 4');
xlabel('Time');
ylabel('Magnitude');
hold on;
plot(abs(y_na(1,:)),'color','green');
plot(abs(y_0(1,:)),'color','red');
legend({'Raw','Non-adaptive','Boresight MVDR'});
grid on;
