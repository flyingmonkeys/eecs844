% Exam 2 Problem 5

clear all;
close all;
hold off;

load p3.mat; % loads X

Lx = length(X);   % number of snapshots
L = 2;           % number of constraints
M = 10;          % number of array elements
numsamps = 20*M; % oversampled angular increments

null_angle = -0.6788; % null angle, radians

MVDR = zeros(numsamps,1); % MVDR power spectrum
p = zeros(numsamps,1); % non-adaptive power spectrum

% Compute correlation matrix based on L snapshots
R = (1/Lx)*X*ctranspose(X);
R_inv = inv(R);

% Form steering vectors
s0 = transpose(exp(-j*0.0*linspace(0,M-1,M))); % boresight
s1 = transpose(exp(-j*null_angle*linspace(0,M-1,M))); % null
g = [1.0 0.0]';

% Form C, Ca
C = [s0 s1];
[U,S,V] = svd(C);
Ca = U(:,L+1:M);

v = inv(ctranspose(C)*C) * g;
wq = C*v;

% Compute wa(optimal)
wa = inv(ctranspose(Ca)*R*Ca)*ctranspose(Ca)*R*wq;

w = wq - Ca*wa;

%% Compute Power spectrums
theta = linspace(-pi,pi,numsamps);


% Compute non-adaptive beamformer response
for idx=1:numsamps
    % Compute steering matrix (electrical angle)
    s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    pp(idx) = ctranspose(s)*wq;
    p(idx) = ctranspose(s)*w;
end


%% Compute MVDR beamformers
MVDR_response = zeros(numsamps,1);
s = transpose(exp(-j*0.0*linspace(0,M-1,M)));
w_mvdr = (R_inv*s) / (ctranspose(s)*R_inv*s);
% Compute beampatterns
for idx=1:numsamps
    s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    MVDR_response(idx) = ctranspose(s)*w_mvdr;
end

%% Plots

% MVDR Power Spectrum (electrical angle)
figure(1);
plot(theta,20*log10(abs(p)));
title('GSC Beampattern');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -50 30]);
hold on;
plot(theta,20*log10(abs(pp)),'color','red');
plot(theta,20*log10(abs(MVDR_response)),'color','green');
legend({'Adaptive+non-adaptive','Non-adaptive','MVDR 0 degrees'});
grid on;


