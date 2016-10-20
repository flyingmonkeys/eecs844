% Exam 2 Problem 6

clear all;
close all;
hold off;

load p3.mat; % loads X

Lx = length(X);  % number of snapshots
L = 6;           % number of constraints
M = 10;          % number of array elements
numsamps = 20*M; % oversampled angular increments

MVDR = zeros(numsamps,1); % MVDR power spectrum
p = zeros(numsamps,1);    % GSC power spectrum
pp =  zeros(numsamps,1);  % non-adaptive power spectrum

% Compute correlation matrix based on Lx snapshots
R = (1/Lx)*X*ctranspose(X);
R_inv = inv(R);

% Form steering vectors
s0 = transpose(exp(-j*0.0*linspace(0,M-1,M))); % boresight
null_angles = (pi/180)*linspace(-50,-15,M); % radians
Rc = zeros(M,M);
for k=1:M
    s = transpose(exp(-j*null_angles(k)*linspace(0,M-1,M)));
    Rc = Rc + s*ctranspose(s);
end
[Q,lambda] = eig(Rc);

% by inspection, there are two principle eigenvalues, but take up to L
% eigenvalues (L constraints)
C = Q(:,M-L+1:M);
C(:,1) = s0;      % boresight constraint

g = zeros(L,1);
g(1) = 1.0; % unity gain constraint

% Form Ca
[U,S,V] = svd(C);
Ca = U(:,L+1:M);

% Compute wq
v = inv(ctranspose(C)*C) * g;
wq = C*v;

% Compute wa(optimal)
wa = inv(ctranspose(Ca)*R*Ca)*ctranspose(Ca)*R*wq;
w = wq - Ca*wa;

%% Compute Power spectrums
theta = linspace(-pi,pi,numsamps);
theta_deg = theta*180/pi;

% Compute non-adaptive beamformer response
for idx=1:numsamps
    % Compute steering matrix (electrical angle)
    s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    pp(idx) = ctranspose(s)*wq; % adaptive power spectrum
    p(idx) = ctranspose(s)*w;   % GSC power spectrum
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
plot(theta_deg,20*log10(abs(p)));
title('GSC Beampattern, 6 contraints');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-180 180 -60 40]);
hold on;
plot(theta_deg,20*log10(abs(pp)),'color','red');
plot(theta_deg,20*log10(abs(MVDR_response)),'color','green');
legend({'Adaptive+non-adaptive (w)','Non-adaptive (wq)','MVDR 0 degrees'});
grid on;


