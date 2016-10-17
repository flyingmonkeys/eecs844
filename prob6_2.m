% Exam 2 Problem 6

clear all;
close all;
hold off;

load p3.mat; % loads X

Lx = length(X);   % number of snapshots
L = 5;           % number of constraints
M = 10;          % number of array elements
numsamps = 20*M; % oversampled angular increments

MVDR = zeros(numsamps,1); % MVDR power spectrum
p = zeros(numsamps,1); % non-adaptive power spectrum

% Compute correlation matrix based on L snapshots
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

% by inspection, there are two principle eigenvalues, so take the two
% eigenvectors of Rc corresponding to those two eigenvalues
C = Q(:,M-L+1:M);
C(:,1) = s0;

g = zeros(L,1);
g(1) = 1.0;

% Form Ca
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


