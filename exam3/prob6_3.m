% Exam 3 Problem 6

clear all;
close all;
hold off;

load p6.mat; % loads X

L = length(X);   % number of snapshots
M = 30;          % number of array elements
numsamps = 20*M; % oversampled angular increments

s = zeros(M,1); % steering vector matrix
MVDR = zeros(numsamps,1); % MVDR power spectrum
MVDR_fb = zeros(numsamps,1); % MVDR power spectrum
p = zeros(numsamps,1); % non-adaptive power spectrum

% Compute correlation matrix based on L snapshots
R = (1/L)*X*ctranspose(X);
R_I = R + eye(length(R));
R_inv = inv(R);
R_I_inv = inv(R_I);

%% Compute Power spectrums
theta = linspace(-pi,pi,numsamps);

% Compute non-adaptive power spectrum
for idx=1:numsamps
    % Compute steering matrix (electrical angle)
    s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    for k=1:L-5
        p(idx) = p(idx) + (ctranspose(s)*X(:,k))^2;
    end
    p(idx) = p(idx) / (M*L);
end

% Form a reflection matrix and forward-backward estimate of R
J = zeros(M,M);
for idx=1:M
    J(M-idx+1,M-idx+1) = 1;
end
Rfb = (1/(2*L))*(X*ctranspose(X) + J*conj(X)*transpose(X)*J);
Rfb_I = Rfb + eye(length(Rfb));
Rfb_inv = inv(Rfb);
Rfb_I_inv = inv(Rfb_I);

% Compute MVDR for both R matrices
for idx=1:numsamps
    % Compute steering matrix (electrical angle)
    s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    MVDR(idx) = 1.0 / ( ctranspose(s)*R_inv*s );
    MVDR_fb(idx) = 1.0 / (ctranspose(s)*Rfb_inv*s );
    MVDR_I(idx) = 1.0 / ( ctranspose(s)*R_I_inv*s );
    MVDR_fb_I(idx) = 1.0 / (ctranspose(s)*Rfb_I_inv*s );
end


%% Plots

% MVDR Power Spectrum (electrical angle)
figure(1);
plot(theta,20*log10(abs(p)));
title('MVDR Power Spectrum');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -40 60]);
hold on;
plot(theta,20*log10(abs(MVDR)),'color','red');
plot(theta,20*log10(abs(MVDR_fb)),'color','green');
plot(theta,20*log10(abs(MVDR_I)),'color','blue');
plot(theta,20*log10(abs(MVDR_fb_I)),'color','yellow');
legend({'Non-adaptive','MVDR R','MVDR Rfb','MVDR R loaded','MVDR Rfb loaded'});
grid on;






