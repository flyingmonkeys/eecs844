% Exam 2 Problem 7

clear all;
close all;
hold off;

load p7.mat; % loads x

K = length(x);   % number of samples
M = 10;          % filter length
numsamps = 20*M; % oversampled angular increments
L = [10 20 K-M+1]; % snapshot lengths

s = zeros(M,1); % steering vector matrix
MVDR = zeros(numsamps,length(L)); % MVDR power spectrums
MVDR_loaded = zeros(numsamps,length(L)); % loaded MVDR power spectrum

% Create X matrix (snapshots)
X = zeros(M,K-M+1);
for k=1:K-M+1
    X(:,k) = flipud(x(k:k+M-1));
end

% Cycle through correlation lengths of 10,20,and all
theta = linspace(-pi,pi,numsamps);
for L_idx=1:3
    % Compute correlation matrix based on L snapshots
    R = (1/L(L_idx))*X(:,1:L(L_idx))*ctranspose(X(:,1:L(L_idx)));
    R_loaded = R + eye(M);
    R_inv = inv(R);
    R_loaded_inv = inv(R_loaded);

    % Compute steering matrix
    for idx=1:numsamps
        % Compute steering matrix (electrical angle)
        s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
        MVDR(idx,L_idx) = 1.0 / ( ctranspose(s)*R_inv*s );
        MVDR_loaded(idx,L_idx) = 1.0 / ( ctranspose(s)*R_loaded_inv*s );
    end
end

% MVDR Power Spectrum (electrical angle)
figure(1);
plot(theta,20*log10(abs(MVDR(:,1))));
title('MVDR Power Spectrum');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -80 20]);
hold on;
plot(theta,20*log10(abs(MVDR(:,2))),'color','red');
plot(theta,20*log10(abs(MVDR(:,3))),'color','green');
legend({'10 snapshots','20 snapshots','All snapshots'});
grid on;

% MVDR Power Spectrum, loaded correlation matrix (electrical angle)
figure(2);
plot(theta,20*log10(abs(MVDR_loaded(:,1))));
title('MVDR Power Spectrum with loaded R matrix');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -30 20]);
hold on;
plot(theta,20*log10(abs(MVDR_loaded(:,2))),'color','red');
plot(theta,20*log10(abs(MVDR_loaded(:,3))),'color','green');
legend({'10 snapshots','20 snapshots','All snapshots'});
grid on;


