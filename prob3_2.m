% Exam 2 Problem 3

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

%% Compute Power spectrums
theta = linspace(-pi,pi,numsamps);

% Compute non-adaptive power spectrum
for idx=1:numsamps
    % Compute steering matrix (electrical angle)
    s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    for k=1:L
        p(idx) = p(idx) + (ctranspose(s)*X(:,k))^2;
    end
    p(idx) = p(idx) / (M*L);
end

% Compute MVDR
for idx=1:numsamps
    % Compute steering matrix (electrical angle)
    s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    MVDR(idx) = 1.0 / ( ctranspose(s)*R_inv*s );
end

%% Compute MVDR beamformers
peak_angles = [-0.5841 -0.3947 0.742 0.0]; % peaks in radians
w = zeros(M,length(peak_angles)); % weights for each of the 3 observed peaks
MVDR_response = zeros(numsamps,length(peak_angles));

for k=1:length(peak_angles)
    s = transpose(exp(-j*peak_angles(k)*linspace(0,M-1,M)));
    w(:,k) = (R_inv*s) / (ctranspose(s)*R_inv*s);
    % Compute beampatterns
    for idx=1:numsamps
        s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
        MVDR_response(idx,k) = ctranspose(s)*w(:,k);
    end
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
legend({'Non-adaptive','MVDR'});
grid on;

% MVDR beamformers
figure(2);
h0(1)=subplot(3,1,1);
plot(theta,20*log10(abs(MVDR_response(:,1))));
grid on;
title('MVDR Beampattern, -0.5841 rad');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -30 30]);
h0(2)=subplot(3,1,2);
plot(theta,20*log10(abs(MVDR_response(:,2))),'color','red');
grid on;
title('MVDR Beampattern, -0.3947 rad');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -30 30]);
h0(2)=subplot(3,1,3);
plot(theta,20*log10(abs(MVDR_response(:,3))),'color','green');
grid on;
title('MVDR Beampattern, 0.742 rad');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -30 30]);
linkaxes(h0,'x');


% Boresight constraint
figure(3);
plot(theta,20*log10(abs(MVDR_response(:,4))));
title('MVDR Beampattern, Boresight at unity');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -50 10]);
legend({'0.0 degrees'});
grid on;



