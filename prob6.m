clear all;
close all;
hold off;

load p6.mat

M = length(w_adap);
s_e = zeros(M,1);
s_s = zeros(M,1);
numsamps = 20*M;
beampattern_w_adap_e = zeros(numsamps,1);
beampattern_w_non_adap_e = zeros(numsamps,1);
beampattern_w_adap_s = zeros(numsamps,1);
beampattern_w_non_adap_s = zeros(numsamps,1);

theta = linspace(-pi,pi,numsamps);
phi   = linspace(-pi/2,pi/2,numsamps);

for idx=1:numsamps
    % Compute beampattern vs electrical angle
    s_e = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    beampattern_w_adap_e(idx) = ctranspose(s_e)*w_adap;
    beampattern_w_non_adap_e(idx) = ctranspose(s_e)*w_non_adap;

    % Compute beampattern vs spatial angle
    s_s = transpose(exp(-j*pi*sin(phi(idx))*linspace(0,M-1,M)));
    beampattern_w_adap_s(idx) = ctranspose(s_s)*w_adap;
    beampattern_w_non_adap_s(idx) = ctranspose(s_s)*w_non_adap;
end

% Plots

% Beampattern (electrical angle)
figure(1);
plot(theta,20*log10(abs(beampattern_w_adap_e)));
title('Adaptive filter response');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -60 5]);
hold on;
plot(theta,20*log10(abs(beampattern_w_non_adap_e)),'color','red','linestyle',':');
legend({'Adaptive','Non-adaptive'});
grid on;

% Beampattern (spatial angle)
figure(2);
plot(theta,20*log10(abs(beampattern_w_adap_s)));
title('Adaptive filter response');
xlabel('Spatial Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -60 5]);
hold on;
plot(theta,20*log10(abs(beampattern_w_non_adap_s)),'color','red','linestyle',':');
legend({'Adaptive','Non-adaptive'});
grid on;