% Exam 3 Problem 2

clear all;
close all;
hold off;

load p2.mat; % loads R,s1,s2

M = 2000;           % number of iterations
N = length(s1);
w = zeros(N,M);     % weight vector, init to zeros
J = zeros(M,1);     % cost function
mu = 1/(300*N);     % step size

numsamps = 20*N;
theta = linspace(-pi,pi,numsamps);

for idx=1:M
    del_J = 2*R*w(:,idx) + ...
            2*s1*(real(ctranspose(w(:,idx))*s1) + imag(ctranspose(w(:,idx))*s1) - 1) + ...
            2*s2*(real(ctranspose(w(:,idx))*s2) + imag(ctranspose(w(:,idx))*s2) - 0.5);
%            2*(ctranspose(w(:,idx))*s1 - 1)*s1 + ...
%            2*(ctranspose(w(:,idx))*s2 - 0.5)*s2;
%        2*s1*(ctranspose(w(:,idx))*s1);
%            2*( real(s1)*real(ctranspose(w(:,idx))*s1) -     real(s1) + imag(s1)*imag(ctranspose(w(:,idx))*s1) ) + ...
%            2*( real(s2)*real(ctranspose(w(:,idx))*s2) - 0.5*real(s2) + imag(s2)*imag(ctranspose(w(:,idx))*s2) );

    w(:,idx+1) = w(:,idx) - 0.5*mu*del_J;
    J(idx) = ctranspose(w(:,idx+1))*R*w(:,idx+1) + ...
             abs(ctranspose(w(:,idx+1))*s1 - 1)^2 + ...
             abs(ctranspose(w(:,idx+1))*s2 - 0.5)^2;
end

beampattern = zeros(numsamps,1);
for idx=1:numsamps
    % Compute beampattern vs electrical angle
    s_e = transpose(exp(-j*theta(idx)*linspace(0,N-1,N)));
    beampattern(idx) = ctranspose(s_e)*w(:,2000);
end

figure(1);
plot(real(J));

figure(2);
plot(theta,20*log10(abs(beampattern)));
title('Adaptive filter response');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -60 5]);
hold on;
grid on;
