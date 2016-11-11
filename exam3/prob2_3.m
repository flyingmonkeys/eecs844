% Exam 3 Problem 2

clear all;
close all;
hold off;

load p2.mat; % loads R,s1,s2

M = 2000;           % number of iterations
N = length(s1);
mu = [1/(10*N) 1/(50*N) 1/(300*N)];     % step size array
J = zeros(M,length(mu));   % cost function for different step sizes

numsamps = 20*N;
beampattern = zeros(numsamps,length(mu)); % Frequency responses for converged filters
theta = linspace(-180,180,numsamps);

%% Run steepest-descent algorithm over all step sizes
for step_size_idx=1:length(mu)
    % Initialize
    w = zeros(N,M);     % weight vector, init to zeros
    
    % Iterate
    for idx=1:M
        del_J = 2*R*w(:,idx) + ...
                2*s1*(real(ctranspose(w(:,idx))*s1) + imag(ctranspose(w(:,idx))*s1) - 1) + ...
                2*s2*(real(ctranspose(w(:,idx))*s2) + imag(ctranspose(w(:,idx))*s2) - 0.5);

        w(:,idx+1) = w(:,idx) - 0.5*mu(step_size_idx)*del_J;
        J(idx,step_size_idx) = ctranspose(w(:,idx+1))*R*w(:,idx+1) + ...
                 abs(ctranspose(w(:,idx+1))*s1 - 1)^2 + ...
                 abs(ctranspose(w(:,idx+1))*s2 - 0.5)^2;
    end

    for idx=1:numsamps
        % Compute beampattern vs electrical angle
        s_e = transpose(exp(-j*theta(idx)*(pi/180)*linspace(0,N-1,N)));
        beampattern(idx,step_size_idx) = ctranspose(s_e)*w(:,2000);
    end
end
%% Plots------------------------
figure(1);
plot(real(J(:,1)));
hold on;
plot(real(J(:,2)),'color','red');
plot(real(J(:,3)),'color','green');
title('Convergence curves');
legend({'mu=1/10N','mu=1/50N','mu=1/300N'});
ylabel('Cost function J');
xlabel('Iteration number');
grid on;

figure(2);
plot(theta,20*log10(abs(beampattern(:,1))));
title('Adaptive filter response');
xlabel('Electrical Theta (deg)');
ylabel('Magnitude (dB)');
axis([-180 180 -60 5]);
hold on;
plot(theta,20*log10(abs(beampattern(:,2))));
plot(theta,20*log10(abs(beampattern(:,3))));
grid on;
legend({'mu=1/10N','mu=1/50N','mu=1/300N'});
