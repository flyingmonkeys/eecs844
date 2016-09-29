clear all;
close all;
hold off;

load p5.mat; % loads X (12x100 array)

M = 12;
N = length(X);

R = (1/N)*X*ctranspose(X);

[V, lambda] = eig(R);
lambda = diag(lambda)
c = cond(R)

[V_invR, lambda_invR] = eig(inv(R));
lambda_invR = diag(lambda_invR)
c = cond(inv(R))

% High condition number means high eigenvalue spread, or
% equivalently a rank-deficient spatial correlation matrix.
% Since only a few eigenvectors are dominant, this means
% that signals are dominant in a few directions, and
% the array is highly correlated

% Condition numbers are the same, for the same reasons as
% cited in problem 2

% Plots
figure(1)
eigval_R = 20*log10(real(lambda));
eigval_Rinv = flipud(20*log10(real(lambda_invR)));
plot(eigval_R,'color','blue');
hold on;
plot(eigval_Rinv,'color','red','linestyle',':');
title('Exam 1 Problem 5');
legend({'eigval\_R','eigval\_Rinv'});
grid on;
xlabel('Eigenvalue index');
ylabel('Eigenvalue (dB)');