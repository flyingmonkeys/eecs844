clear all;
close all;
hold off;

load p5.mat;

M = 12;
K = length(X);
N = K-M+1;

R = (1/N)*X*ctranspose(X);

[V, lambda] = eig(R);
diag(lambda)
c = cond(R)

[V_invR, lambda_invR] = eig(inv(R));
diag(lambda_invR)
c = cond(inv(R))

% High condition number means high eigenvalue spread, or
% equivalently a rank-deficient spatial correlation matrix.
% SInce only a few eigenvectors are dominant, this means
% that signals are dominant in a few directions, and
% the array is highly correlated

% Condition numbers are the same, for the same reasons asctime
% cited in problem 2
