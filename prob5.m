clear all;

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


