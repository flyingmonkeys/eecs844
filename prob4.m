clear all;
hold off;

load p2.mat;

M = 8;
K = length(x);
N = K-M+1;
X = zeros(M,K);
R = zeros(M,M);

for k=1:N
	X(:,k) = flipud(x(k:k+M-1));
end

R = (1/N)*X*ctranspose(X);

[V, lambda] = eig(R);
diag(lambda)
c = cond(R) % same as max(lambda)/min(lambda)
[V_invR, lambda_invR] = eig(inv(R));
diag(lambda_invR) % these should yield 1/diag(lambda)
c_invR = cond(inv(R))

% eigenvalues of inv(R) are 1/eigenvalues of R
% since R is square PSDH, it can be decomposed via eigen decomposition
% since R is Hermitian, lambdas are real-valued

X_new = ctranspose(V)*X;
R_new = (1/N)*X_new*ctranspose(X_new);

lambda_new = eig(R_new)
c_new = cond(R_new)

% R_new shows strong diagonal, pre=multiplying ctranspose(V) to X yields whitened signal
% Calculate Cov(Xnew) = Cov(V'XXV) = V'RoldV
% Then Rold = VRnewV' means Rnew is diag matrix of eigenvalues
% verify this by comparing diagonal of Rnew with lambda_new
% https://theclevermachine.wordpress.com/tag/covariance-matrix/

% Plots
figure(1)
eigval_R_new = 20*log10(real(lambda_new));
plot(eigval_R_new,'color','blue');
title('Exam 1 Problem 4');

