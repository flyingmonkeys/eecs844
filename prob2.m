clear all;

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
c = cond(R)
[V_invR, lambda_invR] = eig(inv(R));
diag(lambda_invR)
c_invR = cond(inv(R))

% eigenvalues of inv(R) are 1/eigenvalues of R
% since R is square PSDH, it can be decomposed via eigen decomposition
% since R is Hermitian, lambdas are real-valued

noise_power = 1;

R_loaded = R + noise_power*eye(M);

lambda_loaded = eig(R_loaded)
c_loaded = cond(R_loaded)
lambda_invR_loaded = eig(inv(R_loaded))

X_new = ctranspose(V)*X;
R_new = (1/N)*X_new*ctranspose(X_new);

lambda_new = eig(R_new)
c_new = cond(R_new)

% R_new shows strong diagonal, pre=multiplying ctranspose(V) to X yields whitened signal
% Calculate Cov(Xnew) = Cov(V’XX’V) = V’RoldV
% Then Rold = VRnewV’ means Rnew is diag matrix of eigenvalues
% verify this by comparing diagonal of Rnew with lambda_new
% https://theclevermachine.wordpress.com/tag/covariance-matrix/


