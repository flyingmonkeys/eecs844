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
% Calculate Cov(Xnew) = Cov(V'XXV) = V'RoldV
% Then Rold = VRnewV' means Rnew is diag matrix of eigenvalues
% verify this by comparing diagonal of Rnew with lambda_new
% https://theclevermachine.wordpress.com/tag/covariance-matrix/

% Plots
figure(1)
eigval_R = 20*log10(real(diag(lambda)));
eigval_Rinv = flipud(20*log10(real(diag(lambda_invR))));
plot(eigval_R,'color','blue');
hold on;
plot(eigval_Rinv,'color','red');
title('Exam 1 Problem 2');
legend({'eigval_R','eigval_Rinv'});

figure(2)
eigval_R_loaded = 20*log10(real(lambda_loaded));
eigval_Rinv_loaded = 20*log10(real(lambda_invR_loaded));
plot(eigval_R_loaded,'color','blue');
hold on;
plot(eigval_Rinv_loaded,'color','red');
title('Exam 1 Problem 3');
legend({'eigval_R','eigval_Rinv'});
