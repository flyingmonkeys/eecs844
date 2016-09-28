clear all;
close all;
hold off;

load p2.mat; % loads x

M = 8;
K = length(x);
N = K-M+1;
X = zeros(M,K);
R = zeros(M,M);

% Create X matrix
for k=1:N
	X(:,k) = flipud(x(k:k+M-1));
end

% Calculate autocorrelation matrix R
R = (1/N)*X*ctranspose(X);

% Get eigenvalues and eigenvectors of R
[V, lambda] = eig(R);
lambda = diag(lambda);
condition_number = cond(R) % same as max(lambda)/min(lambda)
[V_invR, lambda_invR] = eig(inv(R));
lambda_invR = diag(lambda_invR); % these should yield 1/diag(lambda)
condition_number_invR = cond(inv(R));

% eigenvalues of inv(R) are 1/eigenvalues of R
% since R is square PSDH, it can be decomposed via eigen decomposition
% since R is Hermitian, lambdas are real-valued

% Plots
figure(1)
eigval_R = 20*log10(real(lambda));
eigval_Rinv = flipud(20*log10(real(lambda_invR)));
plot(eigval_R,'color','blue');
hold on;
plot(eigval_Rinv,'color','red','linestyle',':');
title('Exam 1 Problem 2');
legend({'eigval\_R','eigval\_Rinv'});
grid on;
xlabel('Eigenvalue index');
ylabel('Eigenvalue (dB)');