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

X_new = ctranspose(V)*X;
R_new = (1/N)*X_new*ctranspose(X_new);

lambda_new = eig(R_new)
condition_number_new = cond(R_new)

% Show diag(R_new) = lambda_new
diag(R_new)

% R_new shows strong diagonal, pre=multiplying ctranspose(V) to X yields whitened signal
% Calculate Cov(Xnew) = Cov(V'XXV) = V'RoldV
% Then Rold = VRnewV' means Rnew is diag matrix of eigenvalues
% verify this by comparing diagonal of Rnew with lambda_new

% Plots
figure(1)
eigval_R_new = 20*log10(real(lambda_new));
h0(1) = subplot(2,1,1);
plot(eigval_R_new,'color','blue');
title('Exam 1 Problem 4');
xlabel('Eigenvalue index');
ylabel('Eigenvalue (dB)');
grid on;
legend({'eigval\_R\_new'});
h0(2) = subplot(2,1,2);
plot(20*log10(real(diag(R_new))),'color','red');
grid on;
legend({'diag\_R\_new'});
xlabel('Eigenvalue index');
ylabel('Eigenvalue (dB)');

