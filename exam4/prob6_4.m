% Problem 7, Exam 4

clear all;

load p6.mat; % loads X1..X4

for a=1:4

switch a
	case 1
		A = X1;
	case 2
		A = X2;
	case 3
		A = X3;
	case 4
		A = X4;
endswitch

[M N] = size(A); % N = number of time samples
p = M; % number of sensors
n = N; % number of time samples

% Calculate autocorrelation matrix R
R = (1/N)*A*ctranspose(A);

% since R is square PSDH, it can be decomposed via eigen decomposition
% since R is Hermitian, lambdas are real-valued
% Get eigenvalues and eigenvectors of R
[V, lambda] = eig(R);
lambda = diag(lambda);
l = flipud(sort(lambda));

BIC = zeros(p,1);

for q=0:p-1
	sum1 = 0; % first summation term
	sum2 = 0; % second summation term
	for i=q+1:p
		sum1 = sum1 + log(l(i));
		sum2 = sum2 + l(i)/(p-q);
	end
	log_lq = n * (sum1 - (p-q)*log(sum2));
	BIC(q+1) = -2*log_lq + (q*((2*p)-q)+1)*log(n);
end

[min_BIC, num_signals] = min(BIC)

% Compute MVDR
R_inv = inv(R);
numsamps = 20*M;
for idx=1:numsamps
    % Compute steering matrix (electrical angle)
    s = transpose(exp(-j*theta(idx)*linspace(0,M-1,M)));
    MVDR(idx) = 1.0 / ( ctranspose(s)*R_inv*s );
end

% MVDR Power Spectrum (electrical angle)
figure(a);
plot(theta,20*log10(abs(MVDR)));
title('MVDR Power Spectrum');
xlabel('Electrical Theta (rad)');
ylabel('Magnitude (dB)');
axis([-pi pi -40 60]);

end


