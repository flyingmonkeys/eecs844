% Exam 3 Problem 4

clear all;
close all;
hold off;

load p4.mat; % loads x

M = length(x);
N = 2*M;
m = 2*M;

%% Compute normalized matched filter (h_nmf)
h_nmf = flipud(conj(x))/(ctranspose(x)*x);

figure(1);
plot(20*log10(abs(conv(x,h_nmf))));
grid on;
title('Normalized Matched Filter');

%% Compute LS mismatched filter (h_mmf)

% Form A
a = zeros(3*M,1);
a(1:M) = x;
A = toeplitz(a,zeros(N+1,1));

% Form elementary vector e
e = zeros(3*M,1);
e(m) = 1;

h_mmf = inv(ctranspose(A)*A)*ctranspose(A)*e;
h_nmmf = h_mmf/(sqrtm(ctranspose(h_mmf)*h_mmf)*sqrtm(ctranspose(x)*x));
mismatch_loss = -20*log10( max(abs(conv(h_nmmf,x)))/max(abs(conv(h_nmf,x))) );
fprintf('Mismatch loss, LS mismatched filter: %.2f dB\n',mismatch_loss);

figure(2);
plot(20*log10(abs(conv(x,h_nmmf))),'color','black');
grid on;
title('LS Mismatched Filter');
hold on;


%% Compute LS mismatched filter with modified A
A_mod = A;
A_mod(m-1,:) = zeros(1,N+1);
A_mod(m+1,:) = zeros(1,N+1);

h_mmf_mod = inv(ctranspose(A_mod)*A_mod)*ctranspose(A_mod)*e;
h_nmmf_mod = h_mmf_mod/(sqrtm(ctranspose(h_mmf_mod)*h_mmf_mod)*sqrtm(ctranspose(x)*x));
mismatch_loss = -20*log10( max(abs(conv(h_nmmf_mod,x)))/max(abs(conv(h_nmf,x))) );
fprintf('Mismatch loss, LS mismatched filter with modified A: %.2f dB\n',mismatch_loss);

%figure(3);
plot(20*log10(abs(conv(x,h_nmmf_mod))));
%grid on;
%title('LS Mismatched Filter, modified A');

%% Compute LS mismatched filter with diagonal loading
h_mmf_diag = inv(ctranspose(A)*A + 0.02*max(eig(ctranspose(A)*A))*eye(N+1,N+1))*ctranspose(A)*e;
h_nmmf_diag = h_mmf_mod/(sqrtm(ctranspose(h_mmf_diag)*h_mmf_diag)*sqrtm(ctranspose(x)*x));
mismatch_loss = -20*log10( max(abs(conv(h_nmmf_diag,x)))/max(abs(conv(h_nmf,x))) );
fprintf('Mismatch loss, LS mismatched filter with diagonal loading: %.2f dB\n',mismatch_loss);

%figure(4);
hold on;
plot(20*log10(abs(conv(x,h_nmmf_diag))),'color','red');
%grid on;
%title('LS Mismatched Filter, diagonal loading');

%% Compute LS mismatched filter with modified A and diagonal loading
h_mmf_mod_diag = inv(ctranspose(A_mod)*A_mod + 0.02*max(eig(ctranspose(A_mod)*A_mod))*eye(N+1,N+1))*ctranspose(A_mod)*e;
h_nmmf_mod_diag = h_mmf_mod_diag/(sqrtm(ctranspose(h_mmf_mod_diag)*h_mmf_mod_diag)*sqrtm(ctranspose(x)*x));
mismatch_loss = -20*log10( max(abs(conv(h_nmmf_mod_diag,x)))/max(abs(conv(h_nmf,x))) );
fprintf('Mismatch loss, LS mismatched filter with modified A and diagonal loading: %.2f dB\n',mismatch_loss);

figure(5);
plot(20*log10(abs(conv(x,h_nmmf_mod_diag))));
grid on;
title('LS Mismatched Filter, modified A, diagonal loading');




