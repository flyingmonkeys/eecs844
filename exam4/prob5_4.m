% Exam 4 Problem 5

clear all;
close all;
hold off;

load p4.mat; % loads x

K = length(x);
p_vals = [4 6 8];
L = 200; % number of PSD bins
w = linspace(0,2*pi,L);

figure(1);

% Cycle through all model orders
for idx=1:length(p_vals)
    p = p_vals(idx);
    M = p+1;
    N = K-M+1;
    
    % Create X snapshot matrix
    X = zeros(M,N);
    for k=1:N
        X(:,k) = flipud(x(k:k+M-1));
    end

    % Calculate autocorrelation matrix R
    R = (1/N)*X*ctranspose(X);
    rxx = R(2:M,1);
    Rxx = R(1:p,1:p);
    
    % Calculate AR coefficients
    a = -inv(Rxx)*rxx;
    
    % Evaluate using freqz
    ar = zeros(M,1);
    ar(1) = 1.0;
    ar(2:M) = conj(a);
    
    % Plot frequency response using freqz
    [H,W] = freqz(1,ar,L,'whole');
    plot(w,20*log10(abs(H)));
    hold on;
end

title('Yule-Walker Power Spectral Density');
xlabel('Normalized Frequency (rad/sample)');
ylabel('Magnitude (dB/rad/sample)');
legend({'4','6','8'});
grid on;
