% Problem 5, Exam 4

clear all;
close all;
hold off;

load p4.mat; % loads x
K = length(x);

p_vals = [4 6 8];

% Cycle through all model orders
for idx=1:length(p_vals)
    p = p_vals(idx);
    M = p+1;
    
    % Create X matrix
    X = zeros(M,K-M+1);
    for k=1:K-M+1
        X(:,k) = flipud(x(k:k+M-1));
    end

    % Calculate R
    R = (1/(K-M+1))*X*ctranspose(X);
    rxx = R(2:M,1);
    Rxx = R(1:p,1:p);
    
    % Calculate AR coefficients
    a = -inv(Rxx)*rxx;
    
    % Evaluate using freqz
    ar = zeros(M,1);
    ar(1) = 1.0;
    ar(2:M) = conj(a);
    
    % Compute AR using MATLAB for comparison
    a = aryule(x,p);
    
    figure(idx);
    freqz(1,ar,200);
end


