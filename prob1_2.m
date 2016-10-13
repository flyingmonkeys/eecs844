clear all;
close all;
hold off;

load p1.mat; % loads x and d

M_MAX = 50;
W = zeros(M_MAX,1);
mse_direct = zeros(M_MAX+1,1);
mse_analytic = zeros(M_MAX+1,1);
K = length(x);
e = zeros(K,1);

sigma_d_squared = (ctranspose(d)*d)/length(d); % for analytic MSE calc

%M = 5;
for M=0:M_MAX
    % Calculate filter weights
    if( M == 0 )
        W(1) = 1.0; % no filter
    else
        % Calculate Wiener filter weights
        
        % Calculate R
        xc = xcorr(x,x,M-1,'unbiased');
        r = xc(M:end);
        R = toeplitz(r,conj(r));
        
        % Create X matrix
        X = zeros(M,K-M+1);
        for k=1:K-M+1
            X(:,k) = flipud(x(k:k+M-1));
        end
        RR = (1/(K-M+1))*X*ctranspose(X);

        % Calculate P
        P = zeros(M,1);
        for k=M:length(d)
            P = P + (flipud(x(k-M+1:k)) * conj(d(k)));
        end
        P = P / (length(d)-M);
        
        % Calculate W
        W = RR\P;
        
        % Calculate MSE (analytic) = Jmin
        sigma_dhat_squared = real(ctranspose(P)*inv(RR)*P);
        mse_analytic(M+1) = sigma_d_squared - sigma_dhat_squared;
    end
    
    % Calculate filter output
    y = filter(conj(W),1,x);
    yy = conv(W,x);
    
    % Calculate estimation error
    e = d - y(1:length(d));
    
    % Calculate MSE (direct)
    mse_direct(M+1) = (ctranspose(e)*e)/length(e);
    if( M == 0 )
        mse_analytic(M+1) = mse_direct(M+1); % no R or P to calc analytic error
    end
end

% Plots
figure(1)
mse_curve_direct = 20*log10(mse_direct);
mse_curve_analytic = 20*log10(mse_analytic);
plot(mse_curve_direct,'color','blue');
hold on;
plot(mse_curve_analytic,'color','red');

title('Exam 2 Problem 1');
legend({'MSE\_direct','MSE\_analytic'});
grid on;
xlabel('Filter length');
ylabel('MSE (dB)');

figure(2)
h0(1)=subplot(2,1,1);
plot(abs(W)); % magnitude
title('Magnitude of W');
ylabel('Magnitude');
xlabel('Filter tap #');
h0(2)=subplot(2,1,2);
plot(angle(W)); % phase
title('Phase');
ylabel('Phase (radians)');
xlabel('Filter tap #');
linkaxes(h0,'x');


figure(3)
freqz(W);
title('Exam 2 Problem 1 Frequency response of W');
