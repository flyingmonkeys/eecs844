% Exam 2 Problem 2

clear all;
close all;
hold off;

d = randn(5000,1) + j*randn(5000,1);
h = [0.7+j*0.7 0.2-j*0.5 -0.35+j*0.9 0.2+j*0.5];
x = filter(conj(h),1,d);

M_MAX = 50;
W            = zeros(M_MAX,1);
mse_direct   = zeros(M_MAX+1,1);
mse_analytic = zeros(M_MAX+1,1);
K = length(x);
e = zeros(K,1);

sigma_d_squared = (ctranspose(d)*d)/length(d); % for analytic MSE calc

% Cycle through all M filter lengths
for M=0:M_MAX
    % Calculate filter weights
    if( M == 0 )
        W(1) = 1.0; % no filter
    else
        % Calculate Wiener filter weights
        
        % Calculate R
%        xc = xcorr(x,x,M-1,'unbiased');
%        r = xc(M:end);
%        R = toeplitz(r,conj(r));
        
        % Create X matrix
        X = zeros(M,K-M+1);
        for k=1:K-M+1
            X(:,k) = flipud(x(k:k+M-1));
        end
        
        % Calculate R
        R = (1/(K-M+1))*X*ctranspose(X);

        % Calculate P
        P = zeros(M,1);
        for k=M:length(d)
            P = P + (flipud(x(k-M+1:k)) * conj(d(k)));
        end
        P = P / (length(d)-M+1);
        
        % Calculate W
        W = R\P;
        
        % Calculate MSE (analytic) = Jmin
        sigma_dhat_squared = real(ctranspose(P)*inv(R)*P);
        mse_analytic(M+1) = sigma_d_squared - sigma_dhat_squared;
    end
    
    % Calculate filter output
    y = filter(conj(W),1,x);
    
    % Calculate estimation error
    e = d - y(1:length(d));
    
    % Calculate MSE (direct)
    mse_direct(M+1) = (ctranspose(e)*e)/length(e);
    if( M == 0 )
        mse_analytic(M+1) = mse_direct(M+1); % no R or P to calc analytic error
    end
end

% Calculate Ryy
%Ryy = ctranspose(h)*R*h;
%Pxy = ctranspose(h)*P;
%ww = inv(Ryy)*Pxy;

%z = filter(conj(ww),1,x);

% Plots
figure(1)
x_axis_label = linspace(0,50,51);
mse_curve_direct = 20*log10(mse_direct);
mse_curve_analytic = 20*log10(mse_analytic);
plot(x_axis_label,mse_curve_direct,'color','blue');
hold on;
plot(x_axis_label,mse_curve_analytic,'color','red');
title('MSE vs. Filter Length');
legend({'MSE\_direct','MSE\_analytic'});
grid on;
xlabel('Filter length');
ylabel('MSE (dB)');

figure(2)
h0(1)=subplot(2,1,1);
plot(abs(W)); % magnitude
title('Magnitude of W, M=50');
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
title('Frequency response of W, M=50');
