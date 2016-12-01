% Exam 4 Problem 3

clear all;
close all;
hold off;

load p1.mat; % loads x,d

M = 60; % Wiener filter order
K = length(x);
N = K-M+1;

lambda_array = [0.99 0.998 0.9995];

%% Wiener filter---------------------

% Create X snapshot matrix
X = zeros(M,K);
for k=1:N
    X(:,k) = flipud(x(k:k+M-1));
end

% Calculate autocorrelation matrix R
Rxx = (1/N)*X*ctranspose(X);

% Calculate Pdx
Pdx = zeros(M,1);
for k=M:K
    Pdx = Pdx + (flipud(x(k-M+1:k)) * conj(d(k)));
end
Pdx = Pdx / N;

% Calculate Wiener filter (wo)
wo = inv(Rxx)*Pdx;

%% RLS---------------------
e_squared   = zeros(K,length(lambda_array)); % squared error
squared_dev = zeros(K,length(lambda_array)); % squared deviation

% Cycle through all values of lambda
for m=1:length(lambda_array)
    P      = eye(M);            % init inverse correlation matrix inv(R)
    w      = zeros(M,1);        % init weights to zero
    lambda = lambda_array(m);   % load lambda for this run
    
    for n=M:K
        xn = flipud(x(n-M+1:n));    % define input vector

        % RLS update
        pi = P*xn;
        k = pi / (lambda + ctranspose(xn)*pi);      % compute gain vector
        y = ctranspose(w)*xn;                       % compute filter output
        e = d(n) - y;                               % compute error
        w = w + k*conj(e);                          % update weight vector
        P = (P/lambda) - (k*ctranspose(xn)*P)/lambda; % update inverse correlation matrix inv(R)
                                                      % via the Ricatti equation

        e_squared(n,m)   = conj(e)*e;               % squared error
        squared_dev(n,m) = ctranspose(w-wo)*(w-wo); % squared deviation
    end
end

%% Plots----------------------------

figure(1);
plot(20*log10(e_squared(:,1)));
hold on
plot(20*log10(e_squared(:,2)),'color','red');
plot(20*log10(e_squared(:,3)),'color','green');
title('Squared Error (RLS)');
xlabel('Iteration number');
ylabel('squared error (dB)');
legend({'0.99','0.998','0.9995'});
grid on;

figure(2);
plot(20*log10(squared_dev(:,1)));
hold on;
plot(20*log10(squared_dev(:,2)),'color','red');
plot(20*log10(squared_dev(:,3)),'color','green');
title('Squared Deviation (RLS)');
xlabel('Iteration number');
ylabel('squared deviation (dB)');
legend({'0.99','0.998','0.9995'});
grid on;

