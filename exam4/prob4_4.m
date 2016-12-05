% Exam 4 Problem 4

clear all;
close all;
hold off;

load p4.mat; % loads x

%segment_lengths = [1 6 12 24 48];
segment_lengths = [12 48];

figure(1);

% Cycle through all segment lengths
for idx=1:length(segment_lengths)
    K = segment_lengths(idx);
    M = floor(length(x)/K); % number of bins in this PSD
    w = linspace(-pi,pi,M);
    degrees = linspace(-180,180,M);
    Pxx = zeros(M,1);       % PSD for each segment
    
    % Cycle through all K segments
    for i=0:K-1
        % Compute ith Pxx(f) over all f
        for f=1:M
            sum = 0;
            for n=0:M-1
                sum = sum + x(n+1+(i*M))*exp(-j*w(f)*n);
            end
            Pxx(f) = Pxx(f) + (sum * conj(sum))/M;
        end
    end
    
    % Normalize by number of segments
    Pxx = Pxx / K;
    
    % Plot Power Spectral Density
    plot(degrees,20*log10(abs(Pxx)));
    hold on;
end

title('Bartlett Power Spectral Density');
xlabel('Frequency (degrees)');
ylabel('Magnitude (dB)');
legend({'1','6','12','24','48'});
grid on;

