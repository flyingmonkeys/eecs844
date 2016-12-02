% Exam 4 Problem 4

clear all;
close all;
hold off;

load p4.mat; % loads x

segment_lengths = [1 6 12 24 48];
%segment_lengths = [24 48];
% PSD for each segment
Pxx = zeros(length(x)/min(segment_lengths),max(segment_lengths));

figure(1);

% Cycle through all segment lengths
for idx=1:length(segment_lengths)
    K = segment_lengths(idx);
    M = floor(length(x)/K); % number of bins in this PSD
    w = linspace(-pi,pi,M);
        
    % Cycle through all K segments
    for i=0:K-1
        % Compute ith Pxx(f) over all f
        for f=1:M
            sum = 0;
            for n=0:M-1
                sum = sum + x(n+1+(i*M))*exp(-j*w(f)*n);
            end
            Pxx(f,i+1) = (1/M) * sum * conj(sum);
        end
    end
    
    % Now take the average of all periodograms to get the Bartlett average
    PxxB = zeros(M,1);
    for f=1:M      % cycle through all frequencies
        sum = 0;
        for i=1:K  % cycle through all segments for this frequency
            sum = sum + Pxx(f,i);
        end
        PxxB(f) = sum / K;
    end
    
    % Power Spectral Density
    plot(w,20*log10(abs(PxxB)));
    hold on;
end

title('Bartlett Power Spectral Density');
xlabel('Frequency (rad)');
ylabel('Magnitude (dB)');
legend({'1','6','12','24','48'});
grid on;

