% Exam 4 Problem 4

clear all;
close all;
hold off;

load p4.mat; % loads x

segment_lengths = [1 6 12 24 48];
R = 200; % number of PSD bins
Pxx = zeros(R,max(segment_lengths)); % PSD for each segment
w = linspace(-pi,pi,R);

figure(1);

% Cycle through all segment lengths
for idx=1:length(segment_lengths)
    K = segment_lengths(idx);
    M = floor(length(x)/K);

    % Cycle through all K segments
    for i=0:K-1
        % Compute ith Pxx(f) over all f
        for f=1:R
            sum = 0;
            for n=0:M-1
                sum = sum + x(n+1+(i*M))*exp(-j*w(f)*n);
            end
            Pxx(f,i+1) = (1/M) * sum * conj(sum);
        end
    end
    
    % Now take the average of all periodograms to get the Bartlett average
    PxxB = zeros(R,1);
    for f=1:R      % cycle through all frequencies
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

