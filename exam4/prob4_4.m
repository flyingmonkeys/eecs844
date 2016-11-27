% Problem 4, Exam 4

clear all;

load p4.mat; % loads x

segment_lengths = [1 6 12 24 48];
R = 200; % number of PSD bins
Pxx = zeros(R,max(segment_lengths));
w = linspace(-pi,pi,R);

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
    for f=1:R
        sum = 0;
        for i=1:K
            sum = sum + Pxx(f,i);
        end
        PxxB(f) = sum / K;
    end
    
    theta = linspace(-pi,pi,R);
%theta = 1:length(PxxB);

    % MVDR Power Spectrum (electrical angle)
    figure(idx);
    plot(theta,20*log10(abs(PxxB)));
    title('MVDR Power Spectrum');
    xlabel('Frequency (rad)');
    ylabel('Magnitude (dB)');

end


