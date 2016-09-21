clear all;
close all;
hold off;

load p6.mat

M = length(w_adap);
s = zeros(M,1);
numsamps = 20*M;

theta_x = linspace(-pi,pi,numsamps);

for idx=1:M
	s(idx) = exp(-j*(idx-1)*theta_x(idx));
end

beampattern = ctranspose(s)*w_adap;