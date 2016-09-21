clear all;

load p6.mat

M = length(w_adap);
s = zeros(M,1);

theta_x = linspace(-pi,pi,20*M);

for idx=1:M
	s(idx) = exp(-j*(idx-1)*theta);
end

beampattern = ctranspose(s)*w_adap;