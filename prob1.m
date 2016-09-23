clear all;
hold off;
close all;

M = 50;
a = 0.99;
w = transpose([1+2i 2-5i 3-3i 4+i]);
b = transpose([1 1 1 1]);
R = diag(ones(length(w),1));
J = zeros(M,1);
dJ0 = zeros(M,4);
x = J;

idx = 1;
for k=-M/2:(M/2)-1
    w(1) = 1+(i*k);
%w(1) = k;
%    J(idx) = ctranspose(w)*R*w;
    J(idx) = 1/(a^(0.5*ctranspose(w)*R*w));
%    J(idx) = norm(w)^3;
%    J(idx) = real(ctranspose(w)*b);
    x(idx) = k;
    idx = idx + 1;
end

figure(1);
plot(x,J);

figure(2);
plot(x(1:M-1),diff(J));

idx = 1;
w_old = transpose([1+((M/2)*i) 2-5i 3-3i 4+i]);
for k=-M/2:(M/2)-1
    w(1) = conj(1+(i*k));
%w(1) = k;
    dw = w - w_old;
    w_old = w;
%    dJ0(idx,:) = R*w*dw(1) + R'*conj(w)*conj(dw(1));
%    dJ0(idx,:) = 1.5*norm(w)*w;
%    dJ0(idx,:) = 0.5*b;
dJ0(idx,:) = 0.25*log(a)*(a^(0.5*ctranspose(w)*R*w))*R*w;
    idx = idx + 1;
end

figure(3);
plot(x,imag(dJ0(:,1)));
