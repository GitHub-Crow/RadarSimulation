function [AF] = compute1D_AF(Am, elementCount, d, f0, f, u0, u)
lambda0 = 0.3/f0; lambda = 0.3/f;
k0 = (2*pi)/lambda0; k = (2*pi)/lambda;
Xm = linspace(1, elementCount, elementCount)*d - (elementCount + 1)/2*d;
AF = zeros(size(u));
for i = 1 : elementCount
    AF = AF + Am(i)*exp(1j*Xm(i)*(k*u - k0*u0));
end
end

