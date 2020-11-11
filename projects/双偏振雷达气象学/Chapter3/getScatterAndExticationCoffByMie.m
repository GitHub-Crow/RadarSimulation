function [ext, sca] = getScatterAndExticationCoffByMie(D, lambda, permittivity)
% Calculate normalized extincttion and scatter cross section
alpha = pi.*D/lambda;
ref = sqrt(permittivity); % refractive index
kext = zeros(size(alpha)); ksca = zeros(size(alpha));
for n = 1 : 50 % times of sum
    syms x f1 f2
    if n == 1 % init parameters
        f1 = sin(x)./x -cos(x);
        f1dot = diff(f1, x, 1);
        f2 =(sin(x) + 1j*cos(x))./x - cos(x) + 1j*sin(x);
        f2dot = diff(f2, x, 1);
    else
        f1 =(pi*x/2)^(1/2).*besselj(n+0.5, x);
        f1dot = (pi*x/2)^(1/2).*besselj(n-0.5, x) - n*(pi*x/2)^(1/2).*besselj(n + 0.5, x)./x;
        f2 = (pi*x/2)^(1/2).*(besselj(n+0.5, x) - 1j*bessely(n+0.5, x));
        f2dot = (pi*x/2)^(1/2).*(besselj(n-0.5, x) - 1j*bessely(n-0.5, x)) - ...
                n/x*(pi*x/2)^(1/2).*(besselj(n+0.5, x) - 1j*bessely(n+0.5, x));
    end
    g1 = matlabFunction(f1); % transform asign function to numeric function
    g1dot = matlabFunction(f1dot);
    g2 = matlabFunction(f2);
    g2dot = matlabFunction(f2dot);
    an = (g1dot(ref*alpha).*g1(alpha) - ref*g1dot(alpha).*g1(ref*alpha))./...
         (g1dot(ref*alpha).*g2(alpha) - ref*g2dot(alpha).*g1(ref*alpha));
    bn = (g1dot(ref*alpha).*g1(alpha) - ref*g1dot(alpha).*g1(ref*alpha))./...
         (ref*g1dot(ref*alpha).*g2(alpha) - g2dot(alpha).*g1(ref*alpha));
    ksca = ksca + (2*n + 1)*(abs(an).^2 + abs(bn).^2);
    kext = kext + (2*n + 1)*real(an + bn);
end
ext = 2./alpha.^2.*kext;
sca = 2./alpha.^2.*ksca;