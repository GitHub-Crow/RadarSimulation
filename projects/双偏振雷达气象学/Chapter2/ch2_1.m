% M-P模型：N(D) = N0*EXP(-lamada*D)假定N0=8*10^6
% 拟合lambda：lambda = (gamma(n+1)*N0/Mn)^(1/(n + 1));
N0 = 8*10^6;
W = 1;
M3 = W/(pi/6*10^-3);
DMAX = 6;
lambda = (gamma(4)*N0/M3)^(1/4);
DSD4 = @(D)(D.^4*N0.*exp(-lambda.*D));
DSD3 = @(D)(D.^3*N0.*exp(-lambda.*D));
DSD1 = @(D)(D.^1*N0.*exp(-lambda.*D));
M4 = integral(DSD4, 0, DMAX);
M1 = integral(DSD1, 0, DMAX);
Dm = M4/M3;
% 利用二分法寻找质量加权平均直径
PRECISION = 1;
left = 0; right = DMAX; M3_0 = 0; mid = (left + right)/2;
while abs(M3_0 - M3/2) > PRECISION
    mid = (left + right)/2;
    M3_0 = integral(DSD3, 0, mid);
    if M3_0 < M3/2
        left = mid;
    else
       right = mid;
    end
end
D0 = mid;
disp(['slop: ', num2str(lambda)]);
disp(['number concentrations: ', num2str(M1)]);
disp(['average diameter: ', num2str(D0)]);
disp(['content weighted diameter: ', num2str(Dm)]);