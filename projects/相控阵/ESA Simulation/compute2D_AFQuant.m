function [AF] = compute2D_AFQuant(amWgts, dx, dy, f, f0, randAmErr, randPhsErr, u, u0, v, v0,  bits)
lambda = 0.3/f;
lambda0 = 0.3/f0;
k = 2*pi/lambda;
k0 = 2*pi/lambda0;
[nxEle, nyEle] = size(amWgts);
AF = zeros(size(u));

LSB = 360/2^bits;
randErr = randAmErr.*exp(1j*randPhsErr*pi/180);
for ii = 1 : nxEle
    xidx = (ii - (nxEle + 1)/2)*dx;
    for jj = 1 : nyEle
        yidx = (jj  - (nyEle + 1)/2)*dy;
        phase = k*u*xidx + k*v*yidx;
        phase0 = k0*u0*xidx + k0*v0*yidx;
        phase0 = pi/180*round(180/pi*phase0/LSB)*LSB;
        AF = AF + amWgts(ii, jj).*randErr(ii, jj).*...
             exp(1j*(phase - phase0));
    end
end

end

