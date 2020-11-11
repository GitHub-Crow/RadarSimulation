function AF = compute2D_AF(amWgts, dx, dy, f, f0, randErrorAm, randErrorPhs, u, u0, v, v0) 
lambda = 0.3/f;
lambda0 = 0.3/f0;
k = 2*pi/lambda;
k0 = 2*pi/lambda0;
[nxEle, nyEle] = size(amWgts);
AF = zeros(size(u));

randErr = randErrorAm.*exp(1j*randErrorPhs*pi/180);
for ii = 1 : nxEle
    xidx = (ii - (nxEle + 1)/2)*dx;
    for jj = 1 : nyEle
        yidx = (jj  - (nyEle + 1)/2)*dy;
        AF = AF + amWgts(ii, jj).*randErr(ii, jj).*...
             exp(1j*xidx*(k*u - k0*u0)).*exp(1j*yidx*(k*v - k0*v0));
    end
end

end