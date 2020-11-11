function [IntGain, IdealGain] = compute2D_IntGain(pat, f, nxEle, nyEle, dx, dy, theta, phi)
thetaN = length(theta); phiN = length(phi);
dtheta = (theta(end) - theta(1))/thetaN;
dphi = (phi(end) - phi(1))/phiN;

theta = theta'*ones(1, phiN);

IntGain = sum(sum(dtheta*dphi*sin(theta)))*pat.^2/...
          sum(sum(dtheta*dphi*sin(theta).*sin(theta)));
area = nxEle*dx*nyEle*dy;
lambda = 0.3/f;
IdealGain = 4*pi*area/lambda^2;
end

