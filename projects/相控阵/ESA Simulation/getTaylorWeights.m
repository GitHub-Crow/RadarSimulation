function [weights] = getTaylorWeights(nele, sll, nbar)
r = 10^(abs(sll)/20);
a = log(r+(r*r-1)^0.5)/pi;
sigma2 = nbar^2/(a^2 + (nbar - 0.5)^2);
F = ones(1, nbar-1);
for m=1:(nbar-1)
    for n=1:(nbar-1)
        F(m) = F(m)*(1-m*m/sigma2/(a*a+(n-0.5)*(n-0.5)));
        if n ~= m
            F(m) = F(m)/(1 - m*m/n/n);
        end
    end
    F(m) = ((-1)^(m+1))/2*F(m);
end
jj = [1:nele]';
xx = (jj-1+0.5)/nele - 1/2;
mm = [1:nbar-1];
W = 1 + 2*cos(2*pi*xx*mm)*F';
WPK = 1 + 2*sum(F);
weights = W / WPK;

% r = 10^(abs(sll)/10);
% a = log(r+(r*r-1)^0.5)/pi;
% deltaSq = nbar^2/(a^2 + (nbar - 0.5)^2);
% SM = ones(1, nbar);
% SM(end) = 0;
% SMweight = 1;
% for m = 1 : nbar - 1
%     SMweight = SMweight*(nbar - m)/(nbar + m - 1);
%     for i = 1 : nbar - 1
%         SM(m) = SM(m)*(1 - (m^2/(deltaSq*(a^2 + (i - 0.5)^2))));
%     end
%     SM(m) = SM(m)*SMweight;    
% end
% 
% p = horzcat([1:1:ceil(nele/2)], [floor(nele/2):-1:1]); %#ok<*NBRAK>
% if mod(nele, 2) == 0
%     p = pi*(2*p - 1)/nele;
% else
%     p = 2*pi*(p - 1)/nele;
% end
% mvec = transpose(linspace(1, nbar, nbar));
% weights = 1 + 2*SM*cos(mvec*p);
end    
    