% function that solves DE to find y(t) through v(t)
function [dy] = finalsystem2(t, y, C, D, g, p, tau, a,b,c, psi)
    psit = psi(t);
    l = D'*psit;
    u1 = p(1) + sign(l(1))*(tau^(1/2)*(a/c)^(-1/4));
    u2 = p(2) + sign(l(1))*2.*(l(2)./l(1)).*(tau^(3/2)*(a/c)^(1/4)*(b/c)^(-1));
    u = [u1; u2];
    dy = C*y+D*u+g;
end