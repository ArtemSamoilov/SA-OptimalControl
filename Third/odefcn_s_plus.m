function dydt = odefcn_s_plus(t,y,alpha)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = alpha - y(1).^2 - sin(y(2)) - y(2).*sin(y(1).^3);
end