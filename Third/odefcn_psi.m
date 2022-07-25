function dydt = odefcn_psi(t,y,alpha)
    dydt = zeros(4,1);
    dydt(1) = y(2);
    dydt(2) = alpha - y(1).^2 - sin(y(2)) - y(2).*sin(y(1).^3);
    dydt(3) = 2*y(4)*y(1) + y(4)*(sin(y(1)^3))^2 + 6*y(4)*y(1)^3*cos(y(1)^3)*sin(y(1)^3);
    dydt(4) = -y(3) + y(4)*cos(y(2));
end