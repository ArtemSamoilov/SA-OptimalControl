function dxdt = fxpsi(x, psi, const)
    ut = u(x, psi, const);
    dxdt = fx([x(1); x(2); x(3)], const, ut);
    dpsi1dt = const.psi0 + psi(1) * (2 * const.k * x(1) - ut) / x(2);
    dpsi2dt = psi(1) * ...
        (-const.k * x(1) * x(1) + ut * (const.l + x(1))) / ...
        (x(2) * x(2));
    dxdt = [dxdt(1); dxdt(2); dxdt(3); dpsi1dt; dpsi2dt];
end