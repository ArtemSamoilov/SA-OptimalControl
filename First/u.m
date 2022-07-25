function res = u(x, psi, const)
    if (x(2) - const.M <= 0)
        res = 0;
        return
    end
    curG = G(x, psi, const);
    if (abs(curG) < const.Gtol)        
        t = zeros(6, 1);
        t(1) = 2 * const.g * const.psi0 * x(2) * x(2);
        t(2) = 2 * const.g * const.k * const.l * psi(1) * x(2);
        t(3) = 6 * const.g * const.k * psi(1) * x(1) * x(2);
        t(4) = -const.k * const.psi0 * x(1) * x(1) * x(2);
        t(5) = -2 * const.k * const.k * const.l * psi(1) * x(1) * x(1);
        t(6) = -2 * const.psi0 * const.k * const.l * x(1) * x(2);
        
        p = zeros(6, 1);
        p(1) = const.g * x(2) * psi(1);
        p(2) = 2 * const.k * const.l * const.l * psi(1);
        p(3) = 6 * const.k * const.l * x(1) * psi(1);
        p(4) = 4 * const.k * x(1) * x(1) * psi(1);
        p(5) = const.l * const.psi0 * x(2);
        p(6) = const.psi0 * x(1) * x(2);
        
        pp = sum(p);
        tt = sum(t);
        if abs(pp) < 1e-5
            pp = 1e-5 * sign(pp);
            if pp == 0
                pp = 1e-5;
            end
        end
        res = tt / pp;
        res = min(res, const.umax);
        res = max(res, const.umin);
    elseif (curG > 0)
        res = const.umax;
    elseif (curG < 0)
        res = const.umin;
    end
end

function res = G(x, psi, const)
    res = psi(1) * (const.l + x(1)) - psi(2) * x(2);
end