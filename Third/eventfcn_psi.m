function [value, isterminal, direction] = eventfcn_psi(t,x)
    value = x(4);
    isterminal = 1;
    direction = 0;
end