function [value, isterminal, direction] = eventfcn(t,x)
    value = x(2);
    isterminal = 1;
    direction = 0;
end