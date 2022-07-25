function [value,isterminal,direction] = stopIntegrating(~, y, const)
    value = (y(2) - const.M - 1e-5) * (y(1) > 0) * ...
        (abs(y(4)) < 1e15) * (abs(y(5)) < 1e15);
    isterminal = 1;
    direction = 0;
end