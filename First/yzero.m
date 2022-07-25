function [value,isterminal,direction] = yzero(~,y)
    value = y(1);
    isterminal = 1;
    direction = -1;
end