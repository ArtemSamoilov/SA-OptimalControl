% function that solves DE to find y(t) through v(t)
function [dy] = finalsystem(t, y, C, D, v, g)
   dy = C*y+D*v+g;
end