function [dxdy] = syss(t, y, A, u, f)
dxdy = A*y + u + f;
end

