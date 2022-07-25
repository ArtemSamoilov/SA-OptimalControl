function [v, h, t] = restFly(v1, h1, t1, const)% fly until v > 0
    C1 = (atan(v1 / const.R1) / const.R2 + t1) / const.M;% calculate constant
    t = linspace(t1, min(const.T, C1 * const.M), const.N)';
    % v = 0 <=> t = C1 * M
    f = @(a) const.R1 * tan(const.R2 * (C1 * const.M - a));
    v = f(t);
    dt = [t; 0] - [0; t];
    dt(1) = 0;
    dt(end) = [];
    dh = ([v; 0] + [0; v]) / 2;
    dh(end) = [];
    h = dh .* dt;
    h(1) = h1;
    h = cumsum(h);
end