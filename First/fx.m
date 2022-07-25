function dxdt = fx(x, const, ut)
    dx1dt = -(const.g + (const.k * x(1) * x(1) - ...
              ut * (const.l + x(1))) / x(2));
    dx2dt = -ut;
    dx3dt = x(1);
    dxdt = [dx1dt; dx2dt; dx3dt];
end