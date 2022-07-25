function [] = performanceproblemsolve(A,f,t0,x0,p,r,a,b,c,x11,x12)
figure();
plot(x0(1), x0(2), 'r*');
hold on;

if (a*abs(x0(1)-x11) + b*abs(x0(2)-x12) <= c)
    disp('Dont need no moving');
else
    tmin = 7;
    t_0 = t0:0.1:tmin;
    n = 50;
    k = 0:(2*pi)/n:2*pi;
    psi01 = cos(k);
    psi02 = sin(k);
    for i = 1:n
        ipsi1 = psi01(i)*exp(-A(1,1)*(t_0-t0)) + psi02(i)*exp(-A(2,1)*(t_0-t0))
        ipsi2 = psi01(i)*exp(-A(1,2)*(t_0-t0)) + psi02(i)*exp(-A(2,2)*(t_0-t0))
        u1 = p(1) + 1 + r.*ipsi1./(9.*sqrt(r*(ipsi1.^2)./9 +r.*(ipsi2.^2)./16));
        u2 = p(2) + r.*ipsi2./(16.*sqrt(r*(ipsi1.^2)./9 +r.*(ipsi2.^2)./16));
        x_mopt = zeros(2, size(u1, 2));
        x_mopt(1,1) = x0(1);
        x_mopt(2,1) = x0(2);
        x_opt = x_mopt;
        for j = 2:size(u1, 2)
            [t, y] = ode45(@(t, y) syss(t, y, A, [u1(j); u2(j)], f), [t0 t_0(j)], x0);
            x_mopt(1, j) = y(end, 1);
            x_mopt(2, j) = y(end, 2); 
            if (a*abs(x_mopt(1,j) - x11)+b*abs(x_mopt(2,j) - x12) <= c)
                break;
            end
        end
        if (j <= size(u1, 2) && t_0(j) < tmin)
            tmin = t_0(j);
            x_opt = x_mopt;
            u_opt = [u1;u2 ]
            jmin = j;
        end
        plot(x_mopt(1, 1:j), x_mopt(2, 1:j), 'black');
        
    end
    plot(x_mopt(1, 1: j), x_mopt(2, 1:j), 'g');
    x = [x11+c/a, x11, x11-c/a, x11];
    y = [x12, x12-c/b,x12, x12+c/b ];
    fill(x,y,'r');
    tmin
    hold off;
end

end
