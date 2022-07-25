% function that triggers the rest
function [] = solvelinear6(A,B,f, t0, xi1, p, rs,xs,ys, a,b,c);
figure();
daspect([1 1 1]);
hold on;

xconv = [xs(1)-rs(1)/2, xs(1)-rs(1)/2, xs(1)+rs(1)/2, xs(1)+rs(1)/2, xs(2)-rs(2)/2, xs(2)-rs(2)/2, xs(2)+rs(2)/2, xs(2)+rs(2)/2, xs(3)-rs(3)/2, xs(3)-rs(3)/2, xs(3)+rs(3)/2, xs(3)+rs(3)/2]';
yconv = [ys(1)-rs(1)/2, ys(1)+rs(1)/2, ys(1)+rs(1)/2, ys(1)-rs(1)/2, ys(2)-rs(2)/2, ys(2)+rs(2)/2, ys(2)+rs(2)/2, ys(2)-rs(2)/2, ys(3)-rs(3)/2, ys(3)+rs(3)/2, ys(3)+rs(3)/2, ys(3)-rs(3)/2]';
global setconv;
setconv = delaunayTriangulation(xconv,yconv);
global Conv;
Conv = convexHull(setconv);
xlabel('x_1');
ylabel('x_2');
axis([-4 1 -2 6]);
%axis([-100 5 -60 10]);
plot(setconv.Points(:,1),setconv.Points(:,2),'ro');
hold on;
fill(setconv.Points(Conv,1),setconv.Points(Conv,2),'y');
hold on;
plot(setconv.Points(Conv,1),setconv.Points(Conv,2),'r');
hold on;
plot(xi1(1), xi1(2), 'b*');
hold on;

N = 100;
tmin = 1;
C = -A;
D = -B;
g = -f;
s0 = -tmin;
s1 = 1;
lin = 1:N;
linpsi = [sin(2*pi*lin./N); cos(2*pi*lin./N)];
tspan = t0:0.01:tmin;
sspan = s0:0.01:0;

maxsquare = @(r,x_1,x_2,l1,l2) (r/2)*(abs(l1)+abs(l2))+(l1*x_1+l2*x_2);
number = 0;
the_one = 1;
the_j = 1;
k = 0.1;
s_theone = 0;
flag = false;
    
check = 1;
for i = 1:N
    psi0 = linpsi(:,i);
    psi = @(t) expm((C').*t)*expm((-C').*t0)*psi0;
 
    MaxStep = 1e-2;
    options = odeset('Events', @stopintegr,'MaxStep', MaxStep);    
    
    tau = 1;
    [t, y, te, xe, ie] = ode45(@(t, y) finalsystem2(t, y, C, D, g, p, tau, a,b,c, psi), [s0 sspan(end)], xi1, options);
    for j = 1:size(t,1)
        ppsi(:,j) = psi(t(j));
    end
    l = D'*ppsi;
    v1 = p(1) + sign(l(1,:))*(tau^(1/2)*(a/c)^(-1/4));
    v2 = p(2) + sign(l(1,:))*2.*(l(2,:)./l(1,:)).*(tau^(3/2)*(a/c)^(1/4)*(b/c)^(-1));
    for j = 1:size(y, 1)
        x_mopt(1, j) = y(j, 1);
        x_mopt(2, j) = y(j, 2);
        if(inpolygon(x_mopt(1,j),x_mopt(2,j),setconv.Points(Conv,1),setconv.Points(Conv,2)))
            flag = 1;
            break;
        end
        if(abs(x_mopt(1,j))>100 || abs(x_mopt(2,j))>100)
            break;
        end
    end
    if(flag)
        if (abs(s_theone) < abs(te))
            %if(s1 < 0)
            te
                s_theone = te;
            %end
            the_t = t'+1;
            the_one = x_mopt;
            the_j = j;
            the_psi = [ppsi(1,end), ppsi(2,end)];
            pppsi = ppsi;
            v_opt = [v1; v2].*k;
        end
    end
    if((j <= size(v1,2)))
        s1 = t(end);
    end
    plot(x_mopt(1, 1:j), x_mopt(2, 1:j), 'black','LineWidth' ,0.1);

end
if(flag)
    plot(the_one(1, 1:the_j), the_one(2, 1:the_j), 'g', 'LineWidth' ,2);
    
    ppsi = -the_psi/norm(the_psi);
    rho(1) = p(1) + sign(ppsi(1))*(tau^(1/2)*(a/c)^(-1/4));
    rho(2) = p(2) + sign(ppsi(1))*2.*(ppsi(2)./ppsi(1)).*(tau^(3/2)*(a/c)^(1/4)*(b/c)^(-1));
    ppsi = -ppsi;
    the_x(1) = the_one(1,the_j);
    the_x(2) = the_one(2,the_j);
    delta = abs(dot(the_x, ppsi-rho));
    
    figure()
    plot(the_t, pppsi(1,end-size(the_t,2)+1:end), the_t, pppsi(2,end-size(the_t,2)+1:end));
    xlabel('t');
    ylabel('\psi(t)');
    legend('\psi_1(t)','\psi_2(t)');
    
    figure()
    plot(the_t, v_opt(1,1:size(the_t,2)), the_t, v_opt(2,1:size(the_t,2)));
    xlabel('t');
    ylabel('u(t)');
    legend('u1(t)','u2(t)');
end
x_mopt(1, j);
x_mopt(2, j);
time = 1+s_theone
1/delta
end

function [value,isterminal,direction] = stopintegr(~,x)
    value = event(x);
    isterminal = 1;
    direction = 0;
end

function res = event(x)
    global Conv;
    global setconv;
    res = inpolygon(x(1),x(2),setconv.Points(Conv,1),setconv.Points(Conv,2));
end