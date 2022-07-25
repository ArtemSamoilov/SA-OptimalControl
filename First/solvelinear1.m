% function that triggers the rest
function [] = solvelinear1(A,B,f, t0, xi1, p, rs,xs,ys, a,b,c);
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
axis([-5 15 0 17]);
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
sspan = s0:0.01:s1;

maxsquare = @(r,x_1,x_2,l1,l2) (r/2)*(abs(l1)+abs(l2))+(l1*x_1+l2*x_2);
number = 0;
the_one = 1;
the_j = 1;
s_theone = -1;
flag = false;
    
check = 1;
for i = 1:N
    psi0 = linpsi(:,i);
 
     MaxStep = 1e-1;
     options = odeset('Events', @stopintegr,'MaxStep', MaxStep);    
    
    [t, psi] = ode45(@(t,psi) iterpsi(t, A, psi), tspan, psi0, options);
    psi = psi';
    ppsi = psi;
    l = D'*psi;
    tau = 1;
    v1 = p(1) + sign(l(1,:))*(tau^(1/2)*(a/c)^(-1/4));
    v2 = p(2) + sign(l(1,:))*2.*(l(2,:)./l(1,:)).*(tau^(3/2)*(a/c)^(1/4)*(b/c)^(-1));
    x_mopt = zeros(2, size(v1, 2));
    x_mopt(1,1) = xi1(1);
    x_mopt(2,1) = xi1(2);
    x_opt = x_mopt;
    number = number +1;
    the_j = j;
    v = zeros(2, size(psi, 1));% fix zero value of B'*psi
    for i = 1 : size(psi, 1)
        t1 = [psi(i, 1); psi(i, 2)];
        t1 = D' * t1;
        v(1, i) = t1(1);
        v(2, i) = t1(2);
    end
    %if (s_theone < s1) && flag
        %if(s1 < 0)
            s_theone = s1;
        %end
        the_one = x_mopt;
        the_j = j;
        the_psi = [psi(1,end), psi(2,end)];
        v_opt = [v1; v2];
    %end
    if((j <= size(v1,2)))
        s1 = sspan(j);
        x_opt = x_mopt;
        jmin = j;
    end
    plot(x_mopt(1, 1:j), x_mopt(2, 1:j), 'black','LineWidth' ,0.1);

end
if(flag)
    plot(the_one(1, 1:the_j), the_one(2, 1:the_j), 'g', 'LineWidth' ,2);
    psi = -the_psi/norm(the_psi);
    rho(1) = p(1) + sign(psi(1))*(tau^(1/2)*(a/c)^(-1/4));
    rho(2) = p(2) + sign(psi(1))*2.*(psi(2)./psi(1)).*(tau^(3/2)*(a/c)^(1/4)*(b/c)^(-1));
    psi = -psi;
    the_x(1) = the_one(1,the_j);
    the_x(2) = the_one(2,the_j);
    delta = abs(dot(the_x, psi-rho));
    figure()
    plot(tspan(1:size(ppsi,2)), ppsi(1,1:end), tspan(1:size(ppsi,2)), ppsi(2,1:end));
    xlabel('t');
    ylabel('\psi(t)');
    legend('\psi_1(t)','\psi_2(t)');
end
x_mopt(1, j);
x_mopt(2, j);
time = tspan(size(ppsi,2))
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