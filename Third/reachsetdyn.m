function reachsetdyn(alpha, t1, t2, N, filename)
    mov(1:N) = struct('cdata',[],'colormap',[]);
    for i=1:N
        tau = t1 + ((t2 - t1).^(i))./N;
        [x, y, x1, y1] = reachset(alpha, tau);
        
        plot(x,y,'Color','b','LineWidth',3);
        
        mov(i) = getframe();
    end
    %mov1 = VideoWriter('a.avi');
    %open(mov);
    %writeVideo(mov1,mov);
    %close(mov1);
end