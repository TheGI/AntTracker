%%
clear all;
pathName = '/Users/GailLee/Documents/exampleVideos/CampoArea_open/num1/colony1/frames/ant1/savedFiles/';
load([pathName 'antData.mat']);
dimension = size(mouseData.coord);

%% shorten variable names
for i = 1:dimension(2)
    x(i,:) = mouseData.coord(:,i,1);
    y(i,:) = mouseData.coord(:,i,2);
    u(i,1) = 0;
    v(i,1) = 0;
    u(i,2:dimension(1)) = diff(mouseData.coord(:,i,1));
    v(i,2:dimension(1)) = diff(mouseData.coord(:,i,2));
    s(i,:) = sqrt(u(i,:).^2 + v(i,:).^2);
    l(i,:) = cumsum(s(i,:));
    seg = (0:l(i,end)/10) * 10;
    ind = arrayfun(@(x) find(l(i,:) > x, 1), seg);
    pu(i,:) = u(i,ind);
    pv(i,:) = v(i,ind);
    for j = 2:size(pu(i,:),2)
        uX1 = pu(i,j-1);
        uX2 = pu(i,j);
        uY1 = pv(i,j-1);
        uY2 = pv(i,j);
        a(i,j-1) = atan2d(uX1*uY2-uY1*uX2,uX1*uX2+uY1*uY2);
    end
    
    mouseData.vel(:,i,1:2) = [u(i,:);v(i,:)]';
    mouseData.speed(:,i,1) = s(i,:);
    mouseData.pathLength(:,i,1) = cumsum(s(i,:));
    mouseData.velPiece(:,i,1:2) = [pu(i,:);pv(i,:)]';
    mouseData.angle(:,i,1) = a(i,:);
end

%% make individual trajectory figures
for i = 1:dimension(2)
    %hold on;
    fig_trajectory = plot(mouseData.coord(:,i,1),mouseData.coord(:,i,2));
    saveas(fig_trajectory,[pathName 'trajectory_ant', num2str(i) '.png']);
    figure;
    %hold off;    
end
%% merge all trajectories into one plot
for i = 1:dimension(2)
    hold on;
    plot(mouseData.coord(:,i,1),mouseData.coord(:,i,2));
    hold off;    
end
saveas(gcf,[pathName 'trajectory_total.png']);

%%
for i = 1:dimension(2)
    
    
    velX = diff(mouseData.coord(:,i,1));
    velY = diff(mouseData.coord(:,i,2));
    speed = sqrt(velX.^2 + velY.^2);
    figure;
    fig_velocity = plot(velX,velY,'.');
    saveas(fig_velocity,[pathName 'velocity_ant',num2str(i),'.png']);
    figure;
    fig_speed = plot(speed);
    saveas(fig_speed,[pathName 'speed_ant',num2str(i),'.png']);
end
%%
for i = 1:dimension(2)
    velX = diff(mouseData.coord(:,i,1));
    velY = diff(mouseData.coord(:,i,2));
    speed = sqrt(velX.^2 + velY.^2);
    hold on;
    plot(velX,velY,'.');
    hold off;
end

%%
for i = 1:dimension(2)
    velX = diff(mouseData.coord(:,i,1));
    velY = diff(mouseData.coord(:,i,2));
    speed = sqrt(velX.^2 + velY.^2);
    hold on;
    plot(speed);
    hold off;
end

%%
for i = 1:dimension(2)
    velX = diff(mouseData.coord(:,i,1));
    velY = diff(mouseData.coord(:,i,2));
    
    pathLength = cumsum(sqrt(velX.^2 + velY.^2));
    piecewise = (0:max(pathLength)/10) * 10;
    index = arrayfun(@(x) find(pathLength > x, 1), piecewise);
    
    pieceX = mouseData.coord(index,i,1);
    pieceY = mouseData.coord(index,i,2);
    piecevX = diff(pieceX);
    piecevY = diff(pieceY);
    
    
    for j = 2:size(piecevX,1)

    u = [piecevX(j-1),piecevY(j-1)];
    v = [piecevX(j),piecevY(j)];
      angle(j-1) = atan2d(piecevX(j-1)*piecevY(j)-piecevY(j-1)*piecevX(j),...
          piecevX(j-1)*piecevX(j)+piecevY(j-1)*piecevY(j));

    end
    %hold on;
    figure;
    plot(angle);
    stdAngle(i) = std(angle);
    %hold off;
end
