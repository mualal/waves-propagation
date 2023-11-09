clear;

% time step to save disp
save_time = 3;
write_video = false;

% lattice geometry
left_x = -100;
right_x = 100;
left_y = -100;
right_y = 100;

% parameters
m_1 = 1;
m_2 = 0.5;
c = 0.75;
omega = 1;
a = 1;
gamma = pi/4.5;

%k_1 = asin(omega/2 * sqrt(m_1/c)) * 2 / a;

syms k;
k_1 = double(vpasolve(m_1*omega^2==4*c*(sin(cos(gamma)*k)^2+sin(sin(gamma)*k)^2), k, 1))*2/a;

num_x = left_x:a:right_x;
num_y = left_y:a:right_y;

mass = zeros(int32((-left_y+right_y)/a),int32((-left_x+right_x)/a));
for j=1:length(num_x)
    for i=1:length(num_y)
        if j <= -left_x / a
            mass(i,j) = m_1;
        else
            mass(i,j) = m_2;
        end
    end
end

% interface slope
% cur_step = 0;
% angle_step = 3;
% for j=-int32(left_x/a):length(num_x)
%     for i=cur_step+1:length(num_y)
%         mass(i,j) = m_1;
%     end
%     cur_step = cur_step + angle_step;
% end

% interface coordinates finder
x_coord = zeros(1,length(num_y)-1);
y_coord = zeros(1,length(num_y)-1);
temp_count=1;
for i=1:length(num_y)
    for j=1:length(num_x)-1
        if mass(i,j) ~= mass(i,j+1)
            x_coord(temp_count) = j*a+left_x;
            y_coord(temp_count) = (i-1)*a+left_y;   
        end
    end
    temp_count = temp_count + 1;
end

disp = zeros(int32((-left_y+right_y+1)/a),int32((-left_x+right_x+1)/a));
vel = zeros(int32((-left_y+right_y+1)/a),int32((-left_x+right_x+1)/a));

beta_x = 0.1;
beta_y = 0.1;
n_0 = -30;
v_drift = 25;
u_0 = 1;
%g_1 = a * sqrt(c/m_1) * cos(k_1*a/2);
g_1 = sqrt(c/m_1)*(a*cos(gamma)*cos(cos(gamma)*k_1*a/2)*sin(cos(gamma)*k_1*a/2)-a*sin(gamma)*cos(sin(gamma)*k_1*a/2)*sin(sin(gamma)*k_1*a/2))/(sqrt((cos(sin(gamma)*k_1*a/2))^2+(sin(cos(gamma)*k_1*a/2))^2));

% initial conditions
for i=1:length(num_x)
    for j=1:length(num_y)
        %disp(j, i) = u_0 * exp(-beta_x^2/2 * (num_x(i) - n_0)^2) * exp(-beta_y^2/2 * (num_y(j))^2);
        %disp(j, i) = disp(j, i) * sin(num_x(i) * k_1);
        disp(j, i) = u_0 * exp(-beta_x^2/2 * (num_x(i)*cos(gamma)+num_y(j)*sin(gamma)-n_0*cos(gamma)+v_drift)^2);
        disp(j, i) = disp(j, i) * exp(-beta_y^2/2 * (-num_x(i)*sin(gamma)+num_y(j)*cos(gamma)+n_0*sin(gamma)+v_drift)^2);
        disp(j, i) = disp(j, i) * sin(k_1*cos(gamma)*num_x(i)+k_1*sin(gamma)*num_y(j));
        vel(j, i) = -u_0 * exp(-beta_x^2/2 * (num_x(i)*cos(gamma)+num_y(j)*sin(gamma)-n_0*cos(gamma)+v_drift)^2);
        vel(j, i) = vel(j, i) * exp(-beta_y^2/2 * (-num_x(i)*sin(gamma)+num_y(j)*cos(gamma)+n_0*sin(gamma)+v_drift)^2);
        vel(j, i) = vel(j, i) * (omega*cos(k_1*cos(gamma)*num_x(i)+k_1*sin(gamma)*num_y(j)) - beta_x^2*g_1/a*(num_x(i)*cos(gamma)+num_y(j)*sin(gamma)-n_0*cos(gamma)+v_drift)*sin(k_1*cos(gamma)*num_x(i)+k_1*sin(gamma)*num_y(j)));
        %vel(j, i) = vel(j, i) * (omega * cos(k_1*(num_x(i))) - beta_x^2*g_1/a*(num_x(i)-n_0)*sin(num_x(i) * k_1));
    end
end

figure('Position',[50,50,850,650]); hold on
[X,Y] = meshgrid(num_x,num_y);
surf(X/a,Y/a,disp,'FaceAlpha',0.9,'EdgeAlpha',0.5);
plot3(x_coord/a,y_coord/a,max(max(disp))*ones(1,length(x_coord)),'LineWidth',2,'Color','Black')
title('Цветовая карта перемещений в момент времени t = 0');
xlabel('Номер частицы по оси Ox');
ylabel('Номер частицы по оси Oy');
colorbar;
pbaspect([1,(right_y-left_y)/(right_x-left_x),1]);
view(17,22);
hold off


dt = 0.01;
t_max = 110;
times = 0:dt:t_max;

result = cell(1, fix(t_max/save_time));

% equations of motion integration
for t=times
    disp_copy = disp;
    %disp(int32(length(num_y)/a), int32(length(num_x)/a)) = sin(omega*t);
    %vel(int32(length(num_y)/a), int32(length(num_x)/a)) = omega*cos(omega*t);
    vel = vel + c./mass.*(circshift(disp,[-1 0])+circshift(disp,[1 0])+circshift(disp,[0 -1])+circshift(disp,[0 1])-4*disp).*dt;
    disp = disp + vel.*dt;
    if rem(t, save_time) == 0
        result{fix(t/save_time)+1} = disp;
    end     
end


figure('Position',[50,50,850,650]); hold on
s1 = surf(X/a,Y/a,disp,'FaceAlpha',0.9,'EdgeAlpha',0.7);
title('Цветовая карта перемещений в момент времени t = '+string(t_max));
xlabel('Номер частицы по оси Ox');
ylabel('Номер частицы по оси Oy');
l1 = plot3(x_coord/a,y_coord/a,10*ones(1,length(x_coord)),'LineWidth',2,'Color','Black');
colorbar;
caxis([-1 1])
%s1.EdgeColor = 'none';
%view(17,22);
pbaspect([1,(right_y-left_y)/(right_x-left_x),1]);
for i = 1:length(result)
    s1.XData = X/a;
    s1.YData = Y/a;
    s1.ZData = result{i};
    l1.ZData = max(max(result{i}))*ones(1,length(x_coord));
    title('Цветовая карта перемещений в момент времени t = '+string((i-1)*save_time));
    if write_video
        filename = 'wave-with-angle-to-interface-1.gif';
        frame = getframe(2); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        if i == 1
           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else
           imwrite(imind,cm,filename,'gif','DelayTime',0.3,'WriteMode','append'); 
        end
    end
    pause(0.3)
end
hold off

