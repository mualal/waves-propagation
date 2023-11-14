% waves propagation (lattice-lattice)

%clc;
clear;

% time step to save data
save_time = 3;
write_gif = false;
gif_filename = '1.gif';
gif_delay = 0.3;
display_delay = 0.3;

% data to display
plot_var = 'energy';
energy_monitor = true;

% lattice geometry boundaries
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

% integration parameters
dt = 0.01;
t_max = 150;

%--------------------------------------------------------------------------

%k_1 = asin(omega/2 * sqrt(m_1/c)) * 2 / a;

syms k;
k_1 = double(vpasolve(m_1*omega^2==4*c*(sin(cos(gamma)*k)^2+...
    sin(sin(gamma)*k)^2), k, 1))*2/a;

num_x = round((left_x:a:right_x)/a);
num_y = round((left_y:a:right_y)/a);

mass = zeros(length(num_y),length(num_x));
for j=1:length(num_x)
    for i=1:length(num_y)
        if num_x(j) < 0
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
            x_coord(temp_count) = j+round(left_x/a)-1;
            y_coord(temp_count) = (i-1)+round(left_y/a);   
        end
    end
    temp_count = temp_count + 1;
end

disp = zeros(length(num_y),length(num_x));
vel = zeros(length(num_y),length(num_x));

beta_x = 0.1;
beta_y = 0.1;
n_0 = -35;
v_0 = -35;
u_0 = 1;
%g_1 = a * sqrt(c/m_1) * cos(k_1*a/2);
g_1 = sqrt(c/m_1)*(a*cos(gamma)*cos(cos(gamma)*k_1*a/2)*...
    sin(cos(gamma)*k_1*a/2)-a*sin(gamma)*cos(sin(gamma)*k_1*a/2)*...
    sin(sin(gamma)*k_1*a/2))/(sqrt((cos(sin(gamma)*k_1*a/2))^2+...
    (sin(cos(gamma)*k_1*a/2))^2));

% initial conditions
for i=1:length(num_x)
    for j=1:length(num_y)
        %disp(j, i) = u_0 * exp(-beta_x^2/2 * (num_x(i) - n_0)^2) * ...
        %    exp(-beta_y^2/2 * (num_y(j))^2);
        %disp(j, i) = disp(j, i) * sin(num_x(i) * k_1 * a);
        disp(j, i) = u_0 * exp(-beta_x^2/2 * (num_x(i)*cos(gamma)+...
            num_y(j)*sin(gamma)-n_0*cos(gamma)-v_0*sin(gamma))^2);
        disp(j, i) = disp(j, i) * exp(-beta_y^2/2 * (-num_x(i)*...
            sin(gamma)+num_y(j)*cos(gamma)+n_0*sin(gamma)-...
            v_0*cos(gamma))^2);
        disp(j, i) = disp(j, i) * sin(k_1*a*cos(gamma)*num_x(i)+...
            k_1*a*sin(gamma)*num_y(j));
        vel(j, i) = -u_0 * exp(-beta_x^2/2 * (num_x(i)*cos(gamma)+...
            num_y(j)*sin(gamma)-n_0*cos(gamma)-v_0*sin(gamma))^2);
        vel(j, i) = vel(j, i) * exp(-beta_y^2/2 * (-num_x(i)*sin(gamma)...
            +num_y(j)*cos(gamma)+n_0*sin(gamma)-v_0*cos(gamma))^2);
        vel(j, i) = vel(j, i) * (omega*cos(k_1*a*cos(gamma)*num_x(i)+...
            k_1*a*sin(gamma)*num_y(j)) - beta_x^2*g_1/a*(num_x(i)*...
            cos(gamma)+num_y(j)*sin(gamma)-n_0*cos(gamma))*...
            sin(k_1*a*cos(gamma)*num_x(i)+k_1*a*sin(gamma)*num_y(j)));
        %vel(j, i) = vel(j, i) * (omega * cos(k_1*(num_x(i))) - ...
        %    beta_x^2*g_1/a*(num_x(i)-n_0)*sin(num_x(i) * k_1 * a));
    end
end


times = 0:dt:t_max;

result_disp = cell(1, fix(t_max/save_time)+1);
result_vel = cell(1, fix(t_max/save_time)+1);
result_e = cell(1, fix(t_max/save_time)+1);

% equations of motion integration
for t=times
    %disp(int32(length(num_y)/a), int32(length(num_x)/a)) = sin(omega*t);
    %vel(int32(length(num_y)/a),int32(length(num_x)/a))=omega*cos(omega*t);
    if rem(t, save_time) == 0
        result_disp{fix(t/save_time)+1} = disp;
        result_vel{fix(t/save_time)+1} = vel;
        result_e{fix(t/save_time)+1} = energy(mass,c,vel,disp);
    end
    
    % leapfrog synchronized form
    acc1 = c./mass.*(circshift(disp,[-1 0])+circshift(disp,[1 0])+...
        circshift(disp,[0 -1])+circshift(disp,[0 1])-4*disp);
    disp = disp + vel*dt + 1/2*acc1*dt^2;
    acc2 = c./mass.*(circshift(disp,[-1 0])+circshift(disp,[1 0])+...
        circshift(disp,[0 -1])+circshift(disp,[0 1])-4*disp);
    vel = vel+1/2*(acc1+acc2)*dt;
    
    % leapfrog kick-drift-kick form
    %acc1 = c./mass.*(circshift(disp,[-1 0])+circshift(disp,[1 0])+...
    %    circshift(disp,[0 -1])+circshift(disp,[0 1])-4*disp);
    %vel_mid = vel + acc1*dt/2;
    %disp = disp + vel_mid*dt;
    %acc2 = c./mass.*(circshift(disp,[-1 0])+circshift(disp,[1 0])+...
    %    circshift(disp,[0 -1])+circshift(disp,[0 1])-4*disp);
    %vel = vel_mid + acc2*dt/2;
    
    % leapfrog self-made form
    %vel = vel + c./mass.*(circshift(disp,[-1 0])+circshift(disp,[1 0])+...
    %    circshift(disp,[0 -1])+circshift(disp,[0 1])-4*disp).*dt;
    %disp = disp + vel.*dt;
end


%--------------------------------------------------------------------------

% set results to show
result = result_disp;
label = "перемещений";
if strcmp(plot_var,'vel')
    temp = cell2mat(result_vel)./(a*omega);
    result = mat2cell(temp, length(num_y),...
        (repmat(length(num_x),1,length(result_vel))));
    %for i=1:length(result_vel)
    %    result{i} = result_vel{i} ./ (a * omega);
    %end
    %result = result_vel;
    label = "скоростей";
end
if strcmp(plot_var,'energy')
    result = result_e;
    label = "энергий";
end

descr_str = sprintf("\n (m_1=%.1f;   m_2=%.1f;   c=%.3f;   a=%.1f;"+...
    "   gamma=%.3f^o;   omega=%f)",m_1,m_2,c,a,gamma*180/pi,omega);

f1 = figure(1); hold on
f1.Position = [50,50,750,650];
[X,Y] = meshgrid(num_x,num_y);
surf(X,Y,result{1},'FaceAlpha',0.9,'EdgeAlpha',0.5);
plot3(x_coord,y_coord,max(max(result{1}))*ones(1,length(x_coord)),...
    'LineWidth',2,'Color','Black')
title('Цветовая карта '+label+' в момент времени t = 0'+descr_str);
xlabel('Номер частицы по оси Ox');
ylabel('Номер частицы по оси Oy');
colorbar;
axis([num_x(1) num_x(end) num_y(1) num_y(end)])
pbaspect([1,(right_y-left_y)/(right_x-left_x),1]);
ax = gca;
ax.Clipping = "off";
view(17,22);
hold off


if ~strcmp(plot_var,'energy')||(strcmp(plot_var,'energy')&&~energy_monitor)
f2=figure(2); hold on
f2.Position = [50,50,750,650];
s1 = surf(X,Y,result{end},'FaceAlpha',0.9,'EdgeAlpha',0.7);
title("Цветовая карта "+label+" в момент времени t = "+string(t_max)+...
    descr_str);
xlabel('Номер частицы по оси Ox');
ylabel('Номер частицы по оси Oy');
l1 = plot3(x_coord,y_coord,10*ones(1,length(x_coord)),'LineWidth',2,...
    'Color','Black');
colorbar;
%caxis([-1 1])
%caxis([-min(min(cell2mat(result))) max(max(cell2mat(result)))])
%s1.EdgeColor = 'none';
view(17,22);
axis([num_x(1) num_x(end) num_y(1) num_y(end)])
pbaspect([1,(right_y-left_y)/(right_x-left_x),1]);
ax = gca;
ax.Clipping = "off";
for i = 1:length(result)
    s1.XData = X;
    s1.YData = Y;
    s1.ZData = result{i};
    l1.ZData = max(max(result{i}))*ones(1,length(x_coord));
    title('Цветовая карта '+label+' в момент времени t = '+...
        string((i-1)*save_time)+descr_str);
    if write_gif
        filename = gif_filename;
        frame = getframe(2); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        if i == 1
           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else
           imwrite(imind,cm,filename,'gif','DelayTime',gif_delay,...
               'WriteMode','append'); 
        end
    end
    pause(display_delay)
end
hold off
end

if strcmp(plot_var,'energy') && ~energy_monitor
f3=figure(3); hold on
f3.Position = [50,50,750,650];
plot(0:save_time:t_max,cell2mat(cellfun(@(x) sum(sum(x)),result_e,...
    'UniformOutput',false)),'LineWidth',1.7,'Color','Black')
plot(0:save_time:t_max,cell2mat(cellfun(@(x) sum(sum(x(1:end,...
    1:-left_x))),result_e,'UniformOutput',false)),'LineWidth',1,...
    'Color','Blue')
plot(0:save_time:t_max,cell2mat(cellfun(@(x) sum(sum(x(1:end,...
    -left_x+1:right_x-left_x+1))),result_e,'UniformOutput',false)),...
    'LineWidth',1,'Color','Red')
title("Зависимость энергий в системе от времени"+descr_str);
xlabel('Время, усл.ед.');
ylabel('Энергия, усл.ед.');
legend('Полная энергия системы', 'Энергия левой решётки',...
    'Энергия правой решётки')
grid on;
grid minor;
end


if strcmp(plot_var,'energy') && energy_monitor
e_sum = cell2mat(cellfun(@(x) sum(sum(x)),result_e,'UniformOutput',false));
e_l = cell2mat(cellfun(@(x) sum(sum(x(1:end,1:-left_x))),result_e,...
    'UniformOutput',false));
e_r=cell2mat(cellfun(@(x) sum(sum(x(1:end,-left_x+1:right_x-left_x+1))),...
    result_e,'UniformOutput',false));
f2=figure(2); hold on
f2.Position = [0,50,1400,650];
til = tiledlayout(1,2,'TileSpacing','Compact');
nexttile; hold on
s1 = surf(X,Y,result{1},'FaceAlpha',0.9,'EdgeAlpha',0.7);
title('Цветовая карта энергий');
xlabel('Номер частицы по оси Ox');
ylabel('Номер частицы по оси Oy');
l1 = plot3(x_coord,y_coord,10*ones(1,length(x_coord)),'LineWidth',2,...
    'Color','Black');
colorbar;
%caxis([-1 1])
%caxis([-min(min(cell2mat(result))) max(max(cell2mat(result)))])
%s1.EdgeColor = 'none';
%view(17,22);
axis([num_x(1) num_x(end) num_y(1) num_y(end)])
pbaspect([1,(right_y-left_y)/(right_x-left_x),1]);
ax = gca;
ax.Clipping = "off";
hold off
nexttile; hold on
p1=plot(0,e_sum(1),'LineWidth',1.7,'Color','Black');
p2=plot(0,e_l(1),'LineWidth',1,'Color','Blue');
p3=plot(0,e_r(1),'LineWidth',1,'Color','Red');
ylim([min([e_sum, e_l, e_r]) 1.1*max([e_sum, e_l, e_r])])
title('Зависимость суммарных энергий в системе от времени');
xlabel('Время, усл.ед.');
ylabel('Энергия, усл.ед.');
legend('Полная энергия системы', 'Энергия левой решётки',...
    'Энергия правой решётки')
grid on;
grid minor;
hold off
til_t = title(til, "Энергия в момент времени t=0"+descr_str);
til_t.FontWeight = 'bold';
for i = 1:length(result)
    set(p1,'XData',(0:i-1)*save_time);
    set(p1,'YData',e_sum(1:i));
    set(p2,'XData',(0:i-1)*save_time);
    set(p2,'YData',e_l(1:i));
    set(p3,'XData',(0:i-1)*save_time);
    set(p3,'YData',e_r(1:i));
    drawnow;
    s1.XData = X;
    s1.YData = Y;
    s1.ZData = result{i};
    l1.ZData = max(max(result{i}))*ones(1,length(x_coord));
    title(til, "Энергия в момент времени t = "+string((i-1)*save_time)+...
        descr_str);
    if write_gif
        filename = gif_filename;
        frame = getframe(2); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        if i == 1
           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else
           imwrite(imind,cm,filename,'gif','DelayTime',gif_delay,...
               'WriteMode','append'); 
        end
    end
    pause(display_delay)
end
hold off
end


function e = energy(m,c,vel,disp)
    e = m./2 .* vel.^2 + c/4 * (circshift(disp,[-1 0])+...
        circshift(disp,[1 0])+circshift(disp,[0 -1])+...
        circshift(disp,[0 1])-4*disp).^2;
end
