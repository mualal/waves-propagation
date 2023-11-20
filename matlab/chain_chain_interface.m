% waves propagation (chain-chain)

%clc;
clear;

%% Auxiliary Parameters

save_time = 10;
write_gif = false;
gif_filename = '1.gif';
anim_delay = 0.3;
gif_delay = 0.3;


%% Input Parameters

% chain geometry boundaries
left_x = -400;
right_x = 600;

% chain-chain parameters
m_1 = 1.0;
m_2 = 0.5;
c_1 = 1.0;
c_2 = 1.0;
c_12 = 3.0;
d_1 = 0;
d_2 = 0;
a = 1;

% wave packet parameters
omega = 1;
beta = 0.03;
n_0 = -150;
u_0 = 1;

% integration parameters
dt = 0.005;
t_max = 350;

fprintf("Input Omega: %.5f.\n\n",omega)
fprintf("Min Omega Chain 1: %.5f.\n", sqrt(d_1/m_1));
fprintf("Max Omega Chain 1: %.5f.\n\n", sqrt((4*c_1+d_1)/m_1));
fprintf("Min Omega Chain 2: %.5f.\n", sqrt(d_2/m_2));
fprintf("Max Omega Chain 2: %.5f.\n", sqrt((4*c_2+d_2)/m_2));


%% Initial Conditions

num = round((left_x:a:right_x)/a);
m = cat(2,m_1*ones(1,sum(num<0)),m_2*ones(1,sum(num>=0)));
c = cat(2,c_1*ones(1,sum(num<-1)),c_12*ones(1,sum(num==-1)),...
    c_2*ones(1,sum(num>-1)));
d = cat(2,d_1*ones(1,sum(num<0)),d_2*ones(1,sum(num>=0)));

k_1 = asin(sqrt(m.*(omega.^2-d./m)./(4*c)))*2/a;
g_1 = a/(2*omega)*sqrt((omega^2-d./m).*((4*c+d)./m-omega.^2));

disp=u_0*exp(-beta^2/2*(num-n_0).^2).*sin(num.*k_1*a);
disp(num>=-1)=0;

vel=-u_0*exp(-beta^2/2*(num-n_0).^2).*...
    (omega*cos(num.*k_1*a)-beta^2*g_1/a.*(num-n_0).*sin(num.*k_1*a));
vel(num>=-1)=0;


%% Solver

times = 0:dt:t_max;

result_disp = cell(1, fix(t_max/save_time)+1);
result_vel = cell(1, fix(t_max/save_time)+1);
result_e = cell(1, fix(t_max/save_time)+1);

for t=times
    if rem(t, save_time) == 0
        result_disp{fix(t/save_time)+1} = disp;
        result_vel{fix(t/save_time)+1} = vel;
        result_e{fix(t/save_time)+1} = energy(m,c,d,vel,disp);
    end
    acc1=(c./m).*(circshift(disp,-1)-disp)+(circshift(c,1)./m).*...
        (circshift(disp,1)-disp)-d./m.*disp;
    disp=disp+vel*dt+1/2*acc1*dt^2;
    acc2=(c./m).*(circshift(disp,-1)-disp)+(circshift(c,1)./m).*...
        (circshift(disp,1)-disp)-d./m.*disp;
    vel=vel+1/2*(acc1+acc2)*dt;
end


%% Plot Results

% energy (initial)
f1=figure(1); hold on
f1.Position = [0,50,1200,650];
plot(num,result_e{1},'Color','blue');
title('Распределение энергии (t=0)');
xlabel('Номер частицы по оси Ox');
ylabel('Энергия, усл.ед.');
grid on;
grid minor;
ylim([0 1])
xline(-1,'Color','Red')
hold off

% energy (solution)
descr_str = sprintf("\n (m_1=%.1f;   m_2=%.1f;   c_1=%.3f;   c_2=%.3f;"+...
    "   c_{12}=%.3f;   d_1=%.3f;   d_2=%.3f;   a=%.1f;   omega=%f)",...
    m_1,m_2,c_1,c_2,c_12,d_1,d_2,a,omega);

saved_times = (0:length(result_e))*save_time;
e_sum = cell2mat(cellfun(@(x) sum(x),result_e,'UniformOutput',false));
e_l = cell2mat(cellfun(@(x) sum(x(1:end,1:sum(num<0))),result_e,...
    'UniformOutput',false));
e_r = cell2mat(cellfun(@(x) sum(x(1:end,sum(num<=0):length(num))),...
    result_e,'UniformOutput',false));
e_l = e_l / max(e_sum);
e_r = e_r / max(e_sum);
e_sum = e_sum / max(e_sum);
f2=figure(2); hold on
f2.Position = [0,50,1400,650];
til = tiledlayout(1,2,'TileSpacing','Compact');
nexttile; hold on
s1 = plot(num,result_e{1},'Color','blue');
title('Распределение энергии');
xlabel('Номер частицы по оси Ox');
ylabel('Энергия, усл.ед.');
grid on;
grid minor;
ylim([0 1])
xline(-1,'Color','Red')
%ax = gca;
%ax.Clipping = "off";
hold off
nexttile; hold on
p1=plot(0,e_sum(1),'LineWidth',1.7,'Color','Black');
p2=plot(0,e_l(1),'LineWidth',1,'Color','Blue');
p3=plot(0,e_r(1),'LineWidth',1,'Color','Red');
ylim([min([e_sum, e_l, e_r]) 1.15*max([e_sum, e_l, e_r])])
title('Зависимость суммарных энергий в системе от времени');
xlabel('Время, усл.ед.');
ylabel('Энергия, усл.ед.');
legend('Полная энергия системы', 'Энергия левой цепочки',...
    'Энергия правой цепочки')
grid on;
grid minor;
hold off
til_title = title(til, "Энергия в момент времени t=0"+descr_str);
til_title.FontWeight = 'bold';
for i = 1:length(result_e)
    set(p1,'XData',saved_times(1:i));
    set(p1,'YData',e_sum(1:i));
    set(p2,'XData',saved_times(1:i));
    set(p2,'YData',e_l(1:i));
    set(p3,'XData',saved_times(1:i));
    set(p3,'YData',e_r(1:i));
    drawnow;
    s1.XData = num;
    s1.YData = result_e{i};
    title(til, "Энергия в момент времени t = "+string(saved_times(i))+...
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
    pause(anim_delay)
end
hold off


%% Energy Function

function e = energy(m,c,d,vel,disp)
    e = m / 2 .* vel.^2 + c / 4 .* (circshift(disp,-1)-disp).^2 + ...
        circshift(c,1) / 4 .* (circshift(disp,1)-disp).^2 + ...
        d / 2 .* disp.^2;
end
