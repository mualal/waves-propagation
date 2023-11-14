%% Lattice dispersion relation

c = 1;
a = 1;
m = 1;

k_x = -2*pi/a:0.01:2*pi/a;
k_y = -2*pi/a:0.01:2*pi/a;
[X, Y] = meshgrid(k_x,k_y);

Z = 2*sqrt(c/m)*sqrt(sin(X*a/2).^2+sin(Y*a/2).^2);

contourf(X,Y,Z,25)
title("Карта изолиний частоты, полученная" +...
    " из дисперсионного соотношения"+...
    sprintf("\n (c=%.2f;   a=%.2f;   m=%.2f)",c,a,m),...
    "FontSize",14)
xlabel('$k_x$','Interpreter',"latex",'FontSize',16)
ylabel('$k_y$','Interpreter',"latex",'FontSize',16)
set(gca,'XTick',-6:1:6)
set(gca,'YTick',-6:1:6)
colorbar
axis equal

%% 
