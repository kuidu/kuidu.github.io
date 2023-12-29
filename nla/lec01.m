%% Unit circles for different norms: 1-,2-,infty-norm

clear; clc; close all;
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[100 100 600 400]);

x1 = linspace(-2,2,500);
x2 = linspace(-2,2,500);
[X,Y] = meshgrid(x1,x2);

Z = abs(X)+abs(Y);
contour(X,Y,Z,[1 1],'r-','LineWidth',2);
axis equal
hold on
axis([-2,2,-2,2])

Z = sqrt(X.^2+Y.^2);
contour(X,Y,Z,[1 1],'g-','LineWidth',2);

Z = max(abs(X),abs(Y));
contour(X,Y,Z,[1 1],'b-','LineWidth',2);

legend('$\|\cdot\|_1$','$\|\cdot\|_2$','$\|\cdot\|_\infty$','interpreter','latex')

title('Unit circles for different norms')

%% Figure for 1-norm
clear; clc; close all;

A = [1 2; 0 2];
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[100 100 1200 400]);
x1 = [-1:0.01:0, 0:0.01:1, 1:-0.01:0, 0:-0.01:-1];
y1 = [0:0.01:1, 1:-0.01:0, 0:-0.01:-1, -1:0.01:0];
x = A*[x1;y1];
subplot(1,2,1)
plot(x1,y1,'r-','LineWidth',2)
axis equal;
axis([-5,5,-5,5])

subplot(1,2,2)
plot(x(1,:),x(2,:),'r-','LineWidth',2)
axis equal;
axis([-5,5,-5,5])
hold on;
X1 = linspace(-5,5,500); Y2 = linspace(-5,5,500);
[X,Y] = meshgrid(X1,Y2); Z = abs(X)+abs(Y);
contour(X,Y,Z,[1,1],'k-','LineWidth',2,'ShowText','on');
contour(X,Y,Z,[2,2],'g-','LineWidth',2,'ShowText','on');
contour(X,Y,Z,[3,3],'m-','LineWidth',2,'ShowText','on');
contour(X,Y,Z,[4,4],'b-','LineWidth',2,'ShowText','on');
plot(2,2,'k.','MarkerSize',20)

%% Figure for 2-norm
clear; clc; close all;

A = [1 2; 0 2];
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[100 100 1200 400]);
theta = linspace(0,2*pi,500);
x1 = sin(theta);
y1 = cos(theta);
x = A*[x1;y1];
subplot(1,2,1)
plot(x1,y1,'g-','LineWidth',2)
axis equal;
axis([-3,3,-3,3])

subplot(1,2,2)
plot(x(1,:),x(2,:),'g-','LineWidth',2)
axis equal;
axis([-3,3,-3,3])
hold on;
X1 = linspace(-3,3,500); Y2 = linspace(-3,3,500);
[X,Y] = meshgrid(X1,Y2); Z = sqrt(X.^2+Y.^2);
contour(X,Y,Z,[2.9208,2.9208],'r-','LineWidth',2,'ShowText','on');
contour(X,Y,Z,[2*sqrt(2),2*sqrt(2)],'b-','LineWidth',2,'ShowText','on');
contour(X,Y,Z,[2,2],'k-','LineWidth',2,'ShowText','on');
plot(2,2,'k.','MarkerSize',20)


%% Figure for infty-norm
clear; clc; close all;

A = [1 2; 0 2];
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[100 100 1200 400]);
v = -1:0.01:1;
x1 = [v, ones(size(v)), -v, -ones(size(v))];
y1 = [ones(size(v)), -v, -ones(size(v)), v];
x = A*[x1;y1];
subplot(1,2,1)
plot(x1,y1,'b-','LineWidth',2)
axis equal;
axis([-4,4,-4,4])

subplot(1,2,2)
plot(x(1,:),x(2,:),'b-','LineWidth',2)
axis equal;
axis([-4,4,-4,4])
hold on;
X1 = linspace(-4,4,500); Y2 = linspace(-4,4,500);
[X,Y] = meshgrid(X1,Y2); Z = max(abs(X),abs(Y));
contour(X,Y,Z,[1,1],'r-','LineWidth',2,'ShowText','on');
contour(X,Y,Z,[2,2],'k-','LineWidth',2,'ShowText','on');
contour(X,Y,Z,[3,3],'g-','LineWidth',2,'ShowText','on');
plot(2,2,'k.','MarkerSize',20)
