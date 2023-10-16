%% Legendre polynomials
clear; clc;
n = 5;
x = sym('x'); 
P = sym('P',[n,1]);
P(1) = 1;
for j = 1:n-1
    P(j+1) = x^j;
    for i = 1:j
        P(j+1) = P(j+1)-int(P(j+1)*P(i),-1,1)/int(P(i)*P(i),-1,1)*P(i);
    end
    P(j+1) = P(j+1)/subs(P(j+1),x,1); 
end
P
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[380 320 800 400]);
fplot(P,[-1,1],'LineWidth',3)
xlabel('$x$','Interpreter','latex','FontSize',28)
axis([-1 1 -1.1 1.1])

%% Discrete Legendre polynomials
clear; clc;
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[380 320 800 400]);

x = (-128:128)'/128;
A = [x.^0 x.^1 x.^2 x.^3];
[Q,R] = qr(A,0);
scale = Q(257,:);
Q = Q*diag(1./scale);
plot(x,Q,'LineWidth',3)
xlabel('$x$','Interpreter','latex','FontSize',28)
axis([-1 1 -1.1 1.1])
