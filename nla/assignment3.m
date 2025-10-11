%% Assignment 3 -- San Zhang

%% Programming 1
%
% <include>mgs.m</include>
%

%% Programming 2 

% Part (i)
clear; clc;
x=sym('x'); w=1+x^2;
P=sym('P',[6,1]);
P(1)=1;
for j=1:5
    P(j+1)=x^j;
    for i=1:j
        P(j+1)=P(j+1)-int(w*P(j+1)*P(i),-1,1)/int(w*P(i)*P(i),-1,1)*P(i);
    end
    P(j+1)=P(j+1)/subs(P(j+1),x,1); 
end
P
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[380 320 800 400]);
fplot(P,[-1,1],'LineWidth',3)
xlabel('$x$','Interpreter','latex','FontSize',28)
axis([-1 1 -1.1 1.1])
%%
% Part (ii)
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[380 320 800 400]);

x = (-128:128)'/128;
A = [x.^0 x.^1 x.^2 x.^3 x.^4 x.^5];
% W = (1+x.^2).^(1/2)*ones(1,6);
% [Q,R] = qr(A.*W,0);
% Q = Q./W;
[Q,R] = qr(diag(sqrt(1+x.^2))*A,0);
Q = diag(1./sqrt(1+x.^2))*Q;
scale = Q(257,:);
Q = Q*diag(1./scale);
plot(x,Q,'LineWidth',3)
xlabel('$x$','Interpreter','latex','FontSize',28)
axis([-1 1 -1.1 1.1])

%%
% close all;


% %%
% % Part (i)
% clear; clc;
% x=sym('x'); w=1;
% P=sym('P',[3,1]);
% P(1)=1;
% for j=1:2
%     P(j+1)=x^j;
%     for i=1:j
%         P(j+1)=P(j+1)-int(w*P(j+1)*P(i),0,1)/int(w*P(i)*P(i),0,1)*P(i);
%     end
%     P(j+1)=P(j+1)/subs(P(j+1),x,1); 
% end
% P
% set(0, 'defaultaxeslinewidth',  1);
% set(0, 'defaultaxesfontsize',   16);
% figure('Position',[380 320 800 400]);
% fplot(P,[-1,1],'LineWidth',3)
% xlabel('$x$','Interpreter','latex','FontSize',28)
% axis([-1 1 -1.1 1.1])
