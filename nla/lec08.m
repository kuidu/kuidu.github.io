%% Assignment 8 -- San Zhang

%% Programming 1
clear; clc; close all;
L=diag([1 2 6 30]); 
S=randn(4); 
A=S*L*inv(S);
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   12);

% power iteration
v = randn(4,1); 
v = v/norm(v);
for k=1:10
    w=A*v;
    v = w/norm(w);
    lambda(k)=v'*A*v;
end        
figure('Position',[100 100 600 500]);
cr1=(6/30).^(1:10);
semilogy(1:10,cr1,'--',1:10,abs(lambda-30),'-','LineWidth',2)
legend({'$(6/30)^k$',...
        '$|\lambda^{(k)}-30|$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of power iteration','Interpreter','latex',...
      'FontSize',20)
  
% inverse iteration
v = randn(4,1); 
v = v/norm(v);
for k=1:10
    w=A\v;
    v = w/norm(w);
    lambda(k)=v'*A*v;
end   
figure('Position',[100 100 600 500]);
cr2=(1/2).^(1:10);
semilogy(1:10,cr2,'--',1:10,abs(lambda-1),'-','LineWidth',2)
legend({'$(1/2)^k$',...
        '$|\lambda^{(k)}-1|$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of inverse iteration','Interpreter','latex',...
      'FontSize',20)
  
  
% rayleigh quotient iteration
v = randn(4,1); 
v = v/norm(v);
lambda(1) = v'*A*v;
for k=1:10
    w=(A-lambda(k)*eye(4))\v;
    v = w/norm(w);
    lambda(k+1)=v'*A*v;
end
[~,i]=min(abs([lambda(11)-1 lambda(11)-2 lambda(11)-6 lambda(11)-30]));
lam=[1 2 6 30];
error = abs(lambda(2:11)-lam(i));
figure('Position',[100 100 600 500]);
semilogy(1:10,error,'-','LineWidth',2)
legend({['$|\lambda^{(k)}-$',num2str(lam(i),'%10.0f'),'$|$']},'Interpreter','latex','FontSize',18);
title('Convergence of rayleigh quotient iteration','Interpreter','latex',...
      'FontSize',20)