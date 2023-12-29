%% Convergence of subspace iteration
clear; clc; format short; format compact;
L=diag([1 2 6 30]);
S=randn(4); %condnumberS=cond(S)
A=S*L/S;
Q=eye(4);
lambda_1cr=zeros(19,1);
lambda_2cr=zeros(19,1);
lambda_3cr=zeros(19,1);
lambda_4cr=zeros(19,1);
Ablock_1=zeros(19,1);
Ablock_2=zeros(19,1);
Ablock_3=zeros(19,1);
for k=1:19
    Z=A*Q;
    [Q,R]=qr(Z);
    Ak=Q'*A*Q;
    lambda_1cr(k)=abs(Ak(1,1)-30);
    lambda_2cr(k)=abs(Ak(2,2)-6);
    lambda_3cr(k)=abs(Ak(3,3)-2);
    lambda_4cr(k)=abs(Ak(4,4)-1);
    Ablock_1(k)=norm(Ak(2:4,1:1));
    Ablock_2(k)=norm(Ak(3:4,1:2));
    Ablock_3(k)=norm(Ak(4:4,1:3));
end
A19=Ak
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[380 320 1000 400]);
subplot(1,2,1)
cr1=(6/30).^(1:19);
semilogy(1:19,cr1,'--',1:19,lambda_1cr,'-','LineWidth',2)
legend({'$|\lambda_2/\lambda_1|^k=(6/30)^k$',...
        '$|{\bf A}^{(k)}_{11}-\lambda_1|$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of $\lambda_1=30$','Interpreter','latex',...
      'FontSize',20)
subplot(1,2,2)
cr2=(2/6).^(1:19);
semilogy(1:19,cr2,'--',1:19,lambda_2cr,'-','LineWidth',2)
legend(...
{'$\max(|\lambda_3/\lambda_2|^k,|\lambda_2/\lambda_1|^k)=(2/6)^k$',...
    '$|{\bf A}^{(k)}_{22}-\lambda_2|$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of $\lambda_2=6$','Interpreter','latex','FontSize',20)
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[380 320 1000 400]);
subplot(1,2,1)
cr3=(1/2).^(1:19);
semilogy(1:19,cr3,'--',1:19,lambda_3cr,'-','LineWidth',2)
legend(...
{'$\max(|\lambda_4/\lambda_3|^k,|\lambda_3/\lambda_2|^k)=(1/2)^k$',...
    '$|{\bf A}^{(k)}_{33}-\lambda_3|$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of $\lambda_3=2$','Interpreter','latex','FontSize',20)
subplot(1,2,2)
cr4=(1/2).^(1:19);
semilogy(1:19,cr4,'--',1:19,lambda_4cr,'-','LineWidth',2)
legend({'$|\lambda_4/\lambda_3|^k=(1/2)^k$',...
    '$|{\bf A}^{(k)}_{44}-\lambda_4|$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of $\lambda_4=1$','Interpreter','latex','FontSize',20)
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[380 320 1000 400]);
subplot(1,2,1)
Acr1=(6/30).^(1:19);
semilogy(1:19,Acr1,'--',1:19,Ablock_1,'-','LineWidth',2)
legend({'$|\lambda_2/\lambda_1|^k=(6/30)^k$',...
    '$\|{\bf A}^{(k)}(2:4,1)\|_2$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of ${\bf A}^{(k)}(2:4,1)$',...
    'Interpreter','latex','FontSize',20)
subplot(1,2,2)
Acr2=(2/6).^(1:19);
semilogy(1:19,Acr2,'--',1:19,Ablock_2,'-','LineWidth',2)
legend({'$|\lambda_3/\lambda_2|^k=(2/6)^k$',...
    '$\|{\bf A}^{(k)}(3:4,1:2)\|_2$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of ${\bf A}^{(k)}(3:4,1:2)$','Interpreter',...
    'latex','FontSize',20)
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[380 320 1000 400]);
subplot(1,2,1)
Acr3=(1/2).^(1:19);
semilogy(1:19,Acr3,'--',1:19,Ablock_3,'-','LineWidth',2)
legend({'$|\lambda_4/\lambda_3|^k=(1/2)^k$',...
    '$\|{\bf A}^{(k)}(4,1:3)\|_2$'},...
        'Location','Best','Interpreter','latex','FontSize',18);
title('Convergence of ${\bf A}^{(k)}(4,1:3)$','Interpreter',...
    'latex','FontSize',20)

%% Convergence of QR with Rayleigh quotient shift
clear; clc; format short e; format compact;
L=diag([1 2 6 30]);
S=randn(4); condnumberS=cond(S)
A=S*L/S;
for k=1:10
    [Q,R]=qr(A-A(4,4)*eye(4));
    A=R*Q+A(4,4)*eye(4);
    A4r=A(4,:)
end
A10=A
