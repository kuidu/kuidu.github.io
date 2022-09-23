% 
%% clgs vs mgs vs matlab qr
clear; clc; % format compact; 
% A = [0.7 0.70711; 0.70001 0.70711];
% [QC,~] = clgs(A); [QM,~] = mgs(A); [QQ,~] = qr(A);
% norm(QC'*QC-eye(2)) 
% norm(QM'*QM-eye(2)) 
% norm(QQ'*QQ-eye(2))

m = 512; 
n = 256;
A = rand(m,n);
[QC,RC] = clgs(A); 
ClGS_error=norm(QC'*QC-eye(n)) 
[QM,RM] = mgs(A); 
MGS_error=norm(QM'*QM-eye(n)) 
% [QQ,RR] = qr(A,0);
% norm(QQ'*QQ-eye(n))

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
%xlabel('$x$','Interpreter','latex','FontSize',28)
axis([-1 1 -1.1 1.1])