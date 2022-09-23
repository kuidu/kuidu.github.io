% test exercise 5 of Assignment 2

clear; clc; close all;
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[100 100 800 400]);

m = 5; n = 5; z = 10;
A = [zeros(n), zeros(n,m-n), diag(1:n); 
     zeros(m-n,m+n);
     diag(1:n), zeros(n,m-n), z*eye(n)]; 
% The first permutation 
p = [1:n, m+1:m+n, n+1:m];
A = A(p,p)
subplot(1,2,1); spy(A,'.',32);
% The second permutation 
p = 1:m+n; p(1:2:2*n-1) = 1:n; p(2:2:2*n) = n+(1:n); 
A = A(p,p)
subplot(1,2,2); spy(A,'.',32);


% figure; pcolor(A); axis ij; axis equal; axis([1 9 1 9]);
% figure; imagesc(A); axis equal;
