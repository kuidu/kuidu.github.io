%%
 load clown.mat, [U,S,V]=svd(X);
   figure(1), clf, image(X), colormap('gray')
   figure(2), clf, for k=[1:10,20:10:200],
     Xk=U(:,1:k)*S(1:k,1:k)*(V(:,1:k))'; ...
     err = norm(X-Xk)/norm(X); compr = k*(200+320+1)/(200*320); ...
     figure(2), image(Xk), colormap('gray'), ...
     title(['k= ',int2str(k),' err= ', num2str(err),' compression=', ...
       num2str(compr)]), ...
    pause, end

%%
figure(3), s = diag(S);
semilogy(1:200,s/s(1),'r',1:200,(1:200)*(200+320+1)/(200*320),'b'),
title('Error in red, compression in blue'), xlabel('rank'), grid

%% test exercise 5 of Assignment 2 

clear; clc; close all;
set(0, 'defaultaxeslinewidth',  1);
set(0, 'defaultaxesfontsize',   16);
figure('Position',[100 100 1200 400]);

m = 12; n = 6; z = 10; % this version requires m>=n
A = [zeros(n), zeros(n,m-n), diag(1:n); 
     zeros(m-n,m+n);
     diag(1:n), zeros(n,m-n), z*eye(n)]; 
subplot(1,3,1); spy(A,'.',32); 
% The first permutation 
p = [1:n, m+1:m+n, n+1:m];
A = A(p,p)
subplot(1,3,2); spy(A,'.',32);
% The second permutation 
p = 1:m+n; p(1:2:2*n-1) = 1:n; p(2:2:2*n) = n+(1:n); 
A = A(p,p)
subplot(1,3,3); spy(A,'.',32);
