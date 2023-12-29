function [xall,tim] = RKRSK(A,B,b,lambda,maxit)
% xall: computed approximate solutions, every m iterations compute once
% tim: computing time, every m iterations record once
% initial y and z are both zero vectors

% Kui Du, 2023.2.15 

[m,l] = size(A); 
[~,n] = size(B);
tim = zeros(maxit+1,1);
xall = zeros(n,maxit+1);
maxit = maxit*m;

y = zeros(l,1);
z = zeros(n,1); x = z;
it = 1;

tstart = tic;
ra = sum(A.^2,2); % pra = ra/sum(ra); % probability for A row chosen
J = randsample(m,maxit,true,ra);
rb = sum(B.^2,2); % prb = rb/sum(rb); % probability for B row chosen 
I = randsample(l,maxit,true,rb);

for k = 1:maxit
    j = J(k);
    i = I(k);
    y = y - (A(j,:)*y-b(j))/ra(j)*A(j,:)';
    z = z - (B(i,:)*x-y(i))/rb(i)*B(i,:)';
    x = max(abs(z)-lambda,0).*sign(z);
    if mod(k,m) == 0
        it = it +1;
        xall(:,it) = x; 
        tim(it) = toc(tstart); % computing time
    end
end
end