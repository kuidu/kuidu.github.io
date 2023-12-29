function [xall,tim] = RKRK(A,B,b,maxit)
% xall: computed approximate solutions, every m iterations compute once
% tim: computing time, every m iterations record once
% initial y and x are both zero vectors

[m,l] = size(A); 
[~,n] = size(B);
tim = zeros(maxit+1,1);
xall = zeros(n,maxit+1);
maxit = maxit*m;

y = zeros(l,1);
x = zeros(n,1); 
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
    x = x - (B(i,:)*x-y(i))/rb(i)*B(i,:)';
    if mod(k,m) == 0
        it = it +1;
        xall(:,it) = x; 
        tim(it) = toc(tstart); % computing time
    end
end
end