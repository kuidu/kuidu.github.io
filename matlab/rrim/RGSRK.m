function [xall,tim] = RGSRK(A,B,b,maxit)
% xall: computed approximate solutions, every m iterations compute once
% tim: computing time, every m iterations record once
% initial y and x are both zero vectors

% Kui Du, 2023.2.15 
[m,l] = size(A); 
[~,n] = size(B);
tim = zeros(maxit+1,1);
xall = zeros(n,maxit+1);
maxit = maxit*m;

y = zeros(l,1); r = b;
x = zeros(n,1); 
it = 1;

tstart = tic;
ra = sum(A.^2); % pra = ra/sum(ra); % probability for A column chosen
J = randsample(l,maxit,true,ra);
rb = sum(B.^2,2); % prb = rb/sum(rb); % probability for B row chosen 
I = randsample(l,maxit,true,rb);

for k = 1:maxit
    j = J(k);
    i = I(k);
    d = A(:,j)'*r/ra(j); y(j) = y(j)+d; r = r - d*A(:,j);
    x = x - (B(i,:)*x-y(i))/rb(i)*B(i,:)';
    if mod(k,m) == 0
        it = it +1;
        xall(:,it) = x;
        tim(it) = toc(tstart); % computing time
    end
end
end