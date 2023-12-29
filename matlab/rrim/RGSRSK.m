function [xall,tim] = RGSRSK(A,B,b,lambda,maxit)
% xall: computed approximate solutions, every m iterations compute once
% tim: computing time, every m iterations record once
% initial y and z are both zero vectors

[m,l] = size(A); 
[~,n] = size(B);
tim = zeros(maxit+1,1);
xall = zeros(n,maxit+1);
maxit = maxit*m;

y = zeros(l,1); r = b;
z = zeros(n,1); x = z;
it = 1;

tstart = tic;
ca = sum(A.^2); % pca = ra/sum(ca); % probability for A column chosen
J = randsample(l,maxit,true,ca);
rb = sum(B.^2,2); % prb = rb/sum(rb); % probability for B row chosen 
I = randsample(l,maxit,true,rb);

for k = 1:maxit
    j = J(k);
    i = I(k);
    d = A(:,j)'*r/ca(j); y(j) = y(j)+d; r = r - d*A(:,j);
    z = z - (B(i,:)*x-y(i))/rb(i)*B(i,:)';
    x = max(abs(z)-lambda,0).*sign(z);
    if mod(k,m) == 0
        it = it +1;
        xall(:,it) = x;
        tim(it) = toc(tstart); % computing time
    end
end
end