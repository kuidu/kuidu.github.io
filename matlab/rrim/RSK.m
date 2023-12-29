function [xall,tim] = RSK(B,y,lambda,maxit)
% xall: computed approximate solutions, every l iterations compute once
% tim: computing time, every l iterations record once
% initial z is zero vector

[l,n] = size(B);
tim = zeros(maxit+1,1);
xall = zeros(n,maxit+1);
maxit = maxit*l;
z = zeros(n,1);
x = max(abs(z)-lambda,0).*sign(z);
it = 1;

tstart = tic;
rb = sum(B.^2,2); % prb = rb/sum(rb); % probability for B row chosen 
I = randsample(l,maxit,true,rb);

for k = 1:maxit
    i = I(k);
    z = z - (B(i,:)*x-y(i))/rb(i)*B(i,:)';
    x = max(abs(z)-lambda,0).*sign(z);
    if mod(k,l) == 0
        it = it +1;
        xall(:,it) = x;
        tim(it) = toc(tstart); % computing time
    end
end
end