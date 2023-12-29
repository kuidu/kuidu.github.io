function [xall,tim] = GERKad(C,b,lambda,maxit)
% xall: computed approximate solutions, every m iterations compute once
% tim: computing time, every m iterations record once
% initial y = b and z = 0

% Kui Du, 2023.2.15 

[m,n] = size(C); 
xall = zeros(n,maxit+1);
tim = zeros(maxit+1,1);
maxit = maxit*m;
y = b;
z = zeros(n,1); x = z;
it = 1;

tstart = tic;
C2 = C.^2;
ra = sum(C2);
J = randsample(n,maxit,true,ra);
rb = sum(C2,2);
I = randsample(m,maxit,true,rb);

for k = 1:maxit
    j = J(k);
    i = I(k);
    y = y - C(:,j)'*y/ra(j)*C(:,j);
    z = z - (C(i,:)*x-b(i)+y(i))/rb(i)*C(i,:)';
    x = max(abs(z)-lambda,0).*sign(z);
    if mod(k,m) == 0
        it = it +1;
        xall(:,it) = x; 
        tim(it) = toc(tstart); % computing time
    end
end
end