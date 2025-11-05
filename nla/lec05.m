%% band matrix LU factorization

n=20,lbw=4,ubw=3; A = tril(triu(randn(n,n),-lbw),ubw);
Ad = A + 100*eye(n); % no pivoting required
[Ld,Ud,Pd]=lu(Ad);
figure(1),
subplot(221),spy(Ad),title('Ad'),subplot(222),spy(Ld),title('Ld'),
subplot(223),spy(Ud),title('Ud'),subplot(224),spy(Pd),title('Pd')
figure(2),
[L,U,P]=lu(A);
subplot(221),spy(A),title('A'),subplot(222),spy(L),title('L'),
subplot(223),spy(U),title('U'),subplot(224),spy(P),title('P')

%% sparse matrix

A = eye(8);A(1,2:8)=.1;A(2:8,1)=.1;
[L,U]=lu(A);A,L,U,
figure(1), clf, spy(A,'k'), pause
spy(L,'b'), hold on, spy(U,'r'), pause
spy(A,'k')

Arev=A(8:-1:1,8:-1:1);
[Lrev,Urev]=lu(Arev); Arev, Lrev, Urev
figure(2), clf, spy(Arev,'k'), pause
spy(Lrev,'b'), hold on, spy(Urev,'r'),pause
spy(Arev,'k')

%%
help sparse