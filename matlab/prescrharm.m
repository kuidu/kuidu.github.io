function [H,r0e1,U,C]=prescrharm(betas,rnrms)
% betas:  An upper triangular matrix with the kth column containing the
% coefficients of the monic polynomial whose roots are the prescribed
% harmonic Ritz values for the kth iteration
% rnrms:   The sequence of prescribed GMRES residual norms
% Output: When GMRES is applied to (H,r0e1) it generates the prescribed
% residual norms and harmonic Ritz values. Also, H=U^(-1)*C*U with C the
% companion matrix for H and U upper triangular with a positive real
% diagonal.
%
% Jurjen Duintjer Tebbens, 10.2.2017. 
n=length(betas(:,1));
U=eye(n);
 U(1,1)=1/rnrms(1);
for k=1:n-1
    e1=zeros(k,1);e1(1)=1;
    v=U(1:k,1:k)*U(1:k,1:k)'*e1;
    delt=angle(betas(1,k));
    if rnrms(k+1)~=rnrms(k)
       m = sqrt(1/rnrms(k+1)^2-1/rnrms(k)^2);
       ud = (1/m)*1/(abs(betas(1,k))*rnrms(k+1)^2);
       U(k+1,k+1)=ud;
       U(1,k+1)=m*exp(delt*sqrt(-1));
       if k > 1
          U(2:k,k+1)=(betas(2:k,k)-v(2:k)/(conj(U(1,k+1)*ud)))*ud;
       end
    else
        U(1,k+1)=0;
        if k > 1
           U(2:k,k+1)=randn(k-1,1)+i*randn(k-1,1);
        end
        U(k+1,k+1)=abs(randn(1,1))+eps;
    end
end
C=zeros(n,n);
C(:,n)=betas(:,n);
C(2:n,1:n-1)=eye(n-1);
H=U^(-1)*C*U;
r0e1=zeros(n,1);
r0e1(1)=rnrms(1);
