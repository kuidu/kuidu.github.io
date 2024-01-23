function[x, varargout] = minares(A, b,varargin)
% An implementation for `Minares`. Minares solve the linear system
%       Ax = b
%
% or the least squares problem
%       min || Ax - b ||_2
%
% minares finds the solution x which is belong to K(A,b).
%
%
% Syntaxes
% --------------
% x = minares(A, b)
% x = minares(A, b, tol)
% x = minares(A, b, tol, maxit)
% x = minares(A, b, tol, maxit, x0)
% x = minares(A, b, tol, maxit, x0, version)
% x = minares(A, b, tol, maxit, x0, version, criteria)
%
% [x, exitflag] = minares(__)
% [x, exitflag, resvec] = minares(__)
% [x, exitflag, resvec, Aresvec] = minares(__)
%
%
% Parameters
% ------------
% A                 m by m matrix
%
% b                 m by 1 right hand vector
%
% tol               stop tolerance
%
% maxit             number of maximum iterations
%
% criteria          stop criteria
%
% Name, Value       Name-Value pairs determine other options. Name should be in
%                   {'x0', 'version','criteria'}. 'x0' means initial guess; 'version'
%                   determines the implementation of rsmar, 1 means simpler
%                   version, i.e, start vector of Krylov subspace is Ar0, 2
%                   means classical version, i.e., start vector is r0.
%                   'criteria' means stop criteria, '1' means ||r||<tol, '2'means 
%                   ||Ar||<tol, '3' means ||r||<tol or ||Ar||<tol.
%                   Corresponding default values are {[], 2, 3}
%
%
% Returns
%-------------
% x                 approximate solution for the above system
%
% exitflag          covergency flag, `1` means failed, `0` means successed
%
% resvec            vectors formed by norm of residuals
%
% Aresvec           vectors formed by norm of Ar, r means residuals
%
% Kui Du, Jia-Jun Fan, and Fang Wang 2024.01.21
%

if nargin == 0
    help minares; return;
end
% default vaules of input parameters
defaultTol = 1e-6;
defaultMaxit = 100;
defaultX0 = [];
defaultVersion = 2;
defaultCriteria = 3;

% check input parameters
checkMatrix = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
checkVector = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'column'});
checkPostive = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
checkVersion = @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 1, '<=', 2});

p = inputParser;
p.CaseSensitive = true;

addRequired(p, 'A', checkMatrix);
addRequired(p, 'b', checkVector);
addOptional(p, 'tol', defaultTol, checkPostive);
addOptional(p, 'maxit', defaultMaxit, checkPostive);
addParameter(p, 'x0', defaultX0, checkVector);
addParameter(p, 'version', defaultVersion, checkVersion);
addOptional(p, 'criteria', defaultCriteria, checkPostive);

parse(p, A, b, varargin{:});
Parameters = p.Results;

tol = Parameters.tol;
maxit = Parameters.maxit;
x0 = Parameters.x0;

version = Parameters.version;
criteria = Parameters.criteria;

if version == 1
    minares_version = @minaresv1;
else
    minares_version = @minaresv2;
end


[x, exitflag, resvec, Aresvec] = minares_version(A, b, tol, maxit, x0, criteria);

if nargout > 1
    varargout{1} = exitflag;
end

if nargout > 2
    varargout{2} = resvec;
end

if nargout > 3
    varargout{3} = Aresvec;
end
end


function [x, exitflag, resvec, Aresvec] = minaresv1(A, b, tol, maxit, x0, criteria)
[n,~] = size(A);

if ~isempty(x0)
    r = b - A * x0;
    x=x0;
else
    x = zeros(n, 1);
    r = b;
end
v_k = A * r;
beta_k = norm(v_k);
v_k = v_k / beta_k;
w_k = r / beta_k;
t_tilde = beta_k;
v_km1 = zeros(n,1);
p_km2 = zeros(n,1);
p_km1 = zeros(n,1);
w_km1 = zeros(n,1);
c = -1;
s = 0;
lambda_tilde = 0;
eta_km2 = 0;
exitflag = 1;
resvec = zeros(maxit, 1); Aresvec = zeros(maxit, 1);
%%
for k=1 : maxit
    v_kp1 = A * v_k - beta_k * v_km1;
    alpha_k = v_k' * v_kp1;
    v_kp1 = v_kp1 - alpha_k * v_k;
    beta_kp1 = norm(v_kp1);
    v_kp1 = v_kp1 / beta_kp1;
    lambda_km1 = c * lambda_tilde + s * alpha_k;
    delta_tilde_k = s * lambda_tilde - c * alpha_k;
    eta_km1 = s * beta_kp1;
    lambda_tilde = -c * beta_kp1;
    delta_k = sqrt(delta_tilde_k^2 + beta_kp1^2);
    c = delta_tilde_k / delta_k;
    s = beta_kp1 / delta_k;
    t_hat_k = c * t_tilde;
    t_tilde = s * t_tilde;
    p_k = (w_k - eta_km2 * p_km2 - lambda_km1 * p_km1)/delta_k;
    x = x + t_hat_k * p_k;
    
    % Aresvec(k)=abs(t_tilde);
    r = b - A * x;

    Aresvec(k)=norm( A * r);
    resvec(k)=norm(r);
    if criteria==1
        if resvec(k) < tol
            exitflag = 0;
            break;
        end
    elseif criteria==2
        if  Aresvec(k) < tol
            exitflag = 0;
            break;
        end
    elseif  criteria==3
        if resvec(k) < tol || Aresvec(k) < tol
            exitflag = 0;
            break;
        end
    end
    w_kp1 = (v_k - alpha_k * w_k - beta_k * w_km1) / beta_kp1;

    % update for next iteraton
    w_km1 = w_k;
    w_k = w_kp1;
    v_km1 = v_k;
    v_k = v_kp1;
    eta_km2 = eta_km1;
    p_km2 = p_km1;
    p_km1 = p_k;
    beta_k = beta_kp1;

end
resvec = resvec(1:k);
Aresvec = Aresvec(1:k);
end


function [x, exitflag, resvec, Aresvec] = minaresv2(A, b, tol, maxit, x0,criteria)
m = size(A, 1);

if ~isempty(x0)
    r = b - A * x0;
    x=x0;
else
    x = zeros(m, 1);
    r = b;
end
beta_k=norm(r);
v_k=r/beta_k;
q_k=A*v_k;
alpha_k=v_k'*q_k;
q_k=q_k-alpha_k*v_k;
beta_kp1=norm(q_k);
v_kp1=q_k/beta_kp1;


zeta_bar=beta_k*beta_kp1;
zeta_bar2=beta_k*alpha_k;

% the short recurrence of ||r|| and ||Ar||
% pai_bar=0;
% pai_bar2=0;
% xi=0;
% theta_bar=0;
% tau_bar2=0;
% psi_bar=1;
% psi_bar2=1;
%
%
% chi_bar=beta_k;
%
% norm_r=chi_bar;
% norm_Ar=sqrt(zeta_bar2^2+zeta_bar^2);
%

resvec=zeros(maxit,1);
Aresvec=resvec;


exitflag = 1;



lambda_bar=alpha_k;
gamma_bar=beta_kp1;


w_kn2=zeros(m,1);
w_kn1=zeros(m,1);
d_kn2=zeros(m,1);
d_kn1=zeros(m,1);

elp_kn2=0;
elp_kn1=0;
gamma_kn1=0;

s_tilde_2kn4=0;
s_tilde_2kn3=0;
s_tilde_2kn2=0;

c_tilde_2kn4=-1;
c_tilde_2kn3=-1;
c_tilde_2kn2=-1;



for k = 1:maxit
    q_kp1=A*v_kp1-beta_kp1*v_k;
    alpha_kp1=v_kp1'*q_kp1;
    q_kp1=q_kp1-alpha_kp1*v_kp1;
    beta_kp2=norm(q_kp1);
    v_kp2=q_kp1/beta_kp2;



    % the first QR factorization
    lambda_k=sqrt(lambda_bar^2+beta_kp1^2);
    c_k=lambda_bar/lambda_k;
    s_k=beta_kp1/lambda_k;
    gamma_k=c_k*gamma_bar+s_k*alpha_kp1;
    elp_k=s_k*beta_kp2;
    lambda_bar=s_k*gamma_bar-c_k*alpha_kp1;
    gamma_bar=-c_k*beta_kp2;


    % the second QR factorization
    rho_kn2=s_tilde_2kn4*lambda_k;
    lambda_hat_k=-c_tilde_2kn4*lambda_k;
    phi_bar_kn1=s_tilde_2kn3*lambda_hat_k;
    phi_kn1=c_tilde_2kn2*phi_bar_kn1+s_tilde_2kn2*gamma_k;
    mu_bar_k=-c_tilde_2kn3*lambda_hat_k;
    gamma_hat_k=s_tilde_2kn2*phi_bar_kn1-c_tilde_2kn2*gamma_k;

    mu_bar2_k=sqrt(mu_bar_k^2+gamma_hat_k^2);
    c_tilde_2kn1=mu_bar_k/mu_bar2_k;
    s_tilde_2kn1=gamma_hat_k/mu_bar2_k;
    mu_k=sqrt(mu_bar2_k^2+elp_k^2);
    c_tilde_2k=mu_bar2_k/mu_k;
    s_tilde_2k=elp_k/mu_k;
    zeta_o_k=c_tilde_2kn1*zeta_bar2+s_tilde_2kn1*zeta_bar;
    zeta_bar2=s_tilde_2kn1*zeta_bar2-c_tilde_2kn1*zeta_bar;
    zeta_bar=s_tilde_2k*zeta_o_k;
    zeta_k=c_tilde_2k*zeta_o_k;

    % update for next iteration
    w_k=(v_k-gamma_kn1*w_kn1-elp_kn2*w_kn2)/lambda_k;
    d_k=(w_k-phi_kn1*d_kn1-rho_kn2*d_kn2)/mu_k;
    x=x+zeta_k*d_k;

    % the short recurrence of ||Ar||
    % norm_Ar=sqrt(zeta_bar2^2+zeta_bar^2);
    % Aresvec(k+1)=norm_Ar;

    r=b-A*x;
    Aresvec(k)=norm(A*r);

    % the short recurrence of ||r||

    % chi_k=c_k*chi_bar;
    % chi_bar=s_k*chi_bar;
    % psi_kn2=sqrt(psi_bar2^2+rho_kn2^2);
    % c_hat_2kn4=psi_bar2/psi_kn2;
    % s_hat_2kn4=rho_kn2/psi_kn2;
    % theta_kn2=c_hat_2kn4*theta_bar+s_hat_2kn4*phi_kn1;
    % ww_kn2=s_hat_2kn4*mu_k;
    % delta_k=s_hat_2kn4*theta_bar-c_hat_2kn4*phi_kn1;
    % eta_k=-c_hat_2kn4*mu_k;
    % tau_kn2=tau_bar2*psi_bar2/psi_kn2;
    % psi_bar2=sqrt(psi_bar^2+delta_k^2);
    % c_hat_2kn3=psi_bar/psi_bar2;
    % s_hat_2kn3=delta_k/psi_bar2;
    % theta_bar=s_hat_2kn3*eta_k;
    % psi_bar=-c_hat_2kn3*eta_k;
    % upsilon_k=s_hat_2kn4*pai_bar2-c_hat_2kn4*chi_k;
    % pai_bar2=c_hat_2kn3*pai_bar+s_hat_2kn3*upsilon_k;
    % pai_bar=s_hat_2kn3*pai_bar-c_hat_2kn3*upsilon_k;
    % tau_bar2=(xi-theta_kn2*tau_kn2)/psi_bar2;
    % xi=zeta_k-ww_kn2*tau_kn2;
    % tau_bar_k=(xi-theta_bar*tau_bar2)/psi_bar;
    % norm_r=sqrt((pai_bar2-tau_bar2)^2+(pai_bar-tau_bar_k)^2+chi_bar^2);
    %  resvec(k+1) = norm_r;
    %
    resvec(k) = norm(r);


    %  stop criteria
    if criteria==1
        if resvec(k) < tol
            exitflag = 0;
            break;
        end
    elseif criteria==2
        if  Aresvec(k) < tol
            exitflag = 0;
            break;
        end
    elseif  criteria==3
        if resvec(k) < tol || Aresvec(k) < tol
            exitflag = 0;
            break;
        end
    end


    % update for next iteration
    v_k=v_kp1;
    v_kp1=v_kp2;
    beta_kp1=beta_kp2;
    gamma_kn1=gamma_k;
    elp_kn2=elp_kn1;
    elp_kn1=elp_k;
    s_tilde_2kn4=s_tilde_2kn2;
    s_tilde_2kn3=s_tilde_2kn1;
    s_tilde_2kn2=s_tilde_2k;
    c_tilde_2kn4=c_tilde_2kn2;
    c_tilde_2kn3=c_tilde_2kn1;
    c_tilde_2kn2=c_tilde_2k;
    w_kn2=w_kn1;
    w_kn1=w_k;
    d_kn2=d_kn1;
    d_kn1=d_k;


end
resvec = resvec(1:k);
Aresvec = Aresvec(1:k);
end





