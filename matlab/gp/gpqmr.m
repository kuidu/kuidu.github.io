function [x, y, varargout] = gpqmr(A, B, b, c, lambda, mu, varargin)
% GPBQMR will find approximate solution of following equation:
% 
%       [\lambda I    A] [x]   [b]
%       [              ] [ ] = [ ]
%       [B        \mu I] [y]   [c]
% 
%
% Syntaxes
% --------------
% [x, y] = gpqmr(A, B, b, c) 
% [x, y] = gpqmr(A, B, b, c, lambda, mu) 
% [x, y] = gpqmr(A, B, b, c, lambda, mu, tol) 
% [x, y] = gpqmr(A, B, b, c, lambda, mu, tol, maxit) 
% [x, y] = gpqmr(A, B, b, c, lambda, mu, tol, maxit)
% [x, y] = gpqmr(A, B, b, c, lambda, mu, tol, maxit, x0, y0)
% 
% [x, y, exitflag] = gpqmr(__)
% [x, y, exitflag, resvec] = gpqmr(__)
% 
% Parameters
% --------------
% A                 A m x n matrix 
%
% B                 A n x m matrix
%
% b, c              Right hand vectors
% 
% lambda, mu        coefficients of matrix
% 
% Name, Value       Name-Value pairs determine other options. Name should be in 
%                   {'tol', 'maxit','x0', 'y0}. 'tol' means
%                   tolerance for convergence, 'maxit' means maximum
%                   iterations, 'x0', 'y0' means initial gauss; Corresponding 
%                   default values are {1e-6, 100, 0, 0}.
%               
%
% Returns
% --------------
% x, y              approximate solution for the above system      
% 
% exitflag          covergence flag, `1` means failed, `0` means successed
% 
% resvec           vectors formed by norm of residuals
%
% Date
% --------------
% Kui Du, Jia-Jun Fan, and Fang Wang, 2024.1.4 


if nargin == 0
    help gpqmr; return;
end

% default vaules of input parameters
defaultTol = 1e-6;
defaultMaxit = 100; 
defaultX0 = [];
defaultY0 = [];

% check input parameters
checkMatrix = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
checkVector = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'column'});
checkPostive = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});
checkScalar = @(x) validateattributes(x, {'numeric'}, {'scalar'});

p = inputParser;
p.CaseSensitive = true;

addRequired(p, 'A', checkMatrix);
addRequired(p, 'B', checkMatrix);
addRequired(p, 'b', checkVector);
addRequired(p, 'c', checkVector);
addRequired(p, 'lambda', checkScalar);
addRequired(p, 'mu', checkScalar);

addOptional(p, 'tol', defaultTol, checkPostive);
addOptional(p, 'maxit', defaultMaxit, checkPostive);
addOptional(p, 'x0', defaultX0, checkVector);
addOptional(p, 'y0', defaultY0, checkVector);

parse(p, A, B, b, c, lambda, mu, varargin{:});
Parameters = p.Results;

tol = Parameters.tol;
maxit = Parameters.maxit;
x0 = Parameters.x0;
y0 = Parameters.y0;

[m, n] = size(A);

q = b; p = b; % p = randn(m, 1);
u = c; v = c; % v = randn(n, 1); 

pdotq = p' * q; 
udotv = u' * v;

eta1 = sqrt(abs(pdotq));
beta1 = pdotq / eta1;

delta1 = sqrt(abs(udotv));
gamma1 = udotv / delta1;

p1 = p / eta1; q1 = q / beta1;
u1 = u / delta1; v1 = v / gamma1;

p0 = zeros(m, 1); q0 = p0;
u0 = zeros(n, 1); v0 = u0;

if isempty(x0) && isempty(y0)
    x = p0; y = u0;
else
    b = b - (lambda * x0 + A * y0);
    c = c - (B * x0 + mu * y0);

    x = x0; y = y0;
end

fx_3 = p0; fx_2 = fx_3;
fx_1 = fx_3; fx0 = fx_3;

fy_3 = u0; fy_2 = fy_3;
fy_1 = fy_3; fy0 = fy_3;

xi_3 = 0; xi_2 = 0;
zeta_2 = 0;

varpi1bar = beta1; 
varpi2bar = delta1;

omega_1bar = 0; 
zeta_1bar = 0;
nu0bar = 0;

exitflag = 1;
resvec = zeros(maxit, 1);

for k = 1:maxit
    p = B' * v1 - delta1 * p0;
    q = A * u1 - gamma1 * q0;
    u = B * q1 - eta1 * u0;
    v = A' * p1 - beta1 * v0;

    alpha = p1' * q; theta = v1' * u;

    p = p - theta * p1;
    q = q - alpha * q1;
    u = u - theta * u1;
    v = v - alpha * v1;

    pdotq = p' * q;
    udotv = u' * v;

    eta2 = sqrt(abs(pdotq));
    beta2 = pdotq / eta2;

    delta2 = sqrt(abs(udotv));
    gamma2 = udotv / delta2;

    if k == 1
        rho1bar = lambda; thetabar = theta;
        nu1bar = alpha; rho2bar = mu;
        omega1bar = 0; nu2bar = eta2;
        zeta1bar = gamma2;

        xi_1 = 0; xi0 = 0;
        zeta0 = 0; zeta_1 = 0;
        omega_1 = 0; omega0 = 0;
        nu0 = 0;
    else
        % apply previous Givens rotations
        % first Givens rotation
        omega_1tilde = c1 * omega_1bar + s1 * theta;
        thetatilde = -s1 * omega_1bar + c1 * theta;
        zeta_1tilde = c1 * zeta_1bar + s1 * mu;
        rho2tilde = -s1 * zeta_1bar + c1 * mu;
        xi_1tilde = s1 * eta2;
        nu2tilde = c1 * eta2;

        % second Givens rotation
        omega_1 = c2 * omega_1tilde + s2 * nu0bar;
        nu0hat = -s2 * omega_1tilde + c2 * nu0bar;
        zeta_1 = c2 * zeta_1tilde;
        omega0hat = -s2 * zeta_1tilde;
        xi_1 = c2 * xi_1tilde;
        zeta0hat = -s2 * xi_1tilde;

        % third Givens rotation
        nu0v = c3 * nu0hat + s3 * thetatilde;
        thetabar = -s3 * nu0hat + c3 * thetatilde;
        omega0v = c3 * omega0hat + s3 * rho2tilde;
        rho2bar = -s3 * omega0hat + c3 * rho2tilde;
        zeta0v = c3 * zeta0hat + s3 * nu2tilde;
        nu2bar = -s3 * zeta0hat + c3 * nu2tilde;

        % forth Givens rotation
        nu0 = c4 * nu0v + s4 * lambda;
        rho1bar = -s4 * nu0v + c4 * lambda;
        omega0 = c4 * omega0v + s4 * alpha;
        nu1bar = -s4 * omega0v + c4 * alpha;
        zeta0 = c4 * zeta0v;
        omega1bar = -s4 * zeta0v;
        xi0 = s4 * gamma2;
        zeta1bar = c4 * gamma2;
    end

    % first Givens rotation
    [G, r] = planerot([rho1bar; delta2]);
    rho1tilde = r(1); c1 = G(1, 1); s1 = G(1, 2);
    nu1tilde = c1 * nu1bar;
    t1 = -s1 * nu1bar;
    
    varpi1tilde = c1 * varpi1bar;
    varpi4tilde = -s1 * varpi1bar;

    % second Givens rotation
    [G, r] = planerot([rho1tilde; thetabar]);
    rho1 = r(1); c2 = G(1, 1); s2 = G(1, 2);
    nu1 = c2 * nu1tilde + s2 * rho2bar;
    rho2hat = -s2 * nu1tilde + c2 * rho2bar;
    
    varpi1 = c2 * varpi1tilde + s2 * varpi2bar;
    varpi2hat = -s2 * varpi1tilde + c2 * varpi2bar;

    % third Givens rotation
    [G, r] = planerot([rho2hat; t1]);
    rho2v = r(1); c3 = G(1, 1); s3 = G(1, 2);

    varpi2v = c3 * varpi2hat + s3 * varpi4tilde;
    varpi4bar = -s3 * varpi2hat + c3 * varpi4tilde;

    % forth Givens rotation
    [G, r] = planerot([rho2v; beta2]);
    rho2 = r(1); c4 = G(1, 1); s4 = G(1, 2);

    varpi2 = c4 * varpi2v;
    varpi3bar = -s4 * varpi2v;
    
    fx1 = (q1 - xi_3 * fx_3 - zeta_2 * fx_2 - omega_1 * fx_1 - nu0 * fx0) / rho1;
    fy1 = -(xi_3 * fy_3 + zeta_2 * fy_2 + omega_1 * fy_1 + nu0 * fy0) / rho1;

    fx2 = -(xi_2 * fx_2 + zeta_1 * fx_1 + omega0 * fx0 + nu1 * fx1) / rho2;
    fy2 = (u1 - xi_2 * fy_2 - zeta_1 * fy_1 - omega0 * fy0 - nu1 * fy1) / rho2;

    dx = varpi1 * fx1 + varpi2 * fx2;
    dy = varpi1 * fy1 + varpi2 * fy2;

    x = x + dx; y = y + dy;

    rx = b - (lambda * x + A * y);
    ry = c - (B * x + mu * y);

    resvec(k) = hypot(norm(rx), norm(ry));

    if resvec(k) < tol
        exitflag = 0;
        break;
    end

    % update for next iteration
    eta1 = eta2;  beta1 = beta2;
    delta1 = delta2; gamma1 = gamma2;

    p0 = p1; p1 = p / eta2;
    q0 = q1; q1 = q / beta2;
    u0 = u1; u1 = u / delta2;
    v0 = v1; v1 = v / gamma2;

    varpi1bar = varpi3bar; 
    varpi2bar = varpi4bar;
    
    omega_1bar = omega1bar; 
    zeta_1bar = zeta1bar;
    nu0bar = nu2bar;

    xi_3 = xi_1; xi_2 = xi0;
    zeta_2 = zeta0;

    fx_3 = fx_1; fx_2 = fx0;
    fx_1 = fx1; fx0 = fx2;

    fy_3 = fy_1; fy_2 = fy0;
    fy_1 = fy1; fy0 = fy2;
end

resvec = resvec(1:k);

if nargout > 2
    varargout{1} = exitflag;
end

if nargout > 3
    varargout{2} = resvec;
end

end