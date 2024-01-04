function [xl, yl, varargout] = gpbilq(A, B, b, c, lambda, mu, varargin)
% GPBILQ will find approximate solution of following equation:
% 
%       [\lambda I    A] [x]   [b]
%       [              ] [ ] = [ ]
%       [B        \mu I] [y]   [c]
% 
%
% Syntaxes
% --------------
% [x, y] = gpbilq(A, B, b, c) 
% [x, y] = gpbilq(A, B, b, c, lambda, mu) 
% [x, y] = gpbilq(A, B, b, c, lambda, mu, tol) 
% [x, y] = gpbilq(A, B, b, c, lambda, mu, tol, maxit) 
% [x, y] = gpbilq(A, B, b, c, lambda, mu, tol, maxit)
% [x, y] = gpbilq(A, B, b, c, lambda, mu, tol, maxit, x0, y0)
% 
% [x, y, exitflag] = gpbilq(__)
% [x, y, exitflag, reslvec] = gpbilq(__)
% [x, y, exitflag, reslvec, rescvec] = gpbilq(__)
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
% xl, yl            approximate solution for the above system obtained by
%                   GPBILQ
% 
% exitflag          covergence flag, `1` means failed, `0` means successed
% 
% reslvec           vectors formed by norm of residuals related to GPBILQ
%
% rescvec           vectors formed by norm of residuals related to GPBICG
%
% Date
% --------------
% Kui Du, Jia-Jun Fan, and Fang Wang, 2024.1.4 


if nargin == 0
    help gpbilq; return;
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

if isempty(x0) && isempty(y0)
    xl = zeros(m, 1); 
    yl = zeros(n, 1);
else
    b = b - (lambda * x0 + A * y0);
    c = c - (B * x0 + mu * y0);
    xl = x0; yl = y0;
end

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

fx1tilde = q1; fx2tilde = xl;
fy1tilde = yl; fy2tilde = u1;

Au = A * u1; Bq = B * q1;
alpha = p1' * Au; theta = v1' * Bq;

p = B' * v1 - theta * p1;
q = Au - alpha * q1;
u = Bq - theta * u1;
v = A' * p1 - alpha * v1;

pdotq = p' * q; 
udotv = u' * v;

eta2 = sqrt(abs(pdotq));
beta2 = pdotq / eta2;

delta2 = sqrt(abs(udotv));
gamma2 = udotv / delta2;

p2 = p / eta2; q2 = q / beta2;
u2 = u / delta2; v2 = v / gamma2;

rho1bar = lambda; alphabar = alpha;
nu2bar = theta; rho2bar = mu;
omega3bar = 0; nu3bar = beta2;
zeta4bar = delta2;

varpi_3 = 0; varpi_2 = 0;
varpi_1 = 0; varpi0 = 0;

nu1 = 0; omega1 = 0; omega2 = 0;
zeta1 = 0; zeta2 = 0; zeta3 = 0;
xi1 = 0; xi2 = 0; xi3 = 0; xi4 = 0;

exitflag = 1;
reslvec = zeros(maxit+1, 1);
rescvec = zeros(maxit+1, 1);

reslvec(1) = hypot(norm(b), norm(c));

if abs(rho1bar * rho2bar - alphabar * nu2bar) > 1e-14
    [G, y] = planerot([rho1bar; alphabar]);
    rho1d = y(1); cc = G(1, 1); ss = G(1, 2);
    nu2d = cc * nu2bar + ss * rho2bar;
    rho2d = -ss * nu2bar + cc * rho2bar;

    varpi3 = beta1 / rho1d;
    varpi4 = (delta1 - nu2d * varpi3) / rho2d;

    aa = cc * varpi3 - ss * varpi4;
    bb = ss * varpi3 + cc * varpi4;

    dxc = aa * fx1tilde + bb * fx2tilde;
    dyc = aa * fy1tilde + bb * fy2tilde;

    xc = xl + dxc; yc = yl + dyc;

    rxc = b - (lambda * xc + A * yc);
    ryc = c - (B * xc + mu * yc);

    rescvec(1) = hypot(norm(rxc), norm(ryc));
else
    rescvec(1) = NaN;
end

for k = 2:maxit
    p = B' * v2 - delta2 * p1;
    q = A * u2 - gamma2 * q1;
    u = B * q2 - eta2 * u1;
    v = A' * p2 - beta2 * v1;

    alpha = p2' * q; theta = v2' * u;

    p = p - theta * p2;
    q = q - alpha * q2;
    u = u - theta * u2;
    v = v - alpha * v2;

    pdotq = p' * q;
    udotv = u' * v;

    eta3 = sqrt(abs(pdotq));
    beta3 = pdotq / eta3;

    delta3 = sqrt(abs(udotv));
    gamma3 = udotv / delta3;
    
    % first Givens rotation
    [G, y] = planerot([rho1bar; gamma2]);
    rho1tilde = y(1); c1 = G(1, 1); s1 = G(1, 2);
    nu2tilde = c1 * nu2bar; t1 = -s1 * nu2bar;
    omega3tilde = c1 * omega3bar + s1 * alpha;
    alphatilde = -s1 * omega3bar + c1 * alpha;
    zeta4tilde = c1 * zeta4bar + s1 * mu;
    rho4tilde = -s1 * zeta4bar + c1 * mu;
    xi5tilde = s1 * beta3;
    nu5tilde = c1 * beta3;

    fx1hat = c1 * fx1tilde; 
    fx4hat = -s1 * fx1tilde;
    fy1hat = c1 * fy1tilde + s1 * u2;
    fy4hat = -s1 * fy1tilde + c1 * u2;

    % second Givens rotation
    [G, y] = planerot([rho1tilde; alphabar]);
    rho1 = y(1); c2 = G(1, 1); s2 = G(1, 2);
    nu2 = c2 * nu2tilde + s2 * rho2bar;
    rho2hat = -s2 * nu2tilde + c2 * rho2bar;
    omega3 = c2 * omega3tilde + s2 * nu3bar;
    nu3hat = -s2 * omega3tilde + c2 * nu3bar;
    zeta4 = c2 * zeta4tilde;
    omega4hat = -s2 * zeta4tilde;
    xi5 = c2 * xi5tilde;
    zeta5hat = -s2 * xi5tilde;

    fx1 = c2 * fx1hat + s2 * fx2tilde;
    fx2hat = -s2 * fx1hat + c2 * fx2tilde;
    fy1 = c2 * fy1hat + s2 * fy2tilde;
    fy2hat = -s2 * fy1hat + c2 * fy2tilde;

    % third Givens rotation
    [G, y] = planerot([rho2hat; t1]);
    rho2v = y(1); c3 = G(1, 1); s3 = G(1, 2);
    nu3v = c3 * nu3hat + s3 * alphatilde;
    alphabar = -s3 * nu3hat + c3 * alphatilde;
    omega4v = c3 * omega4hat + s3 * rho4tilde;
    rho4bar = -s3 * omega4hat + c3 * rho4tilde;
    zeta5v = c3 * zeta5hat + s3 * nu5tilde;
    nu5bar = -s3 * zeta5hat + c3 * nu5tilde;

    fx2v = c3 * fx2hat + s3 * fx4hat;
    fx4tilde = -s3 * fx2hat + c3 * fx4hat;
    fy2v = c3 * fy2hat + s3 * fy4hat;
    fy4tilde = -s3 * fy2hat + c3 * fy4hat;

    % forth Givens rotation
    [G, y] = planerot([rho2v; eta2]);
    rho2 = y(1); c4 = G(1, 1); s4 = G(1, 2);
    nu3 = c4 * nu3v + s4 * lambda;
    rho3bar = -s4 * nu3v + c4 * lambda;
    omega4 = c4 * omega4v + s4 * theta;
    nu4bar = -s4 * omega4v + c4 * theta;
    zeta5 = c4 * zeta5v;
    omega5bar = -s4 * zeta5v;
    xi6 = s4 * delta3; 
    zeta6bar = c4 * delta3;

    fx2 = c4 * fx2v + s4 * q2;
    fx3tilde = -s4 * fx2v + c4 * q2;
    fy2 = c4 * fy2v; 
    fy3tilde = -s4 * fy2v;

    % update solution
    if k == 2
        varpi1 = beta1 / rho1;
        varpi2 = (delta1 - nu2 * varpi1) / rho2;
    elseif k == 3
        varpi1 = -(omega1 * varpi_1 + nu1 * varpi0) / rho1;
        varpi2 = -(zeta2 * varpi_1 + omega2 * varpi0 + nu2 * varpi1) / rho2;
    else
        varpi1 = -(xi1 * varpi_3 + zeta1 * varpi_2 + omega1 * varpi_1 + nu1 * varpi0) / rho1;
        varpi2 = -(xi2 * varpi_2 + zeta2 * varpi_1 + omega2 * varpi0 + nu2 * varpi1) / rho2;
    end

    dxl = varpi1 * fx1 + varpi2 * fx2;
    dyl = varpi1 * fy1 + varpi2 * fy2;

    xl = xl + dxl; yl = yl + dyl;

    rxl = b - (lambda * xl + A * yl);
    ryl = c - (B * xl + mu * yl);

    reslvec(k) = hypot(norm(rxl), norm(ryl));

    if abs(rho3bar * rho4bar - alphabar * nu4bar) > 1e-14
        [G, y] = planerot([rho3bar; alphabar]);
        rho3d = y(1); cc = G(1, 1); ss = G(1, 2);
        nu4d = cc * nu4bar + ss * rho4bar;
        rho4d = -ss * nu4bar + cc * rho4bar;

        varpi3 = -(xi3 * varpi_1 + zeta3 * varpi0 + omega3 * varpi1 + nu3 * varpi2) / rho3d;
        varpi4 = -(xi4 * varpi0 + zeta4 * varpi1 + omega4 * varpi2 + nu4d * varpi3) / rho4d;

        aa = cc * varpi3 - ss * varpi4; 
        bb = ss * varpi3 + cc * varpi4;

        dxc = aa * fx3tilde + bb * fx4tilde;
        dyc = aa * fy3tilde + bb * fy4tilde;

        xc = xl + dxc; yc = yl + dyc;

        rxc = b - (lambda * xc + A * yc);
        ryc = c - (B * xc + mu * yc);

        rescvec(k) = hypot(norm(rxc), norm(ryc));
    else
        rescvec(k) = NaN;
    end

    if reslvec(k) < tol
        exitflag = 0;
        break;
    end

    % update for next iteration
    p1 = p2; p2 = p / eta3;
    q1 = q2; q2 = q / beta3;
    u1 = u2; u2 = u / delta3;
    v1 = v2; v2 = v / gamma3;

    fx1tilde = fx3tilde; fx2tilde = fx4tilde;
    fy1tilde = fy3tilde; fy2tilde = fy4tilde;

    delta2 = delta3; gamma2 = gamma3;
    eta2 = eta3; beta2 = beta3;

    rho1bar = rho3bar; nu2bar = nu4bar;
    rho2bar = rho4bar; omega3bar = omega5bar;
    nu3bar = nu5bar; zeta4bar = zeta6bar;

    varpi_3 = varpi_1; varpi_2 = varpi0;
    varpi_1 = varpi1; varpi0 = varpi2;

    nu1 = nu3; omega1 = omega3; 
    omega2 = omega4; zeta1 = zeta3; 
    zeta2 = zeta4; zeta3 = zeta5;
    xi1 = xi3; xi2 = xi4; 
    xi3 = xi5; xi4 = xi6;
end

reslvec = reslvec(1:k);
rescvec = rescvec(1:k);

if nargout > 2
    varargout{1} = exitflag;
end

if nargout > 3
    varargout{2} = reslvec;
end

if nargout > 4
    varargout{3} = rescvec;
end
end