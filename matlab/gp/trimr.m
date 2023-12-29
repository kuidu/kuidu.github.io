function [x, y, varargout] = trimr(A, b, c, lambda, mu, varargin)
% A implementation for `TRIMR` in the paper:
%   "trimr and trimr: two iterative methods for symmetric qusi-definite systems"
%   
% This method solve the following system:
%       [\lambda I    A] [x]   [b]
%       [              ] [ ] = [ ]
%       [A'       \mu I] [y]   [c]
%
% Syntaxes
% --------------
% [x, y] = trimr(A, b, c) 
% [x, y] = trimr(A, b, c, lambda, mu) 
% [x, y] = trimr(A, b, c, lambda, mu, tol) 
% [x, y] = trimr(A, b, c, lambda, mu, tol, maxit) 
% [x, y] = trimr(A, b, c, lambda, mu, tol, maxit, x0, y0) 
% 
% 
% [x, y, exitflag] = trimr(__)
% [x, y, exitflag, resvec] = trimr(__)
%
%
% Parameters
% ------------ 
% A                 m by n metrix 
% 
% b, c              m by 1 and n by 1 vector, respectively 
% 
% lambda, mu        coefficients of matrix        
% 
% tol               stop tolerance 
% 
% maxit             maximun iterations
% 
% Name, Value       Name-Value pairs determine other options. Name should be in 
%                   {'tol', 'maxit','x0', 'y0}. 'tol' means
%                   tolerance for convergence, 'maxit' means maximum
%                   iterations, 'x0', 'y0' means initial gauss; Corresponding 
%                   default values are {1e-6, 100, [], []}.
% 
%
% Returns
%-------------
% x             approximate solution for the above system
% 
% exitflag      covergency flag, `1` means failed, `0` means successed
% 
% resvec        vectors formed by norm of residuals 

if nargin == 0
    help trimr; return;
end

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
addRequired(p, 'b', checkVector);
addRequired(p, 'c', checkVector);
addRequired(p, 'lambda', checkScalar);
addRequired(p, 'mu', checkScalar);

addOptional(p, 'tol', defaultTol, checkPostive);
addOptional(p, 'maxit', defaultMaxit, checkPostive);
addOptional(p, 'x0', defaultX0, checkVector);
addOptional(p, 'y0', defaultY0, checkVector);

parse(p, A, b, c, lambda, mu, varargin{:});
Parameters = p.Results;

tol = Parameters.tol;
maxit = Parameters.maxit;
x0 = Parameters.x0;
y0 = Parameters.y0;

[m, n] = size(A);

if isempty(x0) && isempty(y0)
    x = zeros(m, 1); y = zeros(n, 1);
    rx = b; ry = c;
else
    rx = b - (lambda * x0 + A * y0);
    ry = c - (A' * x0 + mu * y0);

    x = x0; y = y0;
end

gx_3 = zeros(m, 1); gx_2 = gx_3;
gx_1 = gx_3; gx0 = gx_3;

gy_3 = zeros(n, 1); gy_2 = gy_3;
gy_1 = gy_3; gy0 = gy_3;


v0 = gx_3; u0 = gy_3;

v1 = rx; u1 = ry;
beta1 = norm(v1); gamma1 = norm(u1);
v1 = v1 / beta1; u1 = u1 / gamma1;

% norm_r = sqrt(beta1^2 + gamma1^2);

mu_3 = 0; mu_2 = 0;
lambda_2 = 0;

pi1bar = beta1;
pi2bar = gamma1;

resvec = zeros(maxit, 1);
exitflag = 1;
for k = 1:maxit
    q = A * u1 - gamma1 * v0;
    p = A' * v1 - beta1 * u0;

    alpha = v1' * q;

    v = q - alpha * v1;
    u = p - alpha * u1;
    
    beta2 = norm(v); gamma2 = norm(u);

    if k == 1
        thetabar = alpha; delta1bar = lambda;
        delta2bar = mu; sigma1bar = alpha;
        eta1bar = 0; lambda1bar = gamma2;
        sigma2bar = beta2;

        eta_1 = 0; eta0 = 0;
        lambda_1 = 0; lambda0 = 0;
        mu0 = 0; mu_1 = 0;
        sigma0 = 0;
    else
        % apply previous Givens rotations
        % first Givens rotation
        sigma0tilde = c1 * sigma0bar + s1 * alpha;
        thetatilde = -s1 * sigma0bar + c1 * alpha;
        eta0tilde = mu * s1;
        delta2tilde = mu * c1;
        lambda0tilde = beta2 * s1;
        sigma2tilde = beta2 * c1;

        % second Givens rotation
        eta_1 = c2 * eta_1bar + s2 * sigma0tilde;
        sigma0hat = -s2 * eta_1bar + c2 * sigma0tilde;
        lambda_1 = c2 * lambda_1bar + s2 * eta0tilde;
        eta0hat = -s2 * lambda_1bar + c2 * eta0tilde;
        mu_1 = s2 * lambda0tilde;
        lambda0hat = c2 * lambda0tilde;

        % third Givens rotation
        sigma0o = c3 * sigma0hat + s3 * thetatilde;
        thetabar = -s3 * sigma0hat + c3 * thetatilde;
        eta0o = c3 * eta0hat + s3 * delta2tilde;
        delta2bar = -s3 * eta0hat + c3 * delta2tilde;
        lambda0o = c3 * lambda0hat + s3 * sigma2tilde;
        sigma2bar = -s3 * lambda0hat + c3 * sigma2tilde;

        % forth Givens rotation
        sigma0 = c4 * sigma0o + s4 * lambda;
        delta1bar = -s4 * sigma0o + c4 * lambda;
        eta0 = c4 * eta0o + s4 * alpha;
        sigma1bar = -s4 * eta0o + c4 * alpha;
        lambda0 = c4 * lambda0o;
        eta1bar = -s4 * lambda0o;
        mu0 = s4 * gamma2;
        lambda1bar = c4 * gamma2;
    end

    % first Givens rotation
    [G, r] = planerot([thetabar; gamma2]);
    theta = r(1); c1 = G(1, 1); s1 = G(1, 2);
    delta2tilde = c1 * delta2bar;
    rho1 = -s1 * delta2bar;

    pi2tilde = c1 * pi2bar;
    pi4tilde = -s1 * pi2bar;

    % second Givens rotation
    [G, r] = planerot([delta1bar; theta]);
    delta1 = r(1); c2 = G(1, 1); s2 = G(1, 2);
    sigma1 = c2 * sigma1bar + s2 * delta2tilde;
    delta2hat = -s2 * sigma1bar + c2 * delta2tilde;

    pi1 = c2 * pi1bar + s2 * pi2tilde;
    pi2hat = -s2 * pi1bar + c2 * pi2tilde;

    % third Givens rotation
    [G, r] = planerot([delta2hat; rho1]);
    delta2o = r(1); c3 = G(1, 1); s3 = G(1, 2);

    pi2o = c3 * pi2hat + s3 * pi4tilde;
    pi4bar = -s3 * pi2hat + c3 * pi4tilde;

    % forth Givens rotation
    [G, r] = planerot([delta2o; beta2]);
    delta2 = r(1); c4 = G(1, 1); s4 = G(1, 2);

    pi2 = c4 * pi2o;
    pi3bar = -s4 * pi2o;

    % update approximate solution
    gx1 = (v1 - mu_3 * gx_3 - lambda_2 * gx_2 - eta_1 * gx_1 - sigma0 * gx0) / delta1;
    gx2 = -(mu_2 * gx_2 + lambda_1 * gx_1 + eta0 * gx0 + sigma1 * gx1) / delta2;

    gy1 = -(mu_3 * gy_3 + lambda_2 * gy_2 + eta_1 * gy_1 + sigma0 * gy0) / delta1;
    gy2 = (u1 - mu_2 * gy_2 - lambda_1 * gy_1 - eta0 * gy0 - sigma1 * gy1) / delta2;

    dx = pi1 * gx1 + pi2 * gx2;
    dy = pi1 * gy1 + pi2 * gy2;

    x = x + dx; y = y + dy;

    rx = b - (lambda * x + A * y);
    ry = c - (A' * x + mu * y);

    resvec(k) = hypot(norm(rx), norm(ry));

    if resvec(k) < tol
        exitflag = 0;
        break;
    end

    % update for next iteration
    u0 = u1; u1 = u / gamma2;
    v0 = v1; v1 = v / beta2; 

    beta1 = beta2; gamma1 = gamma2;

    pi1bar = pi3bar; pi2bar = pi4bar;

    gx_3 = gx_1; gx_2 = gx0;
    gx_1 = gx1; gx0 = gx2;

    gy_3 = gy_1; gy_2 = gy0;
    gy_1 = gy1; gy0 = gy2;

    mu_3 = mu_1; mu_2 = mu0;
    lambda_2 = lambda0; % lambda_1 = lambda1;
    % eta_1 = eta1; eta0 = eta2;
    % sigma0 = sigma2;

    sigma0bar = sigma2bar;
    eta_1bar = eta1bar;
    lambda_1bar = lambda1bar;
end

resvec = resvec(1:k);

if nargout > 2
    varargout{1} = exitflag;
end

if nargout > 3
    varargout{2} = resvec;
end
end