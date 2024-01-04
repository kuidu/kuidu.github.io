function [x, y, varargout] = tricg(A, b, c, lambda, mu, varargin)
% A implementation for `TRICG` in the paper:
%   "tricg and trimr: two iterative methods for symmetric qusi-definite systems"
%   
% This method solve the following system:
%       [\lambda I    A] [x]   [b]
%       [              ] [ ] = [ ]
%       [A'       \mu I] [y]   [c]
%
% Syntaxes
% --------------
% [x, y] = tricg(A, b, c) 
% [x, y] = tricg(A, b, c, lambda, mu) 
% [x, y] = tricg(A, b, c, lambda, mu, tol) 
% [x, y] = tricg(A, b, c, lambda, mu, tol, maxit) 
% [x, y] = tricg(A, b, c, lambda, mu, tol, maxit, x0, y0) 
% 
% 
% [x, y, exitflag] = tricg(__)
% [x, y, exitflag, resvec] = tricg(__)
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
%
% Date
% --------------
% Kui Du, Jia-Jun Fan, and Fang Wang, 2024.1.4 


if nargin == 0
    help tricg; return;
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

% initialization
gx_1 = zeros(m, 1); gx0 = gx_1; 
gy_1 = zeros(n, 1); gy0 = gy_1;

v0 = gx0; u0 = gy0;

v1 = rx; u1 = ry;
beta1 = norm(v1); gamma1 = norm(u1);
v1 = v1 / beta1; u1 = u1 / gamma1;

% norm_r = sqrt(beta1^2 + gamma1^2);
d_1 = 0; d0 = 0; sig1 = 0; 
eta1 = 0; lam1 = 0; delta0 = 0;
pi_1 = 0; pi0 = 0;

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
        d1 = lambda;
        delta1 = alpha / d1;
        d2 = mu - delta1^2 * d1;
    else
        sig1 = beta1 / d0;
        eta1 = gamma1 / d_1;
        lam1 = -(eta1 * delta0 * d_1) / d0;
        d1 = lambda - sig1^2 * d0;
        delta1 = (alpha - lam1 * sig1 * d0) / d1;
        d2 = mu - eta1^2 * d_1 - lam1^2 * d0 - delta1^2 * d1;
    end

    if k == 1
        pi1 = beta1 / d1;
        pi2 = (gamma1 - delta1 * beta1) / d2;
    else
        pi1 = -(sig1 * pi0 * d0) / d1;
        pi2 = -(delta1 * pi1 * d1 + lam1 * pi0 * d0 + eta1 * pi_1 * d_1) / d2;
    end

    % update approximate solution
    gx1 = v1 - sig1 * gx0;
    gx2 = -delta1 * gx1 - lam1 * gx0 - eta1 * gx_1;

    gy1 = -sig1 * gy0;
    gy2 = u1 - delta1 * gy1 - lam1 * gy0 - eta1 * gy_1;

    x = x + pi1 * gx1 + pi2 * gx2;
    y = y + pi1 * gy1 + pi2 * gy2;

    % norm_r = sqrt( (gamma2 * (pi1 - delta1 * pi2))^2 + (beta2 * pi2)^2 );
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
    gamma1 = gamma2; beta1 = beta2;
    

    d_1 = d1; d0 = d2;
    pi_1 = pi1; pi0 = pi2;
    delta0 = delta1;

    gx_1 = gx1; gx0 = gx2;
    gy_1 = gy1; gy0 = gy2;
end


if nargout > 2
   varargout{1} = exitflag; 
end

if nargout > 3
    varargout{2} = resvec(1:k);
end

end