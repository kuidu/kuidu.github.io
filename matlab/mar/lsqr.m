function [x, varargout] = lsqr(A, b, varargin)
% An implementation for `LSQR` in the paper:
%       "lsqr: an iterative algorithm for sparse least-squares problems"
% 
% LSQR solve the linear system
%       Ax = b 
% 
% or least square problems 
%       min || Ax - b ||_2 
% 
% LSQR is analytically equivalent to the CG applied to the 
% normal equation 
%       A'Ax = A'b.       
%
% 
% Syntaxes
% --------------
% x = lsqr(A, b) 
% x = lsqr(A, b, tol) 
% x = lsqr(A, b, tol, maxit) 
% x = lsqr(A, b, tol, maxit, x0) 
% 
% 
% [x, exitflag] = lsqr(__)
% [x, exitflag, resvec] = lsqr(__) 
% [x, exitflag, resvec, AResvec] = lsqr(__) 
%
%
% Parameters
% ------------ 
% A                 m by n metrix 
% 
% b                 m by 1 right hand vector      
% 
% tol               stop tolerance 
% 
% maxit             number of maximun iterations
% 
% Name, Value       Name-Value pairs determine other options. Name should be in 
%                   {'tol', 'maxit','x0'}. 'tol' means tolerance for convergence, 
%                   'maxit' means maximum iterations, 'x0' means initial gauss; 
%                   Corresponding default values are {1e-6, 100, []}.
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
% AResvec      vectors formed by norm of Ar, r means residuals

if nargin == 0
    help lsmr; return;
end

% default vaules of input parameters
defaultTol = 1e-6;
defaultMaxit = 100; 
defaultX0 = [];

% check input parameters
checkMatrix = @(x) validateattributes(x, {'numeric'}, {'nonempty'});
checkVector = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'column'});
checkPostive = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});

p = inputParser;
p.CaseSensitive = true;

addRequired(p, 'A', checkMatrix);
addRequired(p, 'b', checkVector);
addOptional(p, 'tol', defaultTol, checkPostive);
addOptional(p, 'maxit', defaultMaxit, checkPostive);
addOptional(p, 'x0', defaultX0, checkVector);

parse(p, A, b, varargin{:});
Parameters = p.Results;

tol = Parameters.tol;
maxit = Parameters.maxit;
x0 = Parameters.x0;

n = size(A, 2);

if isempty(x0)
    x = zeros(n, 1);
    r = b;
else
    x = x0;
    r = b - A * x0;
end

beta = norm(r); u = r / beta;
v = A' * u; alpha = norm(v);
v = v / alpha;
phi_bar = beta;
rho_bar = alpha; 
exitflag = 1;
resvec = zeros(maxit, 1); AResvec = resvec;
w = v;
for j = 1:maxit
    u = A * v - alpha * u;
    beta = norm(u);
    u = u / beta;

    v = A' * u - beta * v;
    alpha = norm(v);
    v = v / alpha;
    [G, y] = planerot([rho_bar; beta]);
    rho = y(1); c = G(1, 1); s = G(1, 2);

    theta = s * alpha; rho_bar = -c * alpha;
    phi = c * phi_bar;
    phi_bar = s * phi_bar;
    % update x and w
    x = x +( phi / rho) * w;
    w = v - (theta / rho) *w;
    % compute norm of residual and ATrans residual
     % resvec(j) = phi_bar;
     % ATransResvec(j) =phi_bar * alpha * abs(c);
    
    % since loss of orthogonality, compute norm of r and A'r exactly
    r = b - A * x;
     AR = A * r;

    resvec(j) = norm(r);
    AResvec(j) = norm(AR);
    if resvec(j) < tol || AResvec(j) < tol
        exitflag = 0;
        break;
    end   
end


if nargout > 1
    varargout{1} = exitflag;
end

if nargout > 2
    varargout{2} = resvec(1:j);
end

if nargout > 3
    varargout{3} = AResvec(1:j);
end
end