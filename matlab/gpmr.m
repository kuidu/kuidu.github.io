function [x, y, varargout] = gpmr(A, B, b, c, lambda, mu, varargin)
% GPMR will find approximate solution of following equation:
% 
%       [\lambda I    A] [x]   [b]
%       [              ] [ ] = [ ]
%       [B        \mu I] [y]   [c]
% 
%
% Syntaxes
% --------------
% [x, y] = gpmr(A, B, b, c) 
% [x, y] = gpmr(A, B, b, c, lambda, mu) 
% [x, y] = gpmr(A, B, b, c, lambda, mu, tol) 
% [x, y] = gpmr(A, B, b, c, lambda, mu, tol, maxit) 
% [x, y] = gpmr(A, B, b, c, lambda, mu, tol, maxit, restart)
% [x, y] = gpmr(A, B, b, c, lambda, mu, tol, maxit, restart, x0, y0)
% 
% [x, y, exitflag] = gpmr(__)
% [x, y, exitflag, resvec] = gpmr(__)
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
%                   {'tol', 'maxit', 'restart','x0', 'y0}. 'tol' means
%                   tolerance for convergence, 'maxit' means maximum
%                   iterations, 'restart' determine whether use restart
%                   tecnique, 'x0', 'y0' means initial gauss; Corresponding 
%                   default values are {1e-6, 100, 0, 0, 0}.
%               
%
% Returns
% --------------
% x, y              approximate solution for the above system
% 
% exitflag          covergence flag, `1` means failed, `0` means successed
% 
% resvec            vectors formed by norm of residuals 

if nargin == 0
    help gpmr; return;
end

% default vaules of input parameters
defaultTol = 1e-6;
defaultMaxit = 100; 
defaultRestart = 0; 
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
addOptional(p, 'restart', defaultRestart, checkPostive);
addOptional(p, 'x0', defaultX0, checkVector);
addOptional(p, 'y0', defaultY0, checkVector);

parse(p, A, B, b, c, lambda, mu, varargin{:});
Parameters = p.Results;

tol = Parameters.tol;
maxit = Parameters.maxit;
restart = Parameters.restart;
x0 = Parameters.x0;
y0 = Parameters.y0;

[m, n] = size(A);

if isempty(x0) && isempty(y0)
    x0 = zeros(m, 1);
    y0 = zeros(n, 1);
end

resvec = [];

if restart > 0
    for i = 1:maxit
        [x, y, exitflag, resvec1] = inner_gpmr(A, B, b, c, lambda, mu, tol, restart, x0, y0);
        resvec = [resvec; resvec1];

        if exitflag == 0
            break;
        end

        x0 = x; y0 = y;
    end
else
   [x, y, exitflag, resvec] = inner_gpmr(A, B, b, c, lambda, mu, tol, maxit, x0, y0);
end

if nargout > 2
    varargout{1} = exitflag;
end

if nargout > 3
    varargout{2} = resvec;
end

end

function [x, y, exitflag, resvec] = inner_gpmr(A, B, b, c, lambda, mu, tol, maxit, x0, y0)
[m, n] = size(A);

b = b - (lambda * x0 + A * y0);
c = c - (B * x0 + mu * y0);

beta = norm(b); 
gamma = norm(c);

v1 = b / beta;
u1 = c / gamma;

exitflag = 1;

V = zeros(m, maxit);
U = zeros(n, maxit);
V(:, 1) = v1;
U(:, 1) = u1;

S = zeros(2*maxit+2, 2*maxit);
gc = zeros(4*maxit, 1);
gs = zeros(4*maxit, 1);
t = zeros(2*maxit+2, 1);
t(1) = beta; t(2) = gamma;

resvec = zeros(maxit, 1);

for k = 1:maxit
    q = A * u1;
    p = B * v1;

    for i = 1:k
        v = V(:, i);
        u = U(:, i);

        h = v' * q;
        f = u' * p;

        q = q - h * v;
        p = p - f * u;

        S(2*i-1, 2*k) = h;
        S(2*i, 2*k-1) = f;

        if i == k
            S(2*i-1, 2*i-1) = lambda;
            S(2*i, 2*i) = mu;
        end
    end

    h = norm(q); f = norm(p);
    % v2 = q / h; u2 = p / f;
    
    S(2*k+1, 2*k) = h;
    S(2*k+2, 2*k-1) = f;

    for i = 1:k-1
        c1 = gc(4*i-3); c2 = gc(4*i-2); 
        c3 = gc(4*i-1); c4 = gc(4*i);

        s1 = gs(4*i-3); s2 = gs(4*i-2); 
        s3 = gs(4*i-1); s4 = gs(4*i);

        S([2*i-1, 2*i+2], 2*k-1:2*k) = [c1, s1; -s1, c1] * S([2*i-1, 2*i+2], 2*k-1:2*k);
        S([2*i-1, 2*i], 2*k-1:2*k) = [c2, s2; -s2, c2] * S([2*i-1, 2*i], 2*k-1:2*k);
        S([2*i, 2*i+2], 2*k-1:2*k) = [c3, s3; -s3, c3] * S([2*i, 2*i+2], 2*k-1:2*k);
        S([2*i, 2*i+1], 2*k-1:2*k) = [c4, s4; -s4, c4] * S([2*i, 2*i+1], 2*k-1:2*k);
    end

    [G, y] = planerot(S([2*k-1, 2*k+2], 2*k-1));
    S([2*k-1, 2*k+2], 2*k-1) = y;
    S([2*k-1, 2*k+2], 2*k) = G * S([2*k-1, 2*k+2], 2*k);

    t([2*k-1, 2*k+2]) = G * t([2*k-1, 2*k+2]);
    c1 = G(1, 1); s1 = G(1, 2);

    [G, y] = planerot(S([2*k-1, 2*k], 2*k-1));
    S([2*k-1, 2*k], 2*k-1) = y;
    S([2*k-1, 2*k], 2*k) = G * S([2*k-1, 2*k], 2*k);

    t([2*k-1, 2*k]) = G * t([2*k-1, 2*k]);
    c2 = G(1, 1); s2 = G(1, 2);

    [G, y] = planerot(S([2*k, 2*k+2], 2*k));
    S([2*k, 2*k+2], 2*k) = y;

    t([2*k, 2*k+2]) = G * t([2*k, 2*k+2]);
    c3 = G(1, 1); s3 = G(1, 2);

    [G, y] = planerot(S([2*k, 2*k+1], 2*k));
    S([2*k, 2*k+1], 2*k) = y;

    t([2*k, 2*k+1]) = G * t([2*k, 2*k+1]);
    c4 = G(1, 1); s4 = G(1, 2);

    gc(4*k-3:4*k) = [c1; c2; c3; c4];
    gs(4*k-3:4*k) = [s1; s2; s3; s4];

    % resvec(k+1) = sqrt(t(2*k+1)^2 + t(2*k+2)^2);
    % for fair comparison purpose, compute every residual explicitly    
    zeta = S(1:2*k, 1:2*k) \ t(1:2*k);
    dx = V(:, 1:k) * zeta(1:2:end);
    dy = U(:, 1:k) * zeta(2:2:end);
    rx = b - (lambda * dx + A * dy);
    ry = c - (B * dx + mu * dy);
    resvec(k) = hypot(norm(rx), norm(ry)); 
    
    if resvec(k) < tol
        exitflag = 0;
        break;
    end

    v1 = q / h; u1 = p / f;
    V(:, k+1) = v1;
    U(:, k+1) = u1;
end

resvec = resvec(1:k);
zeta = S(1:2*k, 1:2*k) \ t(1:2*k);
x = V(:, 1:k) * zeta(1:2:end) + x0;
y = U(:, 1:k) * zeta(2:2:end) + y0;
end