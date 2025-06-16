function [x, varargout] = mgs_gmres(A, b, varargin)
% An implementation for `GMRES`. gmres solve the linear system
%       Ax = b
%
% or least square problems
%       min || Ax - b ||_2
%
% gmres finds the solution x which is belong to K(A,b).
%
%
% Syntaxes
% --------------
% x = gmres(A, b)
% x = gmres(A, b, tol)
% x = gmres(A, b, tol, maxit)
% x = gmres(A, b, tol, maxit, restart)
% x = gmres(A, b, tol, maxit, restart, x0)
%
%
% [x, exitflag] = gmres(__)
% [x, exitflag, resvec] = gmres(__)
% [x, exitflag, resvec, Aresvec] = gmres(__)
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
% restart           number of restart
%
% Name, Value       Name-Value pairs determine other options. Name should be in
%                   {'x0'}. 'x0' means initial gauss; Corresponding default
%                   values are {[]}.
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

if nargin == 0
    help gmres; return;
end

% default vaules of input parameters
defaultTol = 1e-6;
defaultMaxit = 100;
defaultRestart = [];
defaultX0 = [];

% check input parameters
checkMatrix = @(x) validateattributes(x, {'numeric', 'function'}, {'nonempty'});
checkVector = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'column'});
checkPostive = @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'});

p = inputParser;
p.CaseSensitive = true;

addRequired(p, 'A', checkMatrix);
addRequired(p, 'b', checkVector);
addOptional(p, 'tol', defaultTol, checkPostive);
addOptional(p, 'maxit', defaultMaxit, checkPostive);
addOptional(p, 'restart', defaultRestart, checkPostive);
addParameter(p, 'x0', defaultX0, checkVector);

parse(p, A, b, varargin{:});
Parameters = p.Results;

tol = Parameters.tol;
maxit = Parameters.maxit;
restart = Parameters.restart;
x0 = Parameters.x0;

if ~isempty(restart)
    maxit = ceil(maxit / restart);
    resvec = zeros(maxit * restart, 1);
    Aresvec = resvec;
    Atrans_resvec = resvec;
    last_idx = 0;
    for k = 1:maxit
        [x0, exitflag, resvec1, Aresvec1, Atrans_resvec1] = inner_gmres(A, b, tol, restart, x0);

        len = length(resvec1);
        resvec(last_idx+1:len+last_idx) = resvec1;
        Aresvec(last_idx+1:len+last_idx) = Aresvec1;
        Aresvec(last_idx+1:len+last_idx) = Atrans_resvec1;
        last_idx = last_idx + len;

        if exitflag == 0
            x = x0;
            resvec = resvec(1:last_idx);
            Aresvec = Aresvec(1:last_idx);
            Atrans_resvec = Atrans_resvec(1:last_idx);
            break;
        end
    end
else
    [x, exitflag, resvec, Aresvec, Atrans_resvec] = inner_gmres(A, b, tol, maxit, x0);
end

if nargout > 1
    varargout{1} = exitflag;
end

if nargout > 2
    varargout{2} = resvec;
end

if nargout > 3
    varargout{3} = Aresvec;
end

if nargout > 4
    varargout{4} = Atrans_resvec;
end
end

function [x, exitflag, resvec, Aresvec, Atrans_resvec] = inner_gmres(A, b, tol, maxit, x0)
m = size(A, 1);

if ~isempty(x0)
    r0 = b - A * x0;
else
    x0 = zeros(m, 1);
    r0 = b;
end

beta = norm(r0); v1 = r0 / beta;

V = zeros(m, maxit + 1); V(:, 1) = v1;
H = zeros(maxit + 1, maxit);
gc = zeros(maxit, 1); gs = gc;
z = gc; xib = beta;

exitflag = 1; 
resvec = gc;
Aresvec = resvec;
Atrans_resvec = resvec;
for j = 1:maxit
    v = A * v1;

    for i = 1:j
        vi = V(:, i); h = vi' * v; 
        v = v - h * vi;

        H(i, j) = h;
    end

    h = norm(v); H(j + 1, j) = h;

    % apply previous Givens rotations
    for k = 1:j-1
        c = gc(k); s = gs(k);
        H(k:k+1, j) = [c, s; -s, c] * H(k:k+1, j);
    end

    % apply current Givens rotation 
    [G, y] = planerot(H(j:j+1, j));
    H(j:j+1, j) = y; c = G(1, 1); s = G(1, 2);
    gc(j) = c; gs(j) = s;

    z(j) = xib * c; xib = -s * xib;

    y = H(1:j, 1:j) \ z(1:j); dx = V(:, 1:j) * y;
    r = r0 - A * dx;

    resvec(j) = norm(r);
    Aresvec(j) = norm(A*r);
    Atrans_resvec(j) = norm(A'*r);

    % resvec(j) = abs(xib);
    if resvec(j) < tol
        exitflag = 0;
        break;
    end

    % update for next iteration
    v1 = v / h; V(:, j+1) = v1;
end
y = H(1:j, 1:j) \ z(1:j);
x = x0 +  V(:, 1:j) * y;

resvec = resvec(1:j);
Aresvec = Aresvec(1:j);
Atrans_resvec = Atrans_resvec(1:j);
end