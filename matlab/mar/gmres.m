function [x, varargout] = gmres(A, b, varargin)
% An implementation for `GMRES`. gmres solves the linear system
%       Ax = b
%
% or the least squares problem
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
% x = gmres(A, b, tol, maxit, restart, x0, criteria)
%
%
% [x, exitflag] = gmres(__)
% [x, exitflag, resvec] = gmres(__)
% [x, exitflag, resvec, Aresvec] = gmres(__)
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
% restart           number of restart
%
% criteria          stop criteria 
%
% Name, Value       Name-Value pairs determine other options. Name should be in
%                   {'x0'}. 'x0' means initial guess; Corresponding default
%                   value are {[]}.
%                   'criteria' means stop criteria, '1' means ||r||<tol, '2'means 
%                   ||Ar||<tol, '3' means ||r||<tol or ||Ar||<tol.
%                   Corresponding default values are {[], 2, 3}
%
%
% Returns
%-------------
% x                 approximate solution for the above system
%
% exitflag          covergence flag, `1` means failed, `0` means successed
%
% resvec            vectors formed by norm of residuals
%
% Aresvec           vectors formed by norm of Ar, r means residuals
%
% Kui Du, Jia-Jun Fan, and Fang Wang 2024.01.21
%

if nargin == 0
    help gmres; return;
end

% default vaules of input parameters
defaultTol = 1e-6;
defaultMaxit = 100;
defaultRestart = [];
defaultX0 = [];
defaultCriteria = 3;

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
addOptional(p, 'restart', defaultRestart, checkPostive);
addParameter(p, 'x0', defaultX0, checkVector);
addOptional(p, 'criteria', defaultCriteria, checkPostive);

parse(p, A, b, varargin{:});
Parameters = p.Results;

tol = Parameters.tol;
maxit = Parameters.maxit;
restart = Parameters.restart;
x0 = Parameters.x0;
criteria = Parameters.criteria;

if ~isempty(restart)
    maxit = ceil(maxit / restart);
    resvec = zeros(maxit * restart, 1); Aresvec = resvec;
    last_idx = 0;
    for k = 1:maxit
        [x0, exitflag, resvec1, Aresvec1] = inner_gmres(A, b, tol, restart, x0, criteria);

        len = length(resvec1);
        resvec(last_idx+1:len+last_idx) = resvec1;
        Aresvec(last_idx+1:len+last_idx) = Aresvec1;
        last_idx = last_idx + len;

        if exitflag == 0
            x = x0;
            resvec = resvec(1:last_idx);
            Aresvec = Aresvec(1:last_idx);
            break;
        end
    end
else
    [x, exitflag, resvec, Aresvec] = inner_gmres(A, b, tol, maxit, x0, criteria);
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
end

function [x, exitflag, resvec, Aresvec] = inner_gmres(A, b, tol, maxit, x0, criteria)
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

exitflag = 1; resvec = gc; Aresvec = gc;
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
    Ar = A * r;

    resvec(j) = norm(r); Aresvec(j) = norm(Ar);

    % resvec(j) = abs(xib);
    if criteria==1
        if resvec(j) < tol
            exitflag = 0;
            break;
        end
    elseif criteria==2
        if  Aresvec(j) < tol
            exitflag = 0;
            break;
        end
    elseif  criteria==3
        if resvec(j) < tol || Aresvec(j) < tol
            exitflag = 0;
            break;
        end
    end

    % update for next iteration
    v1 = v / h; V(:, j+1) = v1;
end
y = H(1:j, 1:j) \ z(1:j);
x = x0 +  V(:, 1:j) * y;
end