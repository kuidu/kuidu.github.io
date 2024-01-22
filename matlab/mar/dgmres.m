function [x, varargout] = dgmres(A, b, varargin)
% An implementation for `DGMRES`. dgmres solves the linear system
%       Ax = b 
% 
% or the least squares problem 
%       min || Ax - b ||_2 
% 
% dgmres finds the solution x which is belong to K(A, Ab). 
%
% 
% Syntaxes
% --------------
% x = dgmres(A, b) 
% x = dgmres(A, b, tol) 
% x = dgmres(A, b, tol, maxit) 
% x = dgmres(A, b, tol, maxit, restart)
% x = dgmres(A, b, tol, maxit, restart, x0) 
% x = dgmres(A, b, tol, maxit, restart, x0, criteria) 
% 
% 
% [x, exitflag] = dgmres(__)
% [x, exitflag, resvec] = dgmres(__) 
% [x, exitflag, resvec, Aresvec] = dgmres(__) 
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
%                   {'x0'}. 'x0' means initial gauss; Corresponding default 
%                   values are {[], 2}.
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
    help dgmres; return;
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
        [x0, exitflag, resvec1, Aresvec1] = inner_dgmres(A, b, tol, restart, x0, criteria);
        
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
    [x, exitflag, resvec, Aresvec] = inner_dgmres(A, b, tol, maxit, x0, criteria);
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

function [x, exitflag, resvec, Aresvec] = inner_dgmres(A, b, tol, maxit, x0, criteria)
m = size(A, 1);

if ~isempty(x0)
    r0 = b - A * x0;
else
    x0 = zeros(m, 1);
    r0 = b;
end

Ar0 = A * r0;
beta = norm(Ar0); v1 = Ar0 / beta;

v = A * v1; h11 = v1' * v;
v = v - h11 * v1; h21 = norm(v); 
v2 = v / h21;

gc = zeros(maxit + 2, 1); gs = gc;
gct = gc; gst = gc;

V = zeros(m, maxit + 2);
H = zeros(maxit + 2, maxit + 1);
R = zeros(maxit + 1);
z = gc;

[G, y] = planerot([h11; h21]);

gc(1) = G(1, 1); gs(1) = G(1, 2);

% beta1 = beta * h11; beta2 = beta * h21;
xib = gc(1) * beta; xibb = -gs(1) * beta;

H(1, 1) = y(1); R(1, 1) = y(1);
V(:, 1) = v1; V(:, 2) = v2; 
 
exitflag = 1; 
resvec = zeros(maxit, 1); Aresvec = zeros(maxit, 1);
for j = 1:maxit
    v = A * v2;
    for i = 1:j+1
        vi = V(:, i); 
        h = vi' * v; H(i, j+1) = h;

        v = v - h * vi;
    end

    h = norm(v); H(j+2, j+1) = h;
    
    % apply previous Givens rotations
    for k = 1:j
        c = gc(k); s = gs(k);
        H(k:k+1, j+1) = [c, s; -s, c] * H(k:k+1, j+1);
    end

    % apply current Givens rotation 
    [G, y] = planerot(H(j+1:j+2, j+1));
    H(j+1:j+2, j+1) = y; c = G(1, 1); s = G(1, 2);
    gc(j+1) = c; gs(j+1) = s;

    R(1:j+1, j+1) = H(1:j+1, j+1);

    xit = c * xibb; xibb = -s * xibb;
  
    % H(1:j+1, j:j+1) = H(1:j+1, j:j+1) * [c, -s; s, c];

    % apply second previous Givens rotations
    for k = 1:j-1
        c = gct(k); s = gst(k);
        H(k:k+1, j+1) = [c, s; -s, c] * H(k:k+1, j+1);
    end
    
    c = gc(j); s = gs(j);
    H(1:j+1, j:j+1) = H(1:j+1, j:j+1) * [c, -s; s, c];

    % apply second Givens rotation
    [G, y] = planerot(H(j:j+1, j));
    H(j:j+1, j) = y; c = G(1, 1); s = G(1, 2);
    gct(j) = c; gst(j) = s;
    
    H(j:j+1, j+1) = [c, s; -s, c] * H(j:j+1, j+1);
    z(j) = c * xib + s * xit; xib = -s * xib + c * xit;
    
    % compute iterate solution exactly, for obtaining some wanted data
    y = H(1:j, 1:j) \ z(1:j); y = R(1:j, 1:j) \ y;
    dx = V(:, 1:j) * y;
    r = r0 - A * dx;
    Ar = A * r;

    resvec(j) = norm(r);
    Aresvec(j) = norm(Ar);

    % Aresvec(j) = hypot(xib, xibb);
    % if Range(A) = Range(A'), use following stop criteria
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
    v2 = v / h; V(:, j+2) = v2;
end

resvec = resvec(1:j);
Aresvec = Aresvec(1:j);

y = H(1:j, 1:j) \ z(1:j); y = R(1:j, 1:j) \ y;
x = x0 + V(:, 1:j) * y;
end