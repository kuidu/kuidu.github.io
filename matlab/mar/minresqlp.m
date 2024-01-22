function [x, varargout] = minresqlp(A, b, varargin)
% An implementation for `MINRES-QLP` in the paper:
%       "minres-qlp: a krylov subspace method for indefinte OR singular symmetric systems"
% 
% MINRES-QLP solves the linear system
%       Ax = b 
% 
% or the least squares problem 
%       min || Ax - b ||_2 
% 
% MINRES-QLP will find the solution of the least squares problem of minimum
% norm
%
% 
% Syntaxes
% --------------
% x = minresqlp(A, b) 
% x = minresqlp(A, b, tol) 
% x = minresqlp(A, b, tol, maxit) 
% x = minresqlp(A, b, tol, maxit, x0) 
% x = minresqlp(A, b, tol, maxit, x0, criteria) 
% 
% 
% [x, exitflag] = minresqlp(__)
% [x, exitflag, resvec] = minresqlp(__) 
% [x, exitflag, resvec, Aresvec] = minresqlp(__) 
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
% Name, Value       Name-Value pairs determine other options. Name should be in 
%                   {'x0'}. 'x0' means initial gauss; Corresponding default 
%                   values are [].
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
% Aresvec      vectors formed by norm of A'r, r means residuals
%
% Kui Du, Jia-Jun Fan, and Fang Wang 2024.01.21
%

if nargin == 0
    help minresqlp; return;
end

% default vaules of input parameters
defaultTol = 1e-6;
defaultMaxit = 100; 
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
addParameter(p, 'x0', defaultX0, checkVector);
addOptional(p, 'criteria', defaultCriteria, checkPostive);

parse(p, A, b, varargin{:});
Parameters = p.Results;

tol = Parameters.tol;
maxit = Parameters.maxit;
x0 = Parameters.x0;
criteria = Parameters.criteria;

n = size(A, 1);

if isempty(x0)
    x = zeros(n, 1);
    r = b;
else
    x = x0;
    r = b - A * x0;
end

v0 = zeros(n, 1);
beta1 = norm(r);
v1 = r / beta1;

xt = x; 
wd_1 = v0; wdd0 = wd_1;

xi_1 = 0; xi0 = 0; xibar = beta1;
lamb_1 = 0; etab_1 = 0; etab0 = 0;
mu_3 = 0; mu_2 = 0;

deltad0 = 0;

exitflag = 1;
resvec = zeros(maxit, 1);
Aresvec = resvec;
for j = 1:maxit
    q = A * v1 - beta1 * v0;
    alpha = v1' * q;

    v = q - alpha * v1;
    beta2 = norm(v);

    if j == 1
        alphat = alpha; betat2 = beta2;
        lam = 0; eta2 = 0;
    else
        lam = c * betat + s * alpha;
        alphat = -s * betat + c * alpha;

        eta2 = s * beta2; betat2 = c * beta2;
    end

    [G, y] = planerot([alphat; beta2]);
    delta = y(1); c = G(1, 1); s = G(1, 2);

    xi = c * xibar; xibar = -s * xibar;

    if j <= 2
        lamb0 = 0; lamt = lam;
        etab = 0; deltat = delta;

        w_1 = wd_1; wt = v1;
        mu_1 = 0;
    else
        [G, y] = planerot([deltad_1; eta]);
        deltab_1 = y(1); ct1 = G(1, 1); st1 = G(1, 2);

        lamb0 = ct1 * lamd0 + st1 * lam;
        lamt = -st1 * lamd0 + ct1 * lam;

        etab = st1 * delta; deltat = ct1 * delta;

        w_1 = ct1 * wd_1 + st1 * v1;
        wt = -st1 * wd_1 + ct1 * v1;

        mu_1 = (xi_1 - lamb_1 * mu_2 - etab_1 * mu_3) / deltab_1;
    end

    if j == 1
        lamd = 0; deltadd = deltat;

        wd0 = wdd0; wdd = wt;
        mud0 = 0;
    else
        [G, y] = planerot([deltadd0; lamt]);
        deltad0 = y(1); ct = G(1, 1); st = G(1, 2);

        lamd = st * deltat; deltadd = ct * deltat;

        wd0 = ct * wdd0 + st * wt;
        wdd = -st * wdd0 + ct * wt;

        mud0 = (xi0 - lamb0 * mu_1 - etab0 * mu_2) / deltad0;
    end

    % mu_1 = (xi_1 - lamb_1 * mu_2 - etab_1 * mu_3) / deltab_1;
    % mud0 = (xi0 - lamb0 * mu_1 - etab0 * mu_2) / deltad0;

    if abs(deltadd) < 1e-12
        mudd = 0;
    else
        mudd = (xi - lamd * mud0 - etab * mu_1) / deltadd;
    end

    xt = xt + mu_1 * w_1;
    x = xt + mud0 * wd0 + mudd * wdd;
    
    % since loss of orthogonality, compute residual exactly
    r = b - A * x;
    Ar = A * r;

    resvec(j) = norm(r);
    Aresvec(j) = norm(Ar);

    % resvec(j) = abs(xibar);
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
    beta1 = beta2;
    v0 = v1; v1 = v / beta1;

    betat = betat2;
    deltad_1 = deltad0; eta = eta2;
    lamd0 = lamd; deltadd0 = deltadd;

    wd_1 = wd0; wdd0 = wdd;

    xi_1 = xi0; xi0 = xi;
    lamb_1 = lamb0; etab_1 = etab0;
    etab0 = etab;

    mu_3 = mu_2; mu_2 = mu_1;
end
 
if nargout > 1
    varargout{1} = exitflag;
end

if nargout > 2
    varargout{2} = resvec(1:j);
end

if nargout > 3
    varargout{3} = Aresvec(1:j);
end
end