function [x, varargout] = lsmr(A, b, varargin)
% An implementation for `LSMR` in the paper:
%       "lsmr: an iterative algorithm for sparse least-squares problems"
% 
% LSMR solve the linear system
%       Ax = b 
% 
% or least square problems 
%       min || Ax - b ||_2 
% 
% LSMR is analytically equivalent to the MINRES applied to the 
% normal equation 
%       A'Ax = A'b.       
%
% 
% Syntaxes
% --------------
% x = lsmr(A, b) 
% x = lsmr(A, b, tol) 
% x = lsmr(A, b, tol, maxit) 
% x = lsmr(A, b, tol, maxit, x0) 
% 
% 
% [x, exitflag] = lsmr(__)
% [x, exitflag, resvec] = lsmr(__) 
% [x, exitflag, resvec, AResvec] = lsmr(__) 
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

beta = norm(r); u1 = r / beta;
v = A' * u1; alpha = norm(v);
v1 = v / alpha;

alphabar = alpha; xibar = alpha * beta;
rho0 = 1; rho0bar = rho0;
betadd = beta;beta0d = 0; 
tau = 0; theta0tilde = 0; 
xi0 = 0;

h = v1; hbar = zeros(n, 1);

exitflag = 1;
resvec = zeros(maxit, 1); AResvec = resvec;
for j = 1:maxit
    u = A * v1 - alpha * u1;
    beta = norm(u);
    u2 = u / beta;

    v = A' * u2 - beta * v1;
    alpha = norm(v);

    [G, y] = planerot([alphabar; beta]);
    rho = y(1); c = G(1, 1); s = G(1, 2);

    theta = s * alpha; alphabar = c * alpha;

    if j == 1
        thetabar = 0; rhohat = rho;
    else
        thetabar = sbar * rho;
        rhohat = cbar * rho;
    end

    [G, y] = planerot([rhohat; theta]);
    rho1bar = y(1); cbar = G(1, 1); sbar = G(1, 2);

    xi = cbar * xibar; xibar = -sbar * xibar;

    hbar = h - (thetabar * rho) / (rho0 * rho0bar) * hbar;
    x = x + xi / (rho * rho1bar) * hbar;
    
    % compute norm of residual
    betahat = c * betadd;
    betadd = -s * betadd;

    if j == 1
        rhod = rho1bar; betad = betahat;
        taud = xi / rhod; thetatilde = 0;
    else
        [G, y] = planerot([rhod; thetabar]);
        rho0tilde = y(1); ctilde = G(1, 2); 
        stilde = G(1, 2);

        thetatilde = stilde * rho1bar;
        rhod = ctilde * rho1bar;

        betad = -stilde * beta0d + ctilde * betahat;

        tau = (xi0 - theta0tilde * tau) / rho0tilde;
        taud = (xi - thetatilde * tau) / rhod;
    end
    
    % since loss of orthogonality, compute norm of r and Ar exactly
    r = b - A * x;
    AR = A * r;

    resvec(j) = norm(r);
    AResvec(j) = norm(AR);

    % ATransResvec(j) = abs(xibar);
    % resvec(j) = hypot(betad - taud, betadd);

    if resvec(j) < tol || AResvec(j) < tol
        exitflag = 0;
        break;
    end


    % update for next iteration
    u1 = u2; v1 = v / alpha;

    h = v1 - (theta / rho) * h;

    rho0 = rho; rho0bar = rho1bar;
    beta0d = betad;

    theta0tilde = thetatilde;
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