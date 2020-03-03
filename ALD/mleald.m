function [alpha,lambda,kappa,ll] = mleald(x,opt)
%% Initialization
if nargin == 1
    maxiter = 50;
    display = 1;
    tol = 1e-05;
else
    if (~isfield(opt,'maxiter'))
        maxiter = 100;
    else
        maxiter = opt.maxiter;
    end
    
    if (~isfield(opt,'display'))
        display = 1;
    else
        display = opt.display;
    end
    
    if (~isfield(opt,'tol'))
        tol = 1e-05;
    else
        tol = opt.tol;
    end
end
%% Iteration starts
x = x(:);
n = length(x);
kappa = 0.5; lambda = 1; alpha = mean(x);
ll = loglikald(x,alpha,lambda,kappa);
converged = false; t = 1;
while ~converged
    t = t+1;
    % lambda
    tempX = x-alpha;
    inx = (tempX>=0);
    rho = (~inx)*(1-kappa) + inx*kappa;
    lambda = n / sum(abs(tempX).*rho);
    if ~isfinite(lambda)
        lambda = n/1e-6;
    end
    % kappa
    eta = sum(tempX)*lambda;
    sDelta = sqrt(4*n^2+eta^2);
    kappa = (2*n + eta - sDelta) / (2*eta);
    % alpha
    alpha = weight_median(x,rho);
    % log-likelihood
    ll(t) = loglikald(x,alpha,lambda,kappa);
    % whether converge
    if t >= maxiter
        converged = true;
        if display; disp(['Not converged within ',num2str(t),' steps.']);end
    else
        converged = abs(ll(t)-ll(t-1))<tol;
    end
end
ll(1) = [];