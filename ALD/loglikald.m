function ll = loglikald(x,alpha,lambda,kappa)
% Compute log-likelihood for ALD while avoiding numerical underflow.
x = x(:);
logpdf = logpdfald(x,alpha,lambda,kappa);
M = max(logpdf);
logpdf = logpdf - M;
ll = M + log(sum(exp(logpdf)));
if isfinite(ll)
    ll = M;
end
end