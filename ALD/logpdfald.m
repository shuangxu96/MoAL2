function logy = logpdfald(x,alpha,lambda,kappa)
x = x(:);
xa = x-alpha;
xabs = abs(xa);
inx = (xa >= 0);
al_kernal = lambda*xabs.*( kappa*inx + (1-kappa)*~inx );
logy = log(lambda*kappa*(1-kappa)) - al_kernal;