function R = randald(N,alpha,lambda,kappa)
%RANDALD Generate pseudorandom numbers distributed asymmetric Laplace.
%   R = RANDN(N) returns an N-by-1 matrix containing pseudorandom values drawn
%   from the asymmetric Laplace.  

if kappa>=1 | kappa<=0
    error('The given kappa is in (0,1)')
end

e = random('exp',1,N,1);
z = randn(N,1);

R = alpha + ( e*(1-2*kappa)/(kappa*(1-kappa)) + sqrt(2*e/(kappa*(1-kappa))).*z ) / lambda;


