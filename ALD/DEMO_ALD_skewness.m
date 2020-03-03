alpha=0;
Kappa=0.01:0.01:0.99;
lambda=1;
N=100000;
e = random('exp',1,N,1);
z = randn(N,1);

for i=1:length(Kappa)
    kappa=Kappa(i);
    R = alpha + ( e*(1-2*kappa)/(kappa*(1-kappa)) + sqrt(2*e/(kappa*(1-kappa))).*z ) / lambda;
    sk(i) = skewness(R);
end
plot(Kappa,sk)