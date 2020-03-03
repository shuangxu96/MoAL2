function [OutU, OutV, OutW, MoAL, Label, llh] = moal2(W, X, r, varargin)
%MOAL2 Perform low-rank matrix factorization with mixture of asymmetric
%Laplacian (MoAL) noise.
%
%   Positional parameters:
%
%     W                An binary matrix indicating whether the entry is missing
%     Y                An observed matrix
%     r                Rank
%   
%   Optional input parameters:  
%
%     'MaxIter'        The maximum number of outer iterations.(Defualt: 30)
%     'K'              The number of asymmetric Laplacian distributions. (Defualt: 4)
%     'NumIter'        The maximum number of inner iterations in M-step. (Defualt: 100)
%     'TuningK'        Binary scalar. Code automatically tunes K if it
%                      takes true, and does not tune otherwise. (Defualt: true)
%     'IniU','IniV'    Initial U and V.
%     'tol'            Converge tolerence. (Defualt: 1e-5)
%
%   Return values:
%     OutU, OutV       Final estimate of U and V.
%     OutW             The learned weights for each entry.
%     Label            The class label for each entry.
%     llh              The likelihood of each iteration.
%     MoAL             MoAL is a struct that contains information of 
%                      mixture of asymmetric Laplacian distribution. MoAL
%                      constains following fields:
%
%       'alpha'        The location parameter
%       'lambda'       The scale parameter                          
%       'kappa'        The skew parameter                          
%       'pi'           The mixing proportion parameter  

%   References: 
%   [1] Shuang Xu, Chunxia Zhang, Jiangshe Zhang, Adaptive Quantile Low-
%       Rank Matrix Factorization. Pattern Recognition, 2019.

%   Written by Shuang Xu (xu.s@outlook.com; shuangxu@stu.xjtu.edu.cn).

% =========================
% parse the input arguments
% =========================
index = find(W(:)~=0);
if nargin > 2
    [MaxIter, K, NumIter, TuningK, IniU, IniV, tol] = parseArgs_moal2(varargin, X, index, r);
end

% =========================
% Initialize the class label
% =========================
R = Initialize_moal2(X(index)',K);
[~,Label(1,:)] = max(R,[],2);
R = R(:,unique(Label)); % delete blank class

% =========================
% Initialize MoAL parameters
% =========================
MoAL.alpha  = zeros(1,K);
MoAL.lambda = rand(1,K);
MoAL.kappa  = rand(1,K);
MoAL.pi     = mean(R,1);

% =========================
% Initial E-step
% =========================
[TempU, TempV] = deal(IniU,IniV);
TempX = IniU*IniV';
Error = X(:) - TempX(:);
Error = Error(index);
t = 1;
[R, llh(t)] = Estep_moal2(Error', MoAL);

% =========================
% main loop
% =========================
flag = false;
while ~flag
    t = t+1;
    old_U = TempU;
    
    % M-E-M-E step
    [MoAL, rho]           = Mstep_Model_moal2(Error, R, MoAL);
    [R, llh(t)]           = Estep_moal2(Error', MoAL);
    [TempW, TempU, TempV] = Mstep_UV_moal2(W,X,TempU,TempV,R,NumIter,rho,MoAL);
    TempX                 = TempU*TempV';
    Error                 = X(:)-TempX(:);
    Error                 = Error(index);
    [R, llh(t)]           = Estep_moal2(Error', MoAL);
    
    % Tuning the number of ALDs, K
    old_K = K;
    if TuningK
        [~,Label(:)] = max(R,[],2);
        u = unique(Label);   % non-empty components
        if size(R,2) ~= size(u,2)
            R = R(:,u);   % remove empty components
            MoAL.alpha  = MoAL.alpha(:,u);
            MoAL.lambda = MoAL.lambda(:,u);
            MoAL.kappa  = MoAL.kappa(:,u);
            MoAL.pi     = MoAL.pi(:,u);
            K = length(u);
            fprintf('The number of ALDs (K) is reduced to %d. \n ', K)
        end
        
    end
    
    % Converge or not
    if t>=MaxIter
        break
    end
    if K~=old_K
        flag = false;
    elseif sum(MoAL.pi)~=1
        flag = false;
    else
        flag = norm(old_U-TempU)<tol;
    end
    
    % display the iteration number
    if numel(X)>500^2
        fprintf('Iteration %d \n', t)
    end
end

OutU = TempU;
OutV = TempV;
OutW = TempW;
end

%% subfunctions - Mstep_UV_moal2
function [TempW, TempU, TempV] = Mstep_UV_moal2(W,X,TempU,TempV,R,NumIter,rho,MoAL)
TempW = zeros(size(X));
lambda = MoAL.lambda;

TempW(W(:)~=0) = sum(R.*rho.*lambda,2);
TempW = TempW./ (abs(X-TempU*TempV')+1e-6); %for numerical stability

pa = [];
pa.IniU = TempU;
pa.IniV = TempV;
pa.maxiter = NumIter;
pa.display = 0;
[TempU2,TempV2] = EfficientMCL2(X, sqrt(TempW), TempU,TempV, NumIter, 0.00000001);
if ~sum(isnan(TempU2(:))) && ~sum(isnan(TempV2(:)))
    TempU = TempU2;
    TempV = TempV2;
else
    TempU = TempU+eps;
    TempV = TempV+eps;
end
end

%% subfunctions - Mstep_Model_moal2
function [MoAL,rho] = Mstep_Model_moal2(Error,R,MoAL)
kappa = MoAL.kappa;

% pi
nk = sum(R,1);
pi = nk/size(R,1);

% lambda
sign = Error<0;
rho = sign*(1-kappa) + ~sign*kappa;
lambda = nk./ sum(R.*rho.*abs(Error));
indINF = find(~isfinite(lambda));
if ~isempty(indINF)
    lambda(indINF) = nk(indINF)/1e-6;
end
% kappa
eta = sum(R.*Error).*lambda;
sDelta = sqrt(4*nk.^2+eta.^2);
kappa = (2*nk + eta - sDelta) ./ (2*eta);

MoAL.lambda = lambda;
MoAL.kappa  = kappa;
MoAL.pi     = pi;
end

%% subfunctions - Estep_moal2
function [R, llh] = Estep_moal2(Error, MoAL)
% Compute the conditional probability
alpha  = MoAL.alpha;
lambda = MoAL.lambda;
kappa  = MoAL.kappa;
pi     = MoAL.pi;

n = size(Error,2);
k = size(alpha,2);
logRho = zeros(n,k);

for i = 1:k
    logRho(:,i) = logpdfald(Error,alpha(i),lambda(i),kappa(i));
end
logRho = bsxfun(@plus,logRho,log(pi));
T = logsumexp(logRho,2);
llh = sum(T)/n; % loglikelihood
logR = bsxfun(@minus,logRho,T);
R = exp(logR);
end

%% subfunctions - initialize_moal2
function R = Initialize_moal2(X, k)
[~,n] = size(X);
idx = randsample(n,k);
m = X(:,idx);
[~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
[u,~,label] = unique(label);
while k ~= length(u)
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    [u,~,label] = unique(label);
end
R = full(sparse(1:n,label,1,n,k,n));
end

%% subfunctions - parseArgs_moal2
function [MaxIter, K, NumIter, TuningK, IniU, IniV, tol] = parseArgs_moal2(vararginCell, X, index, r)
[vararginCell{:}] = convertStringsToChars(vararginCell{:});

if sum(strcmp(vararginCell, 'IniU'))==0
    IniU = InitializeComponent(X, index, size(X,1), r);
else
    IniU = [];
end
if sum(strcmp(vararginCell, 'IniV'))==0
    IniV = InitializeComponent(X, index, size(X,2), r);
else
    IniV = [];
end

pnames = { 'MaxIter' 'K' 'NumIter' 'TuningK' 'IniU' 'IniV' 'tol'};
dflts  = {     30,   4,    100,      true,    IniU,   IniV,  1e-5};
[MaxIter, K, NumIter, TuningK, IniU, IniV, tol] ...
    = internal.stats.parseArgs(pnames, dflts, vararginCell{:});
end

%% subfunctions - InitializeComponent
function IniU = InitializeComponent(X, index, shape1, shape2)
s = median(abs(X(index)));
s = sqrt(s/shape2);
if min(X(index)) >= 0
    IniU = rand(shape1, shape2)*s;
else
    IniU = rand(shape1, shape2)*s*2-s;
end
end

%% subfunctions - logsumexp
function s = logsumexp(x, dim)
% Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
%   By default dim = 1 (columns).
% Written by Michael Chen (sth4nth@gmail.com).
if nargin == 1
    % Determine which dimension sum will use
    dim = find(size(x)~=1,1);
    if isempty(dim), dim = 1; end
end

% subtract the largest in each column
y = max(x,[],dim);
x = bsxfun(@minus,x,y);
s = y + log(sum(exp(x),dim));
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end
end

%% subfunctions - logpdfald
function logy = logpdfald(x,alpha,lambda,kappa)
x = x(:);
xa = x-alpha;
xabs = abs(xa);
sign = (xa >= 0);
al_kernal = lambda*xabs.*( kappa*sign + (1-kappa)*~sign );
logy = log(lambda*kappa*(1-kappa)) - al_kernal;
end