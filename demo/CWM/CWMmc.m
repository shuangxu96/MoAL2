% The cyclic weighted median algorithm for L1-norm matrix factorization %%%%%
% Solve
% min_{U,V}  ||W.*(X-UV^T)||_1 .
% Please cite: "Deyu Meng, Zongben Xu, Lei Zhang, Ji Zhao. A cyclic
% weighted median method for L1 low-rank matrix factorization with missing entries. AAAI 2013."

% Written by: Deyu Meng, Xi'an Jiaotong University
% Email: dymeng@mail.xjtu.edu.cn
% Homepage: http://dymeng.gr.xjtu.edu.cn

function [U,V] = CWMmc(X,W,r,param)
%%%%%%%%%%%%%   INPUT   %%%%%%%%%%%%%
% X - Input data matrix (d*n)
% W - Indicator matrix (d*n)
% r - Rank of the output matrix
% param.IniU,param.IniV - Initialization of U (d*r) and V (n*r)
% param.maxiter - Maximum number of iterations
%           default = 100
% param.tol - Tolerance for subspace updating
%           default = 1e-10
%param.display - display the iterative process (1 yes; 2 no)
%           default = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%
% U,V - Output U (d*r) and V (n*r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[d n] = size(X);

ind = find(W(:)~=0);
if (~isfield(param,'IniU'))
    MedValX = median(abs(X(ind)));
    MedValX = sqrt(MedValX/r);
    IniU = rand(d,r)*MedValX*2 - MedValX;
else
    IniU = param.IniU;
end
if (~isfield(param,'IniV'))
    MedValX = median(abs(X(ind)));
    MedValX = sqrt(MedValX/r);
    IniV = rand(n,r)*MedValX*2 - MedValX;
else
    IniV = param.IniV;
end
if (~isfield(param,'maxiter'))
    maxiter = 100;
else
    maxiter = param.maxiter;
end
if (~isfield(param,'tol'))
    tol = 1e-10;
else
    tol = param.tol;
end
if (~isfield(param,'display'))
    display = 1;
else
    display = param.display;
end

U = IniU;
V = IniV;
t = sum(sum(abs(X)));
for i = 1:maxiter
    Order = randperm(r);
    for j = Order
        TempX = X - U*V' + U(:,j)*V(:,j)';
        V(:,j) = OPTmc(TempX, W,U(:,j),j);
        U(:,j) = OPTmc(TempX',W',V(:,j),j);
    end
	W_bin = double(W~=0);
    t = [t sum(sum(abs(W_bin.*(X - U*V'))))];
    if display
        disp(['The L1 objective in ',num2str(i),' step is ',num2str(t(end))]);
    end
    if norm(t(end)-t(end-1))<tol
    %if norm(IniU - U) < tol
        break;
    else
        IniU = U;
    end
end

if display
    disp(['The CWM algorithm terminates in ',num2str(i),' steps.']);
end
Nu = sqrt(sum(U.^2))';
Nv = sqrt(sum(V.^2))';
No = diag(Nu.*Nv);
U = U*diag(1./Nu)*sqrt(No);
V = V*diag(1./Nv)*sqrt(No);
