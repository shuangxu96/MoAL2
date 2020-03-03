%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo on data with different kinds of noise                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

m = 40;
n = 20;                     %Data size
r = 4;                      %Rank

load('IND.mat')
for kk = 1:30
    randn('seed',kk),RU = randn(m,r);
    randn('seed',kk),RV = randn(r,n);
    X_Ori = RU * RV;        %Original data
    Ind = IND(kk,:);
    p1 = floor(m*n*0.2);
    W = ones(m,n);
    W(Ind(1:p1)) = 0;       %Indicator matrix
    X_Noi = X_Ori;
    X_Noi = W.*X_Noi;       %Add missing components
    
%     U = rand(1,m*n-p1)-0.5+1e-10;
%     lap = -1.5*sign(U).*log(1-2*abs(U));
% slash = sqrt(2)*randn(1,m*n-p1)/rand(1,m*n-p1);
t1=random('t',2,1,m*n-p1);
% g=randn(1,m*n-p1)*5;
% al = randald(m*n-p1,0,1,0.7)';
% lam=0.7; sigma = 3; delta=lam/(1+lam^2); sn = delta*abs(sigma*randn(1,m*n-p1))+sqrt(1-delta^2)*sigma*randn(1,m*n-p1);
%  for i= 1:(m*n-p1); U = rand(1)-0.5+1e-10;mix1(i) = -sign(U).*log(1-2*abs(U));
%      if rand(1)>0.2; mix1(i)=mix1(i)*2;else mix1(i)=mix1(i)*1;end;end
%  for i= 1:(m*n-p1); U = rand(1)-0.5+1e-10;mix2(i) = -sign(U).*log(1-2*abs(U)); in=rand(1);
%      if in<0.2; mix2(i)=mix2(i)*2;elseif in>=0.2&in<0.5; mix2(i)=mix2(i)*1;else mix2(i)=randn(1)*1;end;end
%  for i= 1:(m*n-p1);  in=rand(1);
%      if in<0.2; mix3(i)=randald(1,0,1,0.8)';
%      elseif in>=0.2&in<0.5; U = rand(1)-0.5+1e-10;mix3(i) = -1*sign(U).*log(1-2*abs(U));
%      else mix3(i)=randn(1)*1;end;end

    X_Noi(Ind(p1+1:end)) = X_Noi(Ind(p1+1:end)) + t1;  %Add noise
    
    aa = median(abs(X_Noi(Ind(p1+1:end))));
    aa = sqrt(aa/r);
    
    U0 = rand(m,r)*aa*2-aa;
    V0 = rand(n,r)*aa*2-aa;
    
    %%%%%%%%%%%%%%%%%% MoAL method %%%%%%%%%%%%%%%%%%%%%%%%
    [OutU, OutV, OutW, MoAL, Label, llh] = moal2(W, X_Noi, r, 'MaxIter', 100, 'IniU', U0, 'IniV', V0, 'K', 4, 'NumIter', 50, 'tol', 1e-3);
    Err1 = X_Ori - OutU*OutV';
    
    [A,B] = CWMmc(X_Noi,W,r,[]);
    Err2 = X_Ori - A*B';
    
    [A,B] = CWM_wei(X_Noi,W,r,[]);
    Err3 = X_Ori - A*B';
    
    E1(kk)=sum(sum(abs((Err1))));
    E2(kk)=sum(sum(abs((Err2))));
    E3(kk)=sum(sum(abs((Err3))));
end

format long
[mean(E1),mean(E2),mean(E3)]/m/n


