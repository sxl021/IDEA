clear;clc
%IDEA framework
%IDEA implements the interval design  for enhanced accuracy (IDEA) algorithm provided in 
%Yuanbo, Cheng., et al. "Interval Design for Signal Parameter Estimation
%from Quantized Data"
%Also, see "The Cram¨¦r¨CRao Bound for Signal Parameter Estimation From
%Quantized Data" in IEEE Signal Processing Magazine, 39(1), 118-125.
%Written by Xiaolei Shang, April,2022.

Umax=5;     % The considered maximum of the threshold.
K=512+1;    % The number of grid
grid=linspace(-Umax,Umax,K).';  % Grid set
N=3;        % The number of quantization bit
Len=2^N -1;
f=@(u,l)(exp(-u.^2/2)-exp(-l.^2/2)).^2./(normcdf(-u,'upper')-normcdf(-l,'upper'))/2/pi;

val=f(grid,-inf);               %Compute and store DP1(u_1)
ind_opt=zeros(K,Len-1);         %Optimal index
ind_tmp=ind_opt;
for k=2:Len
    J=zeros(K,1);
    for m=k:K
        u=grid(m);
        index=grid < u;
        L=sum(index);
        tmp=f(u*ones(L,1),grid(index))+val(index);
        [J(m),I]=max(tmp);
        clear tmp
        ind_opt(m,k-1)=I;
        ind_opt(m,1:k-2)=ind_tmp(I,1:k-2);
    end
    val=J;
    ind_tmp=ind_opt;
end

% Compute DP_A(u_A)
u=inf;
index=grid < u;
L=sum(index);
tmp=f(u*ones(L,1),grid(index))+val(index);
[J,I]=max(tmp);

% Optimal interval
int_opt=[-inf; grid([ind_opt(I,1:end) I]); inf];
J_opt=sum(f(int_opt(2:end),int_opt(1:end-1)));


