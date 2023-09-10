clear; clc
% example 1
eta = 1e-10;
warning off
lam_star = 3.1754;
m=4;n=2; 
A = zeros(2,2,2,2);
A(1,1,1,1) = 4/sqrt(3);A(2,2,2,2)=4/sqrt(3);
A(1,1,1,2)=1;A(1,1,2,1)=1;A(1,2,1,1)=1;A(2,1,1,1)=1;
A(1,2,2,2)=1;A(2,1,2,2)=1;A(2,2,1,2)=1;A(2,2,2,1)=1;
A = tensor(A);

%%
kk = 100;
tt = 2;
x = cell(kk,tt); lambda = zeros(kk,tt); iter = zeros(kk,tt); iter_nt = zeros(kk,tt);
iter_in = zeros(kk,tt); time = zeros(kk,tt); res = zeros(kk,tt); Res = cell(kk,tt);
for k=1:kk
    y=rand(n,1);x0=y/norm(y);
    t=1;
    [x{k,t}, lambda(k,t), iter(k,t), iter_nt(k,t), iter_in(k,t), res(k,t), time(k,t), Res{k,t}]...
        =  Alg2(A,x0);
    t=2;
    [x{k,t}, lambda(k,t), iter(k,t), iter_nt(k,t), iter_in(k,t), res(k,t), time(k,t), Res{k,t}]...
        =  Alg3(A,x0);
end
tt = t;
for t=1:tt
    idx_lam = abs(lambda(:,t) - lam_star)<1e-4;
    occ = sum(idx_lam);
    idx = res(:,t)<eta;
%          plot(log10(Res{40,1}))
    a = mean(iter(idx,t)); b = mean(iter_nt(idx,t));
    c = mean(iter_in(idx,t)); d = mean(time(idx,t));
    e = mean(res(idx,t)); suc = sum(idx);
     fprintf('& %-.1f & %-4.1f & %-5.1f & %-8.4f & %-2.2d & %3d & %3d\\\\\n',...
     a,b,c,d,e,suc,occ)
end

