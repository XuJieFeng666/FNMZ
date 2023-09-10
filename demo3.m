clear;clc
warning off
% example 3
m=3;
N = 10:10:80;
results = cell(length(N),1);
for ni=1:length(N)
    n = N(ni);
    a = (-1).^(1:n)./(1:n);
    a = a'; b = reshape(a, [1 n]); c = reshape(a, [1 1 n]);
    A = a+b+c;   
    maxnum = max(abs(A(:))); 
    A= A/maxnum; 
    A =tensor(A);
    %%
    kk = 100;
    tt = 2;
    x = cell(kk,tt); lambda = zeros(kk,tt); iter = zeros(kk,tt); iter_nt = zeros(kk,tt);
    iter_in = zeros(kk,tt); time = zeros(kk,tt); res = zeros(kk,tt); Res = cell(kk,tt);
    result= zeros(tt,6);
    for k=1:kk
        y=rand(n,1);x0=y/norm(y);
        t=1;
            [x{k,t}, lambda(k,t), iter(k,t), iter_nt(k,t), iter_in(k,t), res(k,t), time(k,t), Res{k,t}] =  Alg2(A, x0);
        t=2;
            [x{k,t}, lambda(k,t), iter(k,t), iter_in(k,t), iter_nt(k,t), res(k,t), time(k,t), Res{k,t}] =  Alg3(A, x0);
    end
    fprintf('%3d ', n)
    for t=1:tt
        idx = res(:,t)<1e-10;
%          plot(log10(Res{40,1}))
        a = mean(iter(idx,t)); b = mean(iter_nt(idx,t));
        c = mean(iter_in(idx,t)); d = mean(time(idx,t));
        e = mean(res(idx,t)); suc = sum(idx);
        result(t, :) = [a b c d e suc];
        fprintf('& %-.1f / %-.1f / %-.1f / %-.4f / %-.2d / %3d',...
            a,b,c,d,e,suc)
    end
    results{ni} = result;
    fprintf('\\\\ \n')
end
