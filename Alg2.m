function [x, lambda, iter, iter_nt, itr_in, res, time, FRes, iter_nt_in] =  Alg2(MT, x0)
m = ndims(MT);
n = size(x0,1);
eta = 10^(-10);maxItr = 300;sigma = 0.01;rho = 0.1;
MT=tensor(MT);
x = x0;
A = ttsv(MT,x,-2);
Ax = A*x;
lambda = x'*Ax;
phi = 1/m * lambda;
Phi = [];
Phi = [Phi lambda];
F = (Ax-lambda*x);
res = norm(F);
FRes = res;
iter = 0;itr_in = 0;iter_nt = 0;iter_nt_in = 0;
tic
while iter<maxItr && res>eta 
    J = (m-1)*A-lambda*eye(n);%-m*x*Ax';
    U = null(x');
    H = U'*J*U;
    b = U'*F;
    q = -(H\b);
%     d = -[J; x']\[F;0];
    d = U*q;
%     d = F;
    flag = 1;
    if F'*d<= 0
        d = F;
        iter_nt=iter_nt+1;
        flag = 0;
    end
    for l=0:10
        alpha = rho^l;
        u = x+alpha*d;
        x1 = u./norm(u); 
        A1 = ttsv(MT,x1,-2);
        Ax1 = A1*x1;
        lambda1 = x1'*Ax1;
        phi1 = 1/m * lambda1;
%             fprintf('itr=%d, isNM=%1d,  l=%d, phi=%3.1g, phi1=%3.1g\n', itr1, isNM, l, phi, phi1);
        if phi1>= phi + alpha*sigma*F'*d
            x = x1;
            A = A1;
            Ax = Ax1;
            lambda = lambda1;
            phi = phi1;
            Phi = [Phi lambda];
            F = (Ax-lambda*x);
            itr_in=itr_in+l;
            if ~flag
               iter_nt_in = iter_nt_in +l;
            end
            break;
        end
    end
    res = norm(F);
    iter = iter + 1;
    FRes = [FRes res]; 
end
time=toc;
end
