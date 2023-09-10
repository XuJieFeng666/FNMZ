function [x, lambda, iter, iter_in, iter_nt, res, time, Res] =Alg3(A,x0)
m = ndims(A);
n = size(A,1);

I=eye(n);Maxiter=50;eps=10^(-10);
alpha=0.073;sigma=0.005;
% alpha=0.1;sigma=0.4;

x=x0;
o=ttsv(A,x,-2);
p=o*x;
q=p'*x;
F=p-q*x;
J=(m-1)*o-q*I-m*x*p';

iter=0;iter_in=0;iter_nt=0;Res=[ ];
tic
res = norm(F);
while iter<Maxiter && res>eps
    U = null(x'); theta=0.5*norm(F)^2;  JT=J'*F;
    H = U'*J*U;    b = -U'*F;    
%     if cond(H)>1e4
%        cond(H) 
%     end
    q1 = H\b;
    d=U*q1;
    if JT'*d>= 0
        d = -JT;
        iter_nt=iter_nt+1;
    end
    
     for i=0:10
        r=x+(alpha^i).*d;
        x_new=r/norm(r);
        o_new=ttsv(A,x_new,-2);
        p_new=o_new*x_new;
        q_new=p_new'*x_new;
        F_new=p_new-q_new*x_new;
        J_new=(m-1)*o_new-q_new*I-m*x_new*(p_new)';
        theta_new=0.5*norm(F_new)^2;
        %w=theta-2*sigma*(alpha^i)*theta;
        w=theta+sigma*(alpha^i)*JT'*d-(alpha^i)*JT'*(x*x')*d;
        if theta_new<=w
%          fprintf('itr=%.2f, new=%.2f, i=%1d, phi=%.3f, phi_new=%.3f, res=%.3g\n',out.iter,out.new,i, 1/m*q,1/m*q_new, norm(F_new))
            x=x_new;
            o=o_new;
            p=p_new;
            q=q_new;
            F=F_new;
            J=J_new;
            iter_in=iter_in+i;
            break
        end
    end
    res=norm(F);
    Res=[Res res];
    iter=iter+1;
end
time=toc;
lambda=q;
%
%fprintf('sigma=%.4f,alpha=%.4f,itr=%.2f,initr=%.2f, new=%.2f, i=%1d, phi=%.3f, phi_new=%.3f, res=%.3g\n',sigma,alpha,out.iter,out.initer,out.new,i, 1/m*q,1/m*q_new, norm(F_new))

% figure (1)
% 
% niter = length(Res);  err=Res(niter)*ones(1,50);
% err(1:niter)=Res(1:niter);
% iter = 0:1:49;
% semilogy(iter,err,'b-d','LineWidth',2,'MarkerSize',6);
% 
% title('case: (m,n)=(4,50)','fontsize',12)
% xlabel('number of iterations','fontsize',12);
% ylabel('residue','fontsize',12);
% set(gca,'YGrid','on');
%   hold on

end
