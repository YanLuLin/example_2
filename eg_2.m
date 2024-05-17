function [A,b,xk] = eg_2(m)
n=m^2;
d=ones(m,1);
e=ones(n,1);
S=spdiags([-d 4*d -d],[-1 0 1],m,m);
M_hat=kron(speye(m,m),S)+spdiags([-e -e],[-m m],n,n);
A=sparse(M_hat+4*speye(n));
x_star=sparse(zeros(n,1)+(-1).^(1:n)');
b=A*x_star-abs(x_star);
xk = sparse(zeros(n,1));
end
