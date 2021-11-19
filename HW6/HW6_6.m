eps = 10^(-6);
Q = [35 19 22 28 16 3 16 6 4 4
    19 43 33 19 5 2 5 4 0 0
    22 33 40 29 12 7 6 2 2 4
    28 19 29 39 16 7 14 6 2 4
    16 5 12 16 12 4 8 2 4 8
    3 2 7 7 4 5 1 0 1 4
    16 5 6 14 8 1 12 2 2 4
    6 4 2 6 2 0 2 4 0 0
    4 0 2 2 4 1 2 0 2 4
    4 0 4 4 8 4 4 0 4 16];
c = [-1; 0; 0; -3; 0; -2; 0; -6; -7; -4];
fun = @(x) 1/2*x'*Q*x-c'*x;
grad_fun = @(x) Q*x-c;
% Conjugate gradient method
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
d0 = -grad_fun(x0);
X = []; F = []; D = []; A=[];
X(:,1) = x0;
F(:,1) = fun(x0);
D(:,1) = d0;
for k=1:10
    if norm(grad_fun(X(:,k))) <= eps
        break;
    else
        A(:,k) = -(grad_fun(X(:,k))'*D(:,k))/(D(:,k)'*Q*D(:,k));
        X(:,k+1) = X(:,k) + A(:,k)*D(:,k);
        F(:,k+1) = fun(X(:,k+1));
        D(:,k+1) = -grad_fun(X(:,k+1))+(grad_fun(X(:,k+1))'*Q*D(:,k))/(D(:,k)'*Q*D(:,k))*D(:,k);
    end
end
F_minus = F-F(:,11);
Rate = [];
for i=2:11
    Rate(:,i-1) = F_minus(:,i)/F_minus(:,i-1);
end
