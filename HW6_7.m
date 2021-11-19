eps = 10^(-6);
cgm_eps = 10^(-4);
theta = 10;
x0 = [10; 20];
fun = @(x) -9*x(1) -10*x(2) + theta*(-log(100-x(1)-x(2))-log(x(1))-log(x(2)) - log(50-x(1)+x(2)));
grad_fun = @(x)[-9 + theta*(1/(100-x(1)-x(2))-1/x(1)+1/(50-x(1)+x(2))); -10 + theta*(1/(100-x(1)-x(2))-1/x(2)-1/(50-x(1)+x(2)))];
% Conjugate gradient method
d0 = -grad_fun(x0);
X = []; F = []; D = []; A=[];
X(:,1) = x0;
F(:,1) = fun(x0);
D(:,1) = d0;
for k=1:1000
    if norm(grad_fun(X(:,k))) <= cgm_eps
        break;
    else
        alpha_L = 0;
        %find alpha_U
        if D(1,k) < 0
            alpha_U_list(1) = -X(1,k)/D(1,k);
        else
            alpha_U_list(1) = 1/0;
        end
        if D(2,k) < 0
            alpha_U_list(2) = -X(2,k)/D(2,k);
        else
            alpha_U_list(2) = 1/0;
        end
        if D(1,k) + D(2,k) > 0
            alpha_U_list(3) = (100-(X(1,k)+X(2,k)))/(D(1,k)+D(2,k));
        else
            alpha_U_list(3) = 1/0;
        end
        if D(1,k) - D(2,k) > 0
            alpha_U_list(4) = (50-(X(1,k)-X(2,k)))/(D(1,k)-D(2,k));
        else
            alpha_U_list(4) = 1/0;
        end
        alpha_U = min(alpha_U_list);
        if alpha_U == 1/0
            disp("error");
        end
        
        while alpha_L ~= alpha_U
            alpha_mid = (alpha_L+alpha_U)/2;
            if abs(grad_fun(X(:,k)+alpha_mid*D(:,k))'*D(:,k)) < eps
                break;
            else
                if grad_fun(X(:,k)+alpha_mid*D(:,k))'*D(:,k) > 0
                    alpha_U = alpha_mid;
                else
                    alpha_L = alpha_mid;
                end
            end
        end
        A(:,k) = min([alpha_mid,alpha_U]);
        X(:,k+1) = X(:,k) + A(:,k)*D(:,k);
        F(:,k+1) = fun(X(:,k+1));
        D(:,k+1) = -grad_fun(X(:,k+1))+(grad_fun(X(:,k+1))'*grad_fun(X(:,k+1)))/(grad_fun(X(:,k))'*grad_fun(X(:,k)))*D(:,k);
    end
end
