function [alpha, F, X] = steepestdescent(fun, grad_fun, eps, x0, choice, stp_eps)
% Write the script 'steepestdescent.m' on your own
% Function: steepestdescent
% Inputs: fun -- the objective function;
%        grad_fun -- the gradient of the objective function;
%        eps -- error tolerance for line search
%        x0 -- initial condition
%        choice -- choice of line search algorithm
%        stp_eps -- error tolerance for steepest descend
% Output: Suppose that the steepest descent has run for n steps to stop
%        alpha: 1-by-n vector, storing the value of \alpha at each step of
%               steepest descent;
%        F: 1-by-n vecotr, stroing the value of the objective function at
%           each step of steepest descent;
%        X: 2-by-n matrix, stroing the value of (x_1, x_2)' at each step of
%           steepest descent.
X(:,1) = x0; %initialize state
D(:,1) = -grad_fun(X(:,1));%initial direction
F(1,1) = fun(X(:,1));
if choice == 1
    for k=1:1000 %allow 1000 steps to converge
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

        while 1
            alpha_mid = (alpha_L+alpha_U)/2;
            if abs(grad_fun(X(:,k)+alpha_mid*D(:,k))'*D(:,k)) < eps
                alpha(1,k) = alpha_mid;
                break;
            else
                if grad_fun(X(:,k)+alpha_mid*D(:,k))'*D(:,k) > 0
                    alpha_U = alpha_mid;
                else
                    alpha_L = alpha_mid;
                end
            end
        end
        X(:,k+1) = X(:,k)+alpha(1,k)*D(:,k);
        D(:,k+1) = -grad_fun(X(:,k+1));
        F(1,k+1) = fun(X(:,k+1));
        if norm(D(:,k+1)) < stp_eps
            break;
        end
    end
else
    epsilon = 0.5;
    sigma = 2;
    for k=1:1000 %allow 1000 steps to converge
        %initial alpha
        if D(1,k) < 0
            alpha_list(1) = -X(1,k)/D(1,k);
        else
            alpha_list(1) = 1/0;
        end
        if D(2,k) < 0
            alpha_list(2) = -X(2,k)/D(2,k);
        else
            alpha_list(2) = 1/0;
        end
        if D(1,k) + D(2,k) > 0
            alpha_list(3) = (100-(X(1,k)+X(2,k)))/(D(1,k)+D(2,k));
        else
            alpha_list(3) = 1/0;
        end
        if D(1,k) - D(2,k) > 0
            alpha_list(4) = (50-(X(1,k)-X(2,k)))/(D(1,k)-D(2,k));
        else
            alpha_list(4) = 1/0;
        end
        alpha0 = min(alpha_list);
        if alpha0 == 1/0
            disp("error");
        end
        alpha_bar = alpha0;
        while F(1,k)-fun(X(:,k)+alpha_bar*D(:,k)) < -alpha_bar*epsilon*grad_fun(X(:,k))'*D(:,k)
            alpha_bar = alpha_bar/sigma;
        end
        alpha(1,k) = alpha_bar;
        X(:,k+1) = X(:,k)+alpha(1,k)*D(:,k);
        D(:,k+1) = -grad_fun(X(:,k+1));
        F(1,k+1) = fun(X(:,k+1));
        if norm(D(:,k+1)) < stp_eps
            break;
        end
    end
end
end

