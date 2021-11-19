clear;clc;close all
x0 = [2, 2]';
fun = @(x) 1/2*(x(1)^2-x(2))^2+1/2*(1-x(1))^2;
%Compute the gradient here:
grad_fun = @(x) [2*x(1)*(x(1)^2-x(2))+1-x(1); -(x(1)^2-x(2))];
%Compute the Heissian here:
heissian_fun = @(x) [6*x(1)^2-2*x(2)-1 -2*x(1); -2*x(1) 1];
% Initials X, D, F, grad_F.
X=[]; D=[]; F=[]; grad_F = [];
X(:,1) = x0;
F(1,1) = fun(x0);
grad_F(:,1) = grad_fun(x0);
D(:,1) = -inv(heissian_fun(x0))*grad_fun(x0);

%Newton's method:
for k=1:1000
    if norm(D(:,k)) < 10^(-5)
        break;
    end
    X(:,k+1) = X(:,k) + D(:,k);
    F(1,k+1) = fun(X(:,k+1));
    grad_F(:,k+1) = grad_fun(X(:,k+1));
    D(:,k+1) = -inv(heissian_fun(X(:,k+1)))*grad_fun(X(:,k+1));
end

% Make table
col_name = ["k","x^k","f(x^k)","f'(x^k)","f(x^k)/f'(x^k)"];
T = table([1:length(F)]',X(:,:)',F',grad_F',F'./grad_F','VariableNames',col_name);
% set desired precision in terms of the number of decimal places
n_decimal = 6;
% create a new table
new_T = varfun(@(x) num2str(x, ['%' sprintf('.%df', n_decimal)]), T);
% preserve the variable names and the row names in the original table
new_T.Properties.VariableNames = T.Properties.VariableNames;
new_T = [T(:,1),new_T(:,2:5)];
new_T
