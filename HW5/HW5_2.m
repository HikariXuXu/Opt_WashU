%% Newton's methods
clear;clc;close all
x0 = 10;
fun = @(x) exp(-x^2)-x;

%% ======================= Part 1==========================================

%**************************************************************************
%Compute the gradient here:
grad_fun = @(x) -2*x*exp(-x^2)-1;
%**************************************************************************


%% ======================= Part 2==========================================
% Suppose the algorithm stops after n tierations.
% Variables: X -- 1-by-n matrix, storing (x_1) at each iteration;
%            D -- 1-by-n matrix, stroing the descent direction at each
%                 iteration;
%            F -- 1-by-n vector, storing the value of objective function
%                 at each iteration;
%            grad_F -- 1-by-n vectore, storing the value of gradient
%                      function at each iteration.
X=[]; D=[]; F=[]; grad_F = [];
X(:,1) = x0;
count = 1;
F(1,1) = fun(x0);
grad_F(:,1) = grad_fun(x0);
D(:,1) = -fun(x0)/grad_fun(x0);


%**************************************************************************
% Write the update rule of Newton's method here:
for k=1:1000
    if norm(D(:,k)) < 10^(-5)
        break;
    end
    X(:,k+1) = X(:,k) + D(:,k);
    count = count + 1;
    F(1,k+1) = fun(X(:,k+1));
    grad_F(:,k+1) = grad_fun(X(:,k+1));
    D(:,k+1) = -fun(X(:,k+1))/grad_fun(X(:,k+1));
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
%**************************************************************************


%% ======================= Part 3==========================================
opt = optimoptions(@fminunc,'Display','iter-detailed','Algorithm','quasi-newton');

%**************************************************************************
% Use fminunc to minize the objective function here:
int_fun = @(x) -sqrt(pi)/2*erf(x)+x^2/2
[x,fval] = fminunc(int_fun,x0,opt)
%**************************************************************************


%% ======================= Part 4==========================================
%**************************************************************************
% Repeat the above process to minimize f(x) = 7x - log x and write the the
% results in your report.
clear;close all
x0 = [1 0 0.01 0.1];
fun = @(x) 7*x-log(x);
%Compute the gradient here:
grad_fun = @(x) 7-1/x;
%Compute the Heissian here:
heissian_fun = @(x) 1/x^2;
% Suppose the algorithm stops after n tierations.
% Variables: X -- 1-by-n matrix, storing (x_1) at each iteration;
%            D -- 1-by-n matrix, stroing the descent direction at each
%                 iteration;
%            F -- 1-by-n vector, storing the value of objective function
%                 at each iteration;
%            grad_F -- 1-by-n vectore, storing the value of gradient
%                      function at each iteration.
X=[]; 
for i=1:4
    D=[]; F=[]; grad_F = [];
    X(i,1) = x0(i);
    count = 1;
    F(1,1) = fun(X(i,1));
    grad_F(:,1) = grad_fun(X(i,1));
    D(:,1) = -grad_fun(X(i,1))/heissian_fun(X(i,1));

    %Newton's method:
    for k=1:1000
        if norm(D(:,k)) < 10^(-5)
            break;
        end
        X(i,k+1) = X(i,k) + D(:,k);
        count = count + 1;
        F(1,k+1) = fun(X(i,k+1));
        grad_F(:,k+1) = grad_fun(X(i,k+1));
        D(:,k+1) = -grad_fun(X(i,k+1))/heissian_fun(X(i,k+1));
    end
end
% Make table
col_name = ["k","x^k1","x^k2","x^k3","x^k4"];
T = table([1:8]',X(1,1:8)',X(2,1:8)',X(3,1:8)',X(4,1:8)','VariableNames',col_name);
%**************************************************************************
