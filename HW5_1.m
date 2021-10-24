clear
close all
tic % start counting time
theta = 10; % change theta here
eps = 10^(-6); % error tolerance for line search
stp_eps = 10^(-4); % error tolerance for steepest descend
opt = optimoptions(@fminunc,'Display','iter');
fun = @(x) -9*x(1) -10*x(2) + theta*(-log(100-x(1)-x(2))-log(x(1))-log(x(2)) - log(50-x(1)+x(2)));

%% ======================= Part 1==========================================

%**************************************************************************
%Compute the gradient here:
grad_fun = @(x)[-9 + theta*(1/(100-x(1)-x(2))-1/x(1)+1/(50-x(1)+x(2))); -10 + theta*(1/(100-x(1)-x(2))-1/x(2)-1/(50-x(1)+x(2)))];
%**************************************************************************




%% ======================= Part 2==========================================
%**************************************************************************

x0 = [8 90]'; % change this for different sub-problems in the homework
choice = 1; % 1 = bisection; 2 = Armijo


[alpha, F, X]= steepestdescent(fun,grad_fun,eps,x0,choice,stp_eps);
%**************************************************************************




%% ======================= Part 3==========================================
%**************************************************************************
% Compute the convergence constant
convconst = (F(1,size(F,2)-1)-F(1,size(F,2)))/(F(1,size(F,2)-2)-F(1,size(F,2)));
%**************************************************************************
time = toc % recording the execution time since `tic'


%% ======================= Part 4==========================================
% In your report, you should compare the number of iterations to achieve 
% the same levle of precision, the overal time the algorithm takes, and 
% the amount of time on each iteration, when choosing different line search 
% methods.
%==========================================================================
