clc;clear all;close all
%% ======================= Part 1 =========================================
fun = @(X) ( X(1)-2 )^4 + ( X(1)-2*X(2) )^2; %define objective function
ep=10^(-6); %tolerance
%-------------------------------------------------------------------------
% Suppose the algorithm stops after n tierations.
% Variables: X -- 2-by-n matrix, storing (x_1, x_2) at each iteration;
%            D -- 2-by-n matrix, stroing the descent direction at each
%                 iteration;
%            A -- 1-by-n vector, storing \alpha at each iteration;
%            F -- 1-by-n vector, storing the value of objective function
%                 at each iteration.
X=[]; D=[]; A=[]; F=[]; %declare state, direction, step, function value
%-------------------------------------------------------------------------

%**************************************************************************
%Compute the gradient here:
grad_fun = @(X) [4*(X(1)-2)^3+2*(X(1)-2*X(2)); -4*(X(1)-2*X(2))];

%**************************************************************************

%**************************************************************************
%Compute the step size here:
step_fun = @(X,D) D'*D/(D'*[12*(X(1)-2)^2+2 -4;-4 8]*D);

%**************************************************************************

X(:,1)=[0;3]; %initialize state
D(:,1)=-grad_fun(X(:,1)); %initial direction
A(1,1)=step_fun(X(:,1),D(:,1)); %initial step size
F(1,1)=fun(X(:,1));
for k=1:1000 %allow 1000 steps to converge
    if(norm(D(:,k))<ep) 
        break; 
    end %quit if tolerance reached
%-------------------------------------------------------------------------
%*************************************************************************
% Write the update rule of steepest descent on your own
    X(:,k+1)=X(:,k)+A(1,k)*D(:,k); %update state
    D(:,k+1)=-grad_fun(X(:,k+1)); %update direction
    A(1,k+1)=step_fun(X(:,k+1),D(:,k+1)); %update step size
    F(1,k+1)=fun(X(:,k+1));
%*************************************************************************
%-------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
%Contour plot
[xx,yy]=meshgrid(-1:.2:5);
zz=( xx-2 ).^4 + ( xx-2*yy ).^2;

figure
[C,h]=contour(xx,yy,zz,100); hold on, plot(X(1,:),X(2,:),'k',X(1,:),X(2,:),'r*')
hold off, xlabel('x'), ylabel('y')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')), colormap cool
set(h,'TextStep',get(h,'LevelStep')), colormap jet

figure
mesh(xx,yy,zz)
zlim([0 200])
grid on

% Make table
col_name = {'k','x_1^k','x_2^k','d_1^k','d_2^k','||d^k||_2','alpha^k','f(x^k)'};
T = table([1:length(F)]',X(1,:)',X(2,:)',D(1,:)',D(2,:)',...
    sqrt(D(1,:)'.^2+D(2,:)'.^2),A',F','VariableNames',col_name);
% set desired precision in terms of the number of decimal places
n_decimal = 6;
% create a new table
new_T = varfun(@(x) num2str(x, ['%' sprintf('.%df', n_decimal)]), T);
% preserve the variable names and the row names in the original table
new_T.Properties.VariableNames = T.Properties.VariableNames;
new_T = [T(:,1),new_T(:,2:8)];
new_T

%% ======================= Part 2 =========================================
opt = optimoptions(@fminunc,'Display','iter-detailed','Algorithm',...
    'quasi-newton','HessUpdate','steepdesc','MaxFunctionEvaluations',2000);

%**************************************************************************
% Define the objective function here:
fun = @(x) ( x(1)-2 )^4 + ( x(1)-2*x(2) )^2;
%**************************************************************************
x0 = [0 3]';

%**************************************************************************
% Use fminunc to minize the objective function here:
[x,fval] = fminunc(fun,x0,opt)
%**************************************************************************


