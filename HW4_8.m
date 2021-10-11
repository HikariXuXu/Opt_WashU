clc;clear all;close all
%% ======================= Part 1 =========================================
Q=[10 4; 4 2]; c=[14;6]; a=20; %define objective function
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

X(:,1)=[0;10]; %initialize state
D(:,1)=-Q*X(:,1)+c; %initial direction
A(1,1)=D(:,1)'*D(:,1)/(D(:,1)'*Q*D(:,1)); %initial step size
F(1,1)=X(:,1)'*Q*X(:,1)/2 - c'*X(:,1) + a;
for k=1:1000 %allow 1000 steps to converge
    if(norm(D(:,k))<ep)
        break; 
    end %quit if tolerance reached
%-------------------------------------------------------------------------
%*************************************************************************
% Write the update rule of steepest descent on your own
    X(:,k+1)=X(:,k)+A(1,k)*D(:,k); %update state
    D(:,k+1)=-Q*X(:,k+1)+c; %update direction
    A(1,k+1)=D(:,k+1)'*D(:,k+1)/(D(:,k+1)'*Q*D(:,k+1)); %update step size
    F(1,k+1)=X(:,k+1)'*Q*X(:,k+1)/2 - c'*X(:,k+1) + a;
%*************************************************************************
%-------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
%Contour plot
[xx,yy]=meshgrid(-10:.2:10);
zz=(Q(1,1)*xx.^2+(Q(2,1)+Q(1,2))*xx.*yy+Q(2,2)*yy.^2)/2-c(1)*xx-c(2)*yy+a;

figure
[C,h]=contour(xx,yy,zz,400); hold on, plot(X(1,:),X(2,:),'k',X(1,:),X(2,:),'r*')
hold off, xlabel('x'), ylabel('y')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')), colormap cool
set(h,'TextStep',get(h,'LevelStep')), colormap jet
axis([-5 5 -2 10])

figure
mesh(xx,yy,zz)
zlim([0 1200])
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
fun = @(x) 5*x(1)^2 + x(2)^2 +4*x(1)*x(2) - 14*x(1) -6*x(2) + 20;
%**************************************************************************
x0 = [0 10]';
%**************************************************************************
% Use fminunc to minize the objective function here:
[x,fval] = fminunc(fun,x0,opt)
%**************************************************************************


%% ======================= Part 3 =========================================
Q=[20 5; 5 2]; c=[14;6]; a=20; %define objective function
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

X(:,1)=[40;-100]; %initialize state
D(:,1)=-Q*X(:,1)+c; %initial direction
A(1,1)=D(:,1)'*D(:,1)/(D(:,1)'*Q*D(:,1)); %initial step size
F(1,1)=X(:,1)'*Q*X(:,1)/2 - c'*X(:,1) + a;
for k=1:1000 %allow 1000 steps to converge
    if(norm(D(:,k))<ep)
        break; 
    end %quit if tolerance reached
%-------------------------------------------------------------------------
%*************************************************************************
% Write the update rule of steepest descent on your own
    X(:,k+1)=X(:,k)+A(1,k)*D(:,k); %update state
    D(:,k+1)=-Q*X(:,k+1)+c; %update direction
    A(1,k+1)=D(:,k+1)'*D(:,k+1)/(D(:,k+1)'*Q*D(:,k+1)); %update step size
    F(1,k+1)=X(:,k+1)'*Q*X(:,k+1)/2 - c'*X(:,k+1) + a;
%*************************************************************************
%-------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
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
new_T = [new_T(1:10,:);new_T(height(new_T)-9:height(new_T),:)];
new_T


%% ======================= Part 4 =========================================
Q=[20 5; 5 16]; c=[14;6]; a=20; %define objective function
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

X(:,1)=[40;-100]; %initialize state
D(:,1)=-Q*X(:,1)+c; %initial direction
A(1,1)=D(:,1)'*D(:,1)/(D(:,1)'*Q*D(:,1)); %initial step size
F(1,1)=X(:,1)'*Q*X(:,1)/2 - c'*X(:,1) + a;
for k=1:1000 %allow 1000 steps to converge
    if(norm(D(:,k))<ep)
        break; 
    end %quit if tolerance reached
%-------------------------------------------------------------------------
%*************************************************************************
% Write the update rule of steepest descent on your own
    X(:,k+1)=X(:,k)+A(1,k)*D(:,k); %update state
    D(:,k+1)=-Q*X(:,k+1)+c; %update direction
    A(1,k+1)=D(:,k+1)'*D(:,k+1)/(D(:,k+1)'*Q*D(:,k+1)); %update step size
    F(1,k+1)=X(:,k+1)'*Q*X(:,k+1)/2 - c'*X(:,k+1) + a;
%*************************************************************************
%-------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
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
