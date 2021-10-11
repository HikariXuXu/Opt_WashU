fun = @(x) - 12 * x(2) + 4 * x(1)^2 + 4 * x(2)^2 + 4 * x(1) * x(2);
x0 = [0 0];
[x,fval] = fminunc(fun,x0)