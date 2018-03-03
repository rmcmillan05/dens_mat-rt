function [t, y] = rk_step(func,t_n,y_n)
global step

k1 = step * feval(func, t_n, y_n);
k2 = step * feval(func, t_n+0.5*step, y_n+0.5*k1);
k3 = step * feval(func, t_n+0.5*step, y_n+0.5*k2);
k4 = step * feval(func, t_n+step, y_n+k3);

y = y_n + 1/6 * (k1 +2*k2 + 2*k3 + k4);
t = t_n + step;