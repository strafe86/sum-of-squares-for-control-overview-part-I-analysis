% app33.m
clear; clc;

gamma_low  = 0;
gamma_high = 20;
tol        = 1e-3;

gamma_best = gamma_low;
%реализация бисекции по gamma
while (gamma_high - gamma_low) > tol
    gamma = 0.5*(gamma_low + gamma_high);
    [ok, V_cur, beta_cur] = is_feasible_gamma(gamma);
    if ok
        gamma_low  = gamma;
        gamma_best = gamma;
        V_best     = V_cur;
        beta_best  = beta_cur;
    else
        gamma_high = gamma;
    end
    fprintf('gamma_low = %.6f, gamma_high = %.6f, ok = %d\n', ...
            gamma_low, gamma_high, ok);
end

gamma = gamma_best
Vn    = V_best
beta  = beta_best