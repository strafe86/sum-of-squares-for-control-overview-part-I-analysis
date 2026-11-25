% --- предполагается, что Vn уже получена из программы 3.1 ---

gamma_low  = 0;
gamma_high = 20;
tol        = 1e-3;

gamma_best = gamma_low;
%реализация бисекции по gamma
while (gamma_high - gamma_low) > tol
    gamma = 0.5*(gamma_low + gamma_high);
    if is_feasible_gamma(gamma, Vn, r)
        gamma_low  = gamma;
        gamma_best = gamma;
    else
        gamma_high = gamma;
    end
    fprintf('gamma_low = %.6f, gamma_high = %.6f\n', gamma_low, gamma_high);
end

gamma = gamma_best