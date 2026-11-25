clear; clc;

r_low  = 0.1;      % небольшое значение r для бисекции (r=0 означало бы, что область G пустая)
r_high = 10;
tol    = 0.001;

V_best = [];
r_best = r_low;
% реализация бисекции
while (r_high - r_low) > tol
    r = (r_low + r_high)/2;
    [ok, V_cur] = is_feasible_r(r);
    if ok
        r_low  = r;
        r_best = r;
        V_best = V_cur;
    else
        r_high = r;
    end
    fprintf('r_low = %.6f, r_high = %.6f, ok = %d\n', r_low, r_high, ok);
end

r = r_best
r^2
Vn = V_best   % полученная функция Ляпунова