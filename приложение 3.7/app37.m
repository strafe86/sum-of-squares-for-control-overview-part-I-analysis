%% 1. Определим независимые переменные (состояния и параметры)
pvar x1 x2 x3        % состояния
pvar p1 p2           % параметрическая неопределенность
x = [x1; x2; x3];
p = [p1; p2];

%% 2. Матрица системы A(p)
A = [ -p1,      1,      -1;
       2-2*p2,  2,      -1;
       3,       1,      -p1*p2 ];

%% 3. Домен неопределенности: прямоугольник p1 в [1.7, 3.7], p2 в [5, 9]
% Представим его как полиномиальные неравенства: g1(p) <= 0, g2(p) <= 0
p1_min = 1.7; p1_max = 3.7;
p2_min = 5;   p2_max = 9;
g1 = (p1 - p1_min)*(p1 - p1_max);   % <=0 внутри интервала
g2 = (p2 - p2_min)*(p2 - p2_max);

%% 4. Инициализация SOS-программы
prog = sosprogram([x; p]);   % независимые переменные: x и p

%% 5. Зависимые переменные
% 5.1 зависящая от параметра матрица P(p) из функции Ляпунова (симметричная)
degP = 1;   % максимальная степень p в P(p) (0,1,2)
monomP = monomials(p, 0:degP);
[prog, P] = sospolymatrixvar(prog, monomP, [3,3], 'symmetric');

% Создаем SOS переменные, используя sossosvar.
[prog, sigma11] = sossosvar(prog, mpmonomials({x, p}, {0:2, 0: degP }));
[prog, sigma12] = sossosvar(prog, mpmonomials({x, p}, {0:2, 0: degP }));
[prog, sigma21] = sossosvar(prog, mpmonomials({x, p}, {0:2, 0: degP }));
[prog, sigma22] = sossosvar(prog, mpmonomials({x, p}, {0:2, 0: degP }));

%% 6. Константы
epsilon = 0.1;   % малая положительная величина для обеспечения строго положительной определенности

%% 7. Ограничения
% 7.1 на функцию: V(x;p) - epsilon*(x'*x) + sigma11*g1 + sigma12*g2  is SOS
V = x' * P * x;
pos_constr = V - epsilon*(x1^2 + x2^2 + x3^2) + sigma11*g1 + sigma12*g2;
prog = sosineq(prog, pos_constr, 'sparsemultipartite', {x, p});

% 7.2 производная: -dV/dt - epsilon*(x'*x) + sigma21*g1 + sigma22*g2  is SOS
% вычислим dV/dt = (dV/dx)*dx/dt
dVdt = x'*(P*A + A'*P)*x;
der_constr = -dVdt - epsilon*(x1^2 + x2^2 + x3^2) + sigma21*g1 + sigma22*g2;
prog = sosineq(prog, der_constr, 'sparsemultipartite', {x, p});

%% 8. Решаем программу
options.solver = 'sedumi';

prog = sossolve(prog, options);

%% 9. Проверяем статус решения
info = prog.solinfo.info;
if info.numerr ~= 0
    error('Численные ошибки (numerr = %d). Решение ненадeжно.', info.numerr);
elseif info.pinf ~= 0
    error('Прямая задача недопустима (primal infeasible).');
elseif info.dinf ~= 0
    error('Двойственная задача недопустима (dual infeasible).');
elseif info.feasratio<0.9
    error('Решение не найдено.');
else
    disp('SDP успешно решена, допустимое решение найдено.');
end

%% 10. Извлекаем и выводим результаты
P_sol = sosgetsol(prog, P);
disp('Полученная симметричная матрица P(p):');
disp(char(P_sol));

% Получим конкретную матрицу для номинальной точки (p1=2.7, p2=7)
P_nom = subs(P_sol, [p1; p2], [2.7; 7]);
disp('P(p_nom) ='); disp(double(P_nom));

% Проверим собственные числа P(p_nom) – должны быть положительными
eig_P = eig(double(P_nom));
disp('Собственные числа в номинальной точке:'); disp(eig_P);

%% 11. Продемонстрируем робастную устойчивость в нескольких точках
% Выберем несколько пар параметров внутри прямоугольника
p_samples = [2.0, 9;  2.5, 8.0; 3.0, 6.5; 3.7, 9.0];
figure; hold on;
for k = 1:size(p_samples,1)
    p1k = p_samples(k,1); p2k = p_samples(k,2);
    A_num = double(subs(A, [p1; p2], [p1k; p2k]));
    % моделирование
    sys = ss(A_num, [], eye(3), []);
    [y,t] = initial(sys, [1; 0.5; -0.5], 20);
    plot(t, y(:,1), 'LineWidth',1.2);
end

xlabel('Время, с','FontName', 'Times New Roman'); ylabel('x_1'); grid; title('Траектории для различных значений неопределенности','FontName', 'Times New Roman');
legend('p1=2.0,p2=9','p1=2.5,p2=8.0','p1=3.0,p2=6.5','p1=3.7,p2=9');

V_expr = x' * P_sol * x;                 % явное выражение V
disp('Функция Ляпунова V(x,p) в явном виде:');
disp(char(V_expr));

