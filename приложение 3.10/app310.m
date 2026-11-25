clear; clc;

%% 1. Определение независимых переменных
pvar x1 x2 p;                 
vars = [x1; x2; p];           

%% 2. Инициализация SOS-программы
prog = sosprogram(vars);

%% 3. Параметры и константы
epsilon = 1e-4;               % малая положительная константа для строгих неравенств

% Мономы для барьерного сертификата (степень 5 по x1,x2, линейно по p)
% Будем искать B(x,p) = B0(x) + p * B1(x)
degB = 5;                     % максимальная степень по (x1,x2)
Z = monomials([x1; x2], 0:degB);   % все мономы от x1,x2 степени 0..5

% Объявление полиномиальных переменных B0 и B1
[prog, B0] = sospolyvar(prog, Z, 'wscoeff');
[prog, B1] = sospolyvar(prog, Z, 'wscoeff');

% Барьерный сертификат
B = B0 + p * B1;

%% 4. Описание множеств
% Небезопасное множество: g_u(x) ≤ 0, где g_u = (x1+1)^2 + (x2+1)^2 - 0.16
g_u = (x1+1)^2 + (x2+1)^2 - 0.16;

% Начальное множество: g_i(x) ≤ 0, где g_i = (x1-1.5)^2 + x2^2 - 0.25
g_i = (x1-1.5)^2 + x2^2 - 0.25;

% Ограничения на параметр p: p-0.9 ≥ 0, 1.1-p ≥ 0
h1 = p - 0.9;
h2 = 1.1 - p;

%% 5. Множители S-процедуры (SOS-переменные)
% Для каждого ограничения введём SOS-множители.

% Степень множителей
degMult = 2;                  % степень множителей по (x1,x2,p)
Zmult = monomials(vars, 0:degMult);

% --- Для условия на небезопасное множество: B - epsilon >= 0 на g_u<=0, h1>=0, h2>=0 ---
[prog, s_u] = sossosvar(prog, Zmult);   % множитель при -g_u (т.к. -g_u >=0)
[prog, t_u] = sossosvar(prog, Zmult);   % множитель при h1
[prog, r_u] = sossosvar(prog, Zmult);   % множитель при h2
% Условие: B - epsilon - s_u*(-g_u) - t_u*h1 - r_u*h2 >= 0
expr_u = B - epsilon - s_u * (-g_u) - t_u * h1 - r_u * h2;
prog = sosineq(prog, expr_u);

% --- Для условия на начальное множество: -B >= 0 на g_i<=0, h1>=0, h2>=0 ---
[prog, s_i] = sossosvar(prog, Zmult);
[prog, t_i] = sossosvar(prog, Zmult);
[prog, r_i] = sossosvar(prog, Zmult);
expr_i = -B - epsilon - s_i * (-g_i) - t_i * h1 - r_i * h2;
prog = sosineq(prog, expr_i);

%% 6. Условие на производную: -d/dt B >= 0 для всех x1,x2 и p \in [0.9,1.1]
% Производная: d/dt B = (dB/dx1)*x2 + (dB/dx2)*(-x1 + (p/3)*x1^3 - x2)
% Вычисляем частные производные
dB_dx1 = diff(B, x1);
dB_dx2 = diff(B, x2);

% Векторное поле
f1 = x2;
f2 = -x1 + (p/3)*x1^3 - x2;

% Производная
Bdot = dB_dx1 * f1 + dB_dx2 * f2;

% Требуем -Bdot >= 0 на h1>=0, h2>=0 (ограничения на p) и для всех x
[prog, lambda1] = sossosvar(prog, Zmult);   % при h1
[prog, lambda2] = sossosvar(prog, Zmult);   % при h2
expr_dot = -Bdot - lambda1 * h1 - lambda2 * h2;
prog = sosineq(prog, expr_dot);

%% 7. Запуск решателя
options.solver = 'sedumi';

prog = sossolve(prog, options);

% проверка того, что задача была решена успешно
info = prog.solinfo.info;
if info.numerr ~= 0
    error('Численные ошибки (numerr = %d). Решение ненадёжно.', info.numerr);
elseif info.pinf ~= 0
    error('Прямая задача недопустима (primal infeasible).');
elseif info.dinf ~= 0
    error('Двойственная задача недопустима (dual infeasible).');
elseif info.feasratio<0.9
    error('Решение не найдено.');
else
    disp('SDP успешно решена, допустимое решение найдено.');
end

%% 8. Получение найденного решения
B0_sol = sosgetsol(prog, B0);
B1_sol = sosgetsol(prog, B1);
B_sol = sosgetsol(prog, B);

fprintf('\n=== Полученный барьерный сертификат ===\n');
fprintf('B0(x1,x2) = \n'); disp(char(B0_sol));
fprintf('B1(x1,x2) = \n'); disp(char(B1_sol));
fprintf('B(x1,x2,p) = B0 + p*B1 = \n'); disp(char(B_sol));


%% --- Times New Roman для всех подписей ---
set(groot, 'defaultAxesFontName',       'Times New Roman');
set(groot, 'defaultTextFontName',       'Times New Roman');
set(groot, 'defaultLegendFontName',     'Times New Roman');
set(groot, 'defaultColorbarFontName',   'Times New Roman');

%% ============================================================
%%  РИСУНОК 1: зоны сертификата + векторное поле (p = 1)
%% ============================================================
coeffs  = B_sol.coefficient;
degmat  = B_sol.degmat;
vnames  = B_sol.varname;
idx_x1  = find(strcmp(vnames, 'x1'));
idx_x2  = find(strcmp(vnames, 'x2'));
idx_p   = find(strcmp(vnames, 'p'));

p_val   = 1;
x1_vec  = linspace(-2.5, 2.5, 200);
x2_vec  = linspace(-2,   2,   200);
[X1, X2] = meshgrid(x1_vec, x2_vec);
B_vals  = zeros(size(X1));
for i = 1:length(coeffs)
    a = degmat(i, idx_x1);
    b = degmat(i, idx_x2);
    c = degmat(i, idx_p);
    B_vals = B_vals + coeffs(i) * (X1.^a) .* (X2.^b) * (p_val^c);
end

fig1 = figure('Name','Барьерный сертификат','Color','w');
ax1  = axes('Parent', fig1); hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');

% --- B = 0 ---
contour(ax1, X1, X2, B_vals, [0 0], 'b--', 'LineWidth', 2, ...
        'DisplayName', '\itB\rm=0');

% --- Начальное и небезопасное множества ---
theta = linspace(0, 2*pi, 200);
fill(ax1, 1.5 + 0.5*cos(theta), 0 + 0.5*sin(theta), [0.5 0.9 0.5], ...
     'FaceAlpha', 0.35, 'EdgeColor', 'k', ...
     'DisplayName', 'Начальное множество');
fill(ax1, -1 + 0.4*cos(theta), -1 + 0.4*sin(theta), [0.95 0.5 0.5], ...
     'FaceAlpha', 0.35, 'EdgeColor', 'k', ...
     'DisplayName', 'Небезопасное множество');

% --- Траектории ---
odefun = @(t, x) [x(2); -x(1) + (1/3)*x(1)^3 - x(2)];
start_points = [1.5, 0; 1.6, 0.2; 1.4, -0.2];
for i = 1:size(start_points,1)
    [~, x] = ode45(odefun, [0 15], start_points(i,:));
    if i == 1
        plot(ax1, x(:,1), x(:,2), 'k-', 'LineWidth', 1, ...
             'DisplayName', 'Траектории');
    else
        plot(ax1, x(:,1), x(:,2), 'k-', 'LineWidth', 1, ...
             'HandleVisibility', 'off');
    end
end

% --- Векторное поле ---
[Xq, Yq] = meshgrid(linspace(-2.5, 2.5, 20), linspace(-2, 2, 20));
U = Yq;
V = -Xq + (1/3)*Xq.^3 - Yq;
L = sqrt(U.^2 + V.^2);  L(L==0) = 1;
U_n = U ./ L;  V_n = V ./ L;
quiver(ax1, Xq, Yq, U_n, V_n, 0.3, 'Color', [0.5 0.5 0.5], ...
       'LineWidth', 0.5, 'AutoScale', 'off', 'HandleVisibility', 'off');

% --- Равновесия ---
plot(ax1, 0, 0, 'go', 'MarkerSize', 9, 'MarkerFaceColor', 'g', ...
     'HandleVisibility', 'off');
plot(ax1, sqrt(3), 0, 'rs', 'MarkerSize', 9, 'MarkerFaceColor', 'r', ...
     'HandleVisibility', 'off');
plot(ax1, -sqrt(3), 0, 'rs', 'MarkerSize', 9, 'MarkerFaceColor', 'r', ...
     'HandleVisibility', 'off');

text(ax1, 0.15, 0.15, 'Устойчивый фокус', ...
     'FontSize', 9, 'Interpreter', 'tex');
text(ax1, sqrt(3)+0.15, 0.15, 'Седло', ...
     'FontSize', 9, 'Interpreter', 'tex');
text(ax1, -sqrt(3)-0.65, 0.15, 'Седло', ...
     'FontSize', 9, 'Interpreter', 'tex');

% --- Оформление ---
xlabel(ax1, '$x_1$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax1, '$x_2$', 'Interpreter', 'latex', 'FontSize', 14);
title (ax1, 'Барьерный сертификат (\itB\rm(\itx\rm)=0) для \itp\rm=1', ...
       'Interpreter', 'tex', 'FontSize', 12);

hLeg1 = legend(ax1, 'show', 'Location', 'best', 'Interpreter', 'tex');
set(hLeg1, 'FontSize', 10, 'Box', 'on');
hold(ax1, 'off');


%% ============================================================
%%  РИСУНОК 2: B(t) и \dot{B}(t) вдоль траекторий (p = 1)
%% ============================================================
% Полином B(x) = B0 + p*B1 при p = 1
pvar x1 x2
if isa(B0_sol, 'dpvar')
    B0_poly = dpvar2poly(B0_sol);
else
    B0_poly = B0_sol;
end
if isa(B1_sol, 'dpvar')
    B1_poly = dpvar2poly(B1_sol);
else
    B1_poly = B1_sol;
end
B_fixed_poly = B0_poly + B1_poly;

dBdx1 = diff(B_fixed_poly, x1);
dBdx2 = diff(B_fixed_poly, x2);
f1 = x2;
f2 = -x1 + (1/3)*x1^3 - x2;
Bdot_poly = dBdx1 * f1 + dBdx2 * f2;

coeffs_B    = B_fixed_poly.coefficient(:);
degmat_B    = B_fixed_poly.degmat;
vnames_B    = B_fixed_poly.varname;
idx_x1_B    = find(strcmp(vnames_B, 'x1'));
idx_x2_B    = find(strcmp(vnames_B, 'x2'));

coeffs_Bdot = Bdot_poly.coefficient(:);
degmat_Bdot = Bdot_poly.degmat;
vnames_Bdot = Bdot_poly.varname;
idx_x1_Bdot = find(strcmp(vnames_Bdot, 'x1'));
idx_x2_Bdot = find(strcmp(vnames_Bdot, 'x2'));

start_points = [1.5, 0; 1.6, 0.2; 1.4, -0.2];
tspan        = [0 15];
colors       = {'r', 'g', 'b'};
labels       = {'Траектория 1', 'Траектория 2', 'Траектория 3'};

fig2 = figure('Name', 'B(t) и dB/dt(t)', 'Color', 'w');
axB  = subplot(2,1,1, 'Parent', fig2); hold(axB,'on'); grid(axB,'on');
axD  = subplot(2,1,2, 'Parent', fig2); hold(axD,'on'); grid(axD,'on');

for i = 1:size(start_points,1)
    [t, x] = ode45(odefun, tspan, start_points(i,:));
    Nt = numel(t);

    Bv = zeros(Nt,1);
    for k = 1:Nt
        s = 0;
        for m = 1:numel(coeffs_B)
            s = s + coeffs_B(m) * x(k,1)^degmat_B(m,idx_x1_B) ...
                            * x(k,2)^degmat_B(m,idx_x2_B);
        end
        Bv(k) = s;
    end

    Dv = zeros(Nt,1);
    for k = 1:Nt
        s = 0;
        for m = 1:numel(coeffs_Bdot)
            s = s + coeffs_Bdot(m) * x(k,1)^degmat_Bdot(m,idx_x1_Bdot) ...
                                  * x(k,2)^degmat_Bdot(m,idx_x2_Bdot);
        end
        Dv(k) = s;
    end

    plot(axB, t, Bv, colors{i}, 'LineWidth', 1.6, ...
         'DisplayName', labels{i});
    plot(axD, t, Dv, colors{i}, 'LineWidth', 1.6, ...
         'DisplayName', labels{i});

    fprintf('Траектория %d: max B = %.3e, max dB/dt = %.3e\n', ...
            i, max(Bv), max(Dv));
end

% --- Верхний график: B(t) ---
xlabel(axB, '$t$',       'Interpreter', 'latex', 'FontSize', 13);
ylabel(axB, '$B(x)$',    'Interpreter', 'latex', 'FontSize', 13);
title (axB, 'Значение барьерного сертификата при \itp\rm=1', ...
       'Interpreter', 'tex', 'FontSize', 12);
hLB = legend(axB, 'show', 'Location', 'best', 'Interpreter', 'tex');
set(hLB, 'FontSize', 10);

% --- Нижний график: \dot{B}(t) ---
xlabel(axD, '$t$',         'Interpreter', 'latex', 'FontSize', 13);
ylabel(axD, '$\dot{B}$',   'Interpreter', 'latex', 'FontSize', 13);
title (axD, 'Производная \itB\rm вдоль траекторий', ...
       'Interpreter', 'tex', 'FontSize', 12);
hLD = legend(axD, 'show', 'Location', 'best', 'Interpreter', 'tex');
set(hLD, 'FontSize', 10);
