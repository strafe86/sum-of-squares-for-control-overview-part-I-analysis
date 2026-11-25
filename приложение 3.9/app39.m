%% Единый скрипт: множество достижимости + ЛМН-эллипсоид + SOS аппроксимация множества достижимости
clear; clc; close all;

% --- Система ---
A = [0 1; -1 -2];
D = [0; 1];
n = 2; m = 1;
alpha = 0.57;

%% === ЧАСТЬ 1: ЛМН-эллипсоид для ||w||_2 <= 1 ===
Q = sdpvar(n, n, 'symmetric');
M_lmi = [Q*A' + A*Q + alpha*Q, D;
D', -alpha*eye(m)];
C_lmi = [Q >= 1e-3*eye(n), M_lmi <= -1e-3*eye(n+m)];
opts_y = sdpsettings('solver','sedumi','verbose',0);
sol_lmi = optimize(C_lmi, trace(Q), opts_y);

if sol_lmi.problem ~= 0
	error('ЛМН не решены: %s', sol_lmi.info);
end

Q_val = (value(Q) + value(Q)')/2;
[Vq_lmi, Dq_lmi] = eig(Q_val);
axes_lmi = sqrt(diag(Dq_lmi));

th = linspace(0, 2*pi, 500);
pts_lmi = Vq_lmi * diag(axes_lmi) * [cos(th); sin(th)];
fprintf('ЛМН: полуоси = [%.4f, %.4f]\n', axes_lmi(1), axes_lmi(2));

%% === ЧАСТЬ 2: SOS-множество {V<=1}, вложенное в ЛМН-эллипсоид ===
P_lmi = inv(Q_val);

syms x1s x2s w1s y1s y2s y3s y4s real
vars     = [x1s; x2s];
vars_aug = [x1s; x2s; w1s];

options_sos = struct();
options_sos.solver  = 'sedumi';
options_sos.verbose = 0;

eps_val = 1e-3;

prog = sosprogram([vars_aug; y1s; y2s; y3s; y4s]);

% V — полином степени 1..6
[prog, V] = sospolyvar(prog, monomials(vars, 1:6));

% V2 = x'*P2*x, P2 = [p11 p12; p12 p22].
% sospolyvar(prog, sym(1)) создаёт константу
one_mon = sym(1);
[prog, p11] = sospolyvar(prog, one_mon);
[prog, p12] = sospolyvar(prog, one_mon);
[prog, p22] = sospolyvar(prog, one_mon);

% V2 — в явном виде
V2 = p11*x1s^2 + 2*p12*x1s*x2s + p22*x2s^2;

% Xi — вспомогательная матрица для [Xi I; I P2] >= 0
[prog, xi11] = sospolyvar(prog, one_mon);
[prog, xi12] = sospolyvar(prog, one_mon);
[prog, xi22] = sospolyvar(prog, one_mon);

% Условие 1
prog = sosineq(prog, V  - eps_val*(x1s^2 + x2s^2));
prog = sosineq(prog, V2 - eps_val*(x1s^2 + x2s^2));

% mu
[prog, mu] = sossosvar(prog, monomials(vars_aug, 2));

% Условие 2: инвариантность {V<=1}
Vdot = jacobian(V, vars) * (A*vars + D*w1s);
prog = sosineq(prog, -Vdot + alpha*(1 - V) ...
    - mu*(1 - w1s^2) - eps_val*(x1s^2 + x2s^2));

% Условие 3: вложенность
prog = sosineq(prog, V - V2);
prog = sosineq(prog, V2 - vars'*P_lmi*vars);

% LMI [Xi I; I P2] >= 0  <=>  y'*[Xi I; I P2]*y >= 0 для всех y \in R^4
LMI_quad = xi11*y1s^2 + 2*xi12*y1s*y2s + xi22*y2s^2 ...
         + 2*(y1s*y3s + y2s*y4s) ...
         + p11*y3s^2 + 2*p12*y3s*y4s + p22*y4s^2;
prog = sosineq(prog, LMI_quad);

% Целевая функция: tr(Xi) = xi11 + xi22 -> min
prog = sossetobj(prog, xi11 + xi22);

prog = sossolve(prog, options_sos);

% --- Проверка решения ---
info = prog.solinfo.info;
ok = isfield(info,'numerr') && info.numerr==0 && ...
     isfield(info,'pinf')   && info.pinf==0 && ...
     isfield(info,'feasratio') && info.feasratio>=0.9 && ...
     isfield(info,'dinf')   && info.dinf==0;
if ~ok
    error('SOS не решена: %s', info.info);
end
fprintf('SOS-программа решена\n');

% --- Получение решения ---
V_sol = sosgetsol(prog, V);
V2_sol = sosgetsol(prog, V2);

disp('V(x) = '); disp(V_sol);

%% === ЧАСТЬ 3: достижимое множество для ||w||_infty <= 1 ===
nTheta = 400;
theta_p = linspace(0, 2*pi, nTheta+1);
theta_p(end) = [];
T = 50;
phi = zeros(size(theta_p));

for k = 1:numel(theta_p)
	c = [cos(theta_p(k)); sin(theta_p(k))];
	f = @(t) abs(D' * expm(A' * t) * c);
	phi(k) = integral(f, 0, T, 'ArrayValued', true, 'AbsTol',1e-10, 'RelTol',1e-8);
end

Mbig = 20;
Pp = [-Mbig -Mbig; Mbig -Mbig; Mbig Mbig; -Mbig Mbig];

for k = 1:numel(theta_p)
	c = [cos(theta_p(k)); sin(theta_p(k))];
	Pp = clipPolygon(Pp, c, phi(k));
	if isempty(Pp), error('Пустой многоугольник'); end
end

fprintf('Полигон R: %d вершин\n', size(Pp,1));

%% === ЧАСТЬ 4: вывод всех множеств на одном графике ===
V_fun  = matlabFunction(V_sol,  'Vars', {x1s, x2s});
V_fun2 = matlabFunction(V2_sol, 'Vars', {x1s, x2s});

figure('Color','w'); hold on; grid on; axis equal;

h_poly = fill(Pp(:,1), Pp(:,2), [1.00 0.90 0.70], ...
    'EdgeColor',[0.80 0.40 0.00], 'LineWidth',2, 'FaceAlpha',0.35);

h_lmi = fill(pts_lmi(1,:), pts_lmi(2,:), [0.85 0.90 1.00], ...
    'EdgeColor',[0.20 0.30 0.80], 'LineWidth',2, 'FaceAlpha',0.30);

h_sos = fimplicit(@(a,b) V_fun(a,b) - 1, [-2 2 -2 2], ...
    'Color',[0.10 0.55 0.10], 'LineWidth', 2.5);

h_sos1 = fimplicit(@(a,b) V_fun2(a,b) - 1, [-2 2 -2 2], ...
    'Color',[0 0 0], 'LineWidth', 2.5);

plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 1.5);

% --- Оси ---
xlabel('$x_1$', 'Interpreter','latex', ...
       'FontName','Times New Roman', 'FontSize', 13);
ylabel('$x_2$', 'Interpreter','latex', ...
       'FontName','Times New Roman', 'FontSize', 13);

% --- Заголовок ---
title('Множество достижимости, ЛМН-эллипсоид и SOS-множество \itV\rm', ...
      'Interpreter','tex', 'FontName','Times New Roman', 'FontSize', 11);

% --- Легенда ---
hLeg = legend([h_poly, h_lmi, h_sos, h_sos1], ...
    { 'Множество достижимости', ...
      'ЛМН-эллипсоид', ...
      'SOS: \{\itV\rm(x) \leq 1\}', ...
      '\{\itV\rm_2(x) \leq 1\}' }, ...
    'Interpreter', 'tex', ...
    'FontName','Times New Roman', 'FontSize', 10, ...
    'Location', 'best');
set(hLeg, 'Box', 'on');