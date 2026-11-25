% маятник с демпфированием 
%
%   theta_dot = omega
%   omega_dot = -(g/l) sin(theta) - c*omega
%
% Преобразование:
%   z1 = sin(theta), z2 = cos(theta), zeta2 = z2 - 1 % zeta вводится для того, чтобы точкой равновесия расширенной системы была точка (0,0)
%   theta_dot  = omega
%   omega_dot  = -(g/l) z1 - c*omega
%   z1_dot     = (zeta2 + 1) * omega
%   zeta2_dot  = -z1 * omega
%
% Ограничение: G1 = z1^2 + z2^2 - 1 = 0
% Домен:
%   D1 = {theta^2 + omega^2 < r^2}
%   D2 = {z1^2 + zeta2^2 < rho^2}     (поскольку (z2-1)^2 = zeta2^2)
% с r <= 2*arcsin(rho/2)
%
%   G_D = [r^2 - theta^2 - omega^2 ;
%          rho^2 - z1^2 - zeta2^2]

clear; clc;

%% ---------------- Параметры ----------------
g_val = 9.81; l_val = 1;      % маятник
c_damp = 0.5;                 % демпфирование

r   = 0.2;                    % радиус D1
rho = 0.25;                    % радиус D2
assert(r <= 2*asin(rho/2) + 1e-9, ...
       'Ошибка: требуется r <= 2*arcsin(rho/2)');

eps_val = 1e-4;                
deg_V   = 2;                  % степень функции Ляпунова
deg_lam = 2;                  % степень lambda_1, lambda_2
deg_mu  = 2;                  % степень SOS-множителей mu_ij 
deg_phi = 2;				  % степень phi

%% ---------------- Переменные ----------------
syms theta omega z1 zeta2 real
vars = [theta; omega; z1; zeta2];

%% ---------------- Преобразованная динамика ----------------
z2 = zeta2 + 1;
f1 = [omega; -(g_val/l_val)*z1 - c_damp*omega];
f2 = [z2*omega; -z1*omega];

%% ---------------- Ограничения ----------------
G1  = z1^2 + z2^2 - 1;          % = z1^2 + (zeta2+1)^2 - 1
GD1 = r^2   - theta^2 - omega^2;
GD2 = rho^2 - z1^2 - zeta2^2;

%% ---------------- SOSTOOLS program ----------------
Program = sosprogram(vars);

% создаем переменную V на мономах степени от 1 до deg_V
[Program, V] = sospolyvar(Program, monomials(vars, [1 deg_V]));

% создаем SOS-переменную phi на мономах степени от 0 до deg_phi
[Program, phi] = sospolyvar(Program, monomials(vars, [0 deg_phi]));
% добавлем в программу условие phi - eps*(||x||^2 + ||z-z0||^2) in SOS
Program = sosineq(Program, ...
    phi - eps_val*(theta^2 + omega^2 + z1^2 + zeta2^2));

% создаем переменные lambda_1, lambda_2 на мономах степени от 0 до deg_lam
[Program, lam1] = sospolyvar(Program, monomials(vars, [0 deg_lam]));
[Program, lam2] = sospolyvar(Program, monomials(vars, [0 deg_lam]));

%  создаем SOS-переменные mu_1 = (mu11, mu12), mu_2 = (mu21, mu22) на мономах степени от 0 до deg_mu
[Program, mu11] = sossosvar(Program, monomials(vars, [0 deg_mu/2]));
[Program, mu12] = sossosvar(Program, monomials(vars, [0 deg_mu/2]));
[Program, mu21] = sossosvar(Program, monomials(vars, [0 deg_mu/2]));
[Program, mu22] = sossosvar(Program, monomials(vars, [0 deg_mu/2]));

% добавляем первое ограничение: V - lam1*G1 - mu1'*GD - phi in SOS
Program = sosineq(Program, ...
    V - lam1*G1 - mu11*GD1 - mu12*GD2 - phi);

% добавляем второе ограничение: -dV/dt - lam2*G1 - mu2'*GD in SOS
dV_dt = diff(V, theta)*f1(1) + diff(V, omega)*f1(2) ...
      + diff(V, z1)*f2(1)    + diff(V, zeta2)*f2(2);
Program = sosineq(Program, ...
    -dV_dt - lam2*G1 - mu21*GD1 - mu22*GD2);

%% ---------------- решаем SOS-программу ----------------
Program = sossolve(Program, struct('solver', 'mosek')); %если у вас нет mosek, замените на sedumi. Однако при этом следует помнить, что вы получите результат, который будет отличаться от представленного в статье

info = Program.solinfo.info;
if info.numerr ~= 0 || info.pinf ~= 0 || info.dinf ~= 0 || info.feasratio < 0.8 % проверяем корректность полученного решения
   error('Решение не найдено.');
end

%% ---------------- получаем найденное решение V(theta, omega, z1, zeta2) ----------------
Vsol   = sosgetsol(Program, V);
phisol = sosgetsol(Program, phi);

%z2 в workspace — это zeta2+1, поэтому создаём отдельный символ
z2_sym  = sym('z2');
Vsol_z2 = expand(subs(Vsol, zeta2, z2_sym - 1));

% подстановка z1 = sin(theta), zeta2 = cos(theta) - 1
V0 = subs(Vsol, [z1, zeta2], [sin(theta), cos(theta) - 1]);
V0 = simplify(V0);

%% ---------------- вывод на экран ----------------
tol = 1e-5;
nd  = 4;

str_V    = fmt_poly(Vsol,    tol, nd);
str_V_z2 = fmt_poly(Vsol_z2, tol, nd);

s_aux = sym('s_aux');
c_aux = sym('c_aux');
V0_aux = subs(V0, [sin(theta), cos(theta)], [s_aux, c_aux]);
str_V0 = fmt_poly(V0_aux, tol, nd);
str_V0 = strrep(str_V0, 's_aux', 'sin(theta)');
str_V0 = strrep(str_V0, 'c_aux', 'cos(theta)');

fprintf('\nV(theta, omega, z1, zeta2) =\n%s\n', str_V);
fprintf('\nV(theta, omega, z1, z2) =\n%s\n', str_V_z2);
fprintf('\nV0(theta, omega) = %s\n', str_V0);
