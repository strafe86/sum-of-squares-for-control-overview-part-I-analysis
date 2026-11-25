%% Параметры системы и результаты SOS-программ
mu = 1;                         % параметр осциллятора
r = 1.6791;                 % радиус множества G (из программы 3.1)
gamma = 8.2074;                   % значение уровня Vn (из программы 3.2)

%% Определение функции Ляпунова Vn (из программы 3.1)
Vn = @(x,y) 0.34446*x^4 + 1.9852*x^3*y + 2.0345*x^2*y^2 - 2.186*x*y^3 + 1.194*y^4 ...
  + 1.2143e-11*x^3 + 5.8294e-11*x^2*y - 8.5159e-11*x*y^2 - 1.5434e-11*y^3 ...
  + 4.0263*x^2 - 6.6736*x*y + 4.3431*y^2 + 1.2619e-14*x - 2.3044e-14*y ...
  + 4.3556e-11;

%% Векторное поле - правая часть системы
% dx/dt = -y
% dy/dt = -mu*(1 - x^2)*y + x
f = @(x,y) [-y; -mu*(1 - x.^2).*y + x];

%% Построение фазового портрета
xLim = [-2.5, 2.5];    % диапазон по оси X
yLim = [-2.5, 2.5];    % диапазон по оси Y
nPts = 30;             % число точек по каждому направлению

x = linspace(xLim(1), xLim(2), nPts);
y = linspace(yLim(1), yLim(2), nPts);
[X, Y] = meshgrid(x, y);

% Вычисляем компоненты векторного поля
U = -Y;                              % dx/dt
V = -mu*(1 - X.^2).*Y + X;           % dy/dt

% нормировка векторов
speed = sqrt(U.^2 + V.^2);
speed(speed == 0) = 1;   % избегаем деления на ноль
U = U ./ speed;
V = V ./ speed;

scale = 0.2;
U = scale * U;
V = scale * V;

figure;
hold on;

% 1. Векторное поле
quiver(X, Y, U, V, 'b', 'LineWidth', 1, 'AutoScale', 'off', ...
       'DisplayName', '${f}(x,y)$');

% 2. Множество G – окружность радиуса r
theta = linspace(0, 2*pi, 200);
xG = r * cos(theta);
yG = r * sin(theta);
plot(xG, yG, 'r--', 'LineWidth', 2, ...
     'DisplayName', sprintf('$x^{2}+y^{2} = %.2f$', r^2));

% 3. Линия уровня Vn = gamma
fContour = @(x,y) Vn(x,y) - gamma;
fimplicit(fContour, [xLim, yLim], 'g-', 'LineWidth', 2, ...
          'DisplayName', sprintf('$V_n = \\gamma = %.2f$', gamma));

% Оформление
xlabel('$x$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('$y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Фазовый портрет осциллятора Ван дер Поля в обратном времени', ...
      'FontName', 'Times New Roman', 'FontSize', 11);
grid on;
axis equal;

% Легенда с LaTeX
hLeg = legend('show', 'Location', 'best', 'Interpreter', 'latex');
set(hLeg, 'FontName', 'Times New Roman', 'FontSize', 12);

hold off;
