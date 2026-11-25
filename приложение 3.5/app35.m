%% Параметры системы
J1 = 2; J2 = 1; J3 = 1;
k1 = 1; k2 = 1; k3 = 1;

% Правая часть замкнутой системы
f = @(w) [ (-k1/J1)*w(1) + ((J2-J3)/J1)*w(2)*w(3);
           (-k2/J2)*w(2) + ((J3-J1)/J2)*w(3)*w(1);
           (-k3/J3)*w(3) + ((J1-J2)/J3)*w(1)*w(2) ];

%% Функция Ляпунова (без членов < 1e-10)
Vn = @(w1,w2,w3) 1.6449*w1.^4 + 0.6785*w1.^2.*w2.^2 + 0.6785*w1.^2.*w3.^2 + 1.0095*w2.^4 + 2.019*w2.^2.*w3.^2 + 1.0095*w3.^4;

%% 1. Векторное поле
lim = 2.0;
n = 8;
w1 = linspace(-lim, lim, n);
w2 = linspace(-lim, lim, n);
w3 = linspace(-lim, lim, n);
[W1, W2, W3] = meshgrid(w1, w2, w3);

% Вычисляем производные элементов вектора состояния системы
U = zeros(size(W1));
V = zeros(size(W2));
S = zeros(size(W3));
for i = 1:numel(W1)
    w = [W1(i); W2(i); W3(i)];
    dw = f(w);
    U(i) = dw(1);
    V(i) = dw(2);
    S(i) = dw(3);
end

% Нормировка и масштаб
speed = sqrt(U.^2 + V.^2 + S.^2);
speed(speed == 0) = 1;
U = U ./ speed;
V = V ./ speed;
S = S ./ speed;
scale = 0.3;   % длина стрелок
U = scale * U;
V = scale * V;
S = scale * S;

%% 2. Изоповерхности (мелкая сетка для гладкости)
nFine = 25;
wf1 = linspace(-lim, lim, nFine);
wf2 = linspace(-lim, lim, nFine);
wf3 = linspace(-lim, lim, nFine);
[WF1, WF2, WF3] = meshgrid(wf1, wf2, wf3);
Vvals = Vn(WF1, WF2, WF3);

% Уровни функции Ляпунова
vLevels = [1.5, 3.0, 6.0];
% Цвета
colors = [0.85 0.85 0.4;  % жёлтый
          0.5 0.8 0.6;    % мятный
          0.4 0.6 0.9];   % голубой

%% Построение графика
figure;
hold on;

% Стрелки
quiver3(W1, W2, W3, U, V, S, ...
        'Color', [0.3 0.3 0.8], 'LineWidth', 0.6, ...
        'AutoScale', 'off', 'HandleVisibility', 'off');

% Изоповерхности
for i = 1:length(vLevels)
    p = patch(isosurface(WF1, WF2, WF3, Vvals, vLevels(i)));
    p.FaceColor = colors(i,:);
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.45;
    p.SpecularStrength = 0.2;
    p.AmbientStrength = 0.4;
    p.DiffuseStrength = 0.8;
    isonormals(WF1, WF2, WF3, Vvals, p);
    set(p, 'HandleVisibility', 'off');
end

% Освещение
camlight('headlight');
lighting gouraud;

xlabel('$\omega_1$', 'Interpreter', 'latex', ...
       'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('$\omega_2$', 'Interpreter', 'latex', ...
       'FontName', 'Times New Roman', 'FontSize', 14);
zlabel('$\omega_3$', 'Interpreter', 'latex', ...
       'FontName', 'Times New Roman', 'FontSize', 14);

% Заголовок
title('Фазовый портрет системы Эйлера с демпфированием', ...
      'FontName', 'Times New Roman', 'FontSize', 10);
grid on;
axis equal;
view(45, 30);   % ракурс

% Легенда
for i = 1:length(vLevels)
    patch(NaN, NaN, NaN, colors(i,:), 'FaceAlpha', 0.45, ...
          'DisplayName', sprintf('$V_n = %.1f$', vLevels(i)));
end

hLeg = legend('show', 'Location', 'best', ...
              'Interpreter', 'latex', ...
              'FontName', 'Times New Roman', 'FontSize', 11);
set(hLeg, 'Box', 'on');

hold off;
