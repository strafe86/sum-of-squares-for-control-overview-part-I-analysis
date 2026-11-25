%% Параметрическая матрица A(p1,p2) и расчёт устойчивости
p1_range = linspace(0, 6, 80);
p2_range = linspace(3, 12, 80);
[P1, P2] = meshgrid(p1_range, p2_range);

maxRe = zeros(size(P1));
for i = 1:numel(P1)
    p1 = P1(i);
    p2 = P2(i);
    A_num = [-p1,      1,      -1;
             2-2*p2,  2,      -1;
             3,       1,      -p1*p2];
    maxRe(i) = max(real(eig(A_num)));
end

% Определяем устойчивость
threshold = 0;
stable = maxRe < threshold;

% Прямоугольник области, полученной с помощью SOS (p1 in [1.7, 3.7], p2 in [5, 9])
rect_x = [1.7, 3.7, 3.7, 1.7, 1.7];
rect_y = [5, 5, 9, 9, 5];

%% 1. Двумерное изображение области устойчивости
figure('Name', 'Область устойчивости');
set(gcf, 'DefaultAxesFontName', 'Times New Roman');
set(gcf, 'DefaultTextFontName', 'Times New Roman');

h_img = imagesc(p1_range, p2_range, stable);
set(h_img, 'HandleVisibility', 'off');
axis xy;
xlabel('$p_1$', 'Interpreter', 'latex', ...
       'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('$p_2$', 'Interpreter', 'latex', ...
       'FontName', 'Times New Roman', 'FontSize', 14);
title('Область устойчивости системы', 'FontName', 'Times New Roman', 'FontSize', 14);
colormap([1 0 0; 0 0.7 0]);   % красный / зелёный
hcb = colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'неуст.', 'уст.'});
set(hcb, 'FontName', 'Times New Roman', 'FontSize', 10);
grid on;
hold on;

% Прямоугольник, полученный с помощью SOS (показан в легенде)
plot(rect_x, rect_y, 'k--', 'LineWidth', 2, 'DisplayName', 'SOS-результат');
legend('Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 10);
hold off;

%% 2. график maxRe с нулевой линией
figure('Name', 'Область устойчивости');
set(gcf, 'DefaultAxesFontName', 'Times New Roman');
set(gcf, 'DefaultTextFontName', 'Times New Roman');

% Заливка контуров – нет в легенде
contourf(P1, P2, maxRe, 20, 'LineStyle', 'none', 'HandleVisibility', 'off');
xlabel('$p_1$', 'Interpreter', 'latex', ...
       'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('$p_2$', 'Interpreter', 'latex', ...
       'FontName', 'Times New Roman', 'FontSize', 14);
title('max Re(A)', 'FontName', 'Times New Roman', 'FontSize', 14);
hcb = colorbar;
set(hcb, 'FontName', 'Times New Roman', 'FontSize', 10);
hold on;
% Граница устойчивости (maxRe = 0) – чёрная линия
contour(P1, P2, maxRe, [0 0], 'k-', 'LineWidth', 2, 'DisplayName', 'граница устойчивости');
% Прямоугольник SOS-результата – красный пунктир
plot(rect_x, rect_y, 'r--', 'LineWidth', 2, 'DisplayName', 'SOS-результат');

legend('Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 10);
grid on;
hold off;
