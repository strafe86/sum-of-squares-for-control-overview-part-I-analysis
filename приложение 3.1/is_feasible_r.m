function [ok, Vn] = is_feasible_r(r)
    pvar x y
    rho = 1; % параметр осциллятора
    mu = 0.0001;

    f = [-y; -rho*(1-x^2)*y + x];       % динамика системы
    g = r^2 - (x^2 + y^2);             % область G

    prog = sosprogram([x y]); %создаем SOS программу
    Z2 = monomials([x y],0:2); %создаем мономы от нулевой до второй степени

    [prog, V] = sossosvar(prog, Z2); % добавляем в программу новую SOS переменную V, причем она будет иметь степень Z2*2, то есть 4
    prog = soseq(prog, subs(V, [x;y], [0;0])); %создаем условие-равенство V(0,0)=0

    [prog, s] = sossosvar(prog, Z2); % добавляем в программу новую SOS переменную s
    prog = sosineq(prog, V - s*g - mu*(x^2 + y^2)); %гарантируем радиальную неограниченность V на G - первое условие из теоремы 3.4

    nablaV = [diff(V,x); diff(V,y)]; % находим градиент V(x,y)

    [prog, t] = sossosvar(prog, Z2); % добавляем в программу SOS-переменную t, ее максимальный порядок – 2*2=4
    prog = sosineq(prog, -nablaV'*f - t*g - mu*(x^2 + y^2)); % второе условие из теоремы 3.4

    options.solver = 'sedumi';

    ok = false;
    Vn = [];

    try
        prog = sossolve(prog, options); %решаем созданную программу
    catch
        return;
    end

    info = prog.solinfo.info;
    if info.numerr ~= 0 || info.pinf ~= 0 || info.dinf ~= 0 || info.feasratio < 0.9 % проверяем, что получено корректное решение
        return;
    end

    ok = true;
    Vn = sosgetsol(prog, V); % получаем найденную функцию Ляпунова
end