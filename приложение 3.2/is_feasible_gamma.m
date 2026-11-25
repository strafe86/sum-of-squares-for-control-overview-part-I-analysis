function ok = is_feasible_gamma(gamma, Vn, r)
    pvar x y

    Vg = gamma - Vn;              % линия уровня gamma функции Ляпунова Vn
    g  = r^2 - (x^2 + y^2);       % область G

    prog = sosprogram([x y]); % создаем программу
    Z2 = monomials([x y], 0:2); % создаем мономы со степенями от нулевой до второй

    [prog, s] = sossosvar(prog, Z2); %добавляем в программу SOS-переменную, максимальная степень которой равна 2*2=4
    prog = sosineq(prog, g - s*Vg); %добавляем в программу неравенство из программы 3.2

    options.solver = 'sedumi';

    ok = false;
    try
        prog = sossolve(prog, options); %решаем программу
    catch
        return;
    end

    info = prog.solinfo.info;
    if info.numerr ~= 0 || info.pinf ~= 0 || info.dinf ~= 0 || info.feasratio < 0.9 %проверяем решение на корректность
        return;
    end

    ok = true;
end