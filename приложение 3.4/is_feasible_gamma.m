% is_feasible_gamma.m 
function [ok, Vn, beta_out] = is_feasible_gamma(gamma)
    pvar w1 w2 w3 %создаем переменные – элементы вектора состояния системы

    J1 = 2; J2 = 1; J3 = 1; %задаем параметры системы
    k1 = 1; k2 = 1; k3 = 1; % задаем параметры законов управления
    epsilon = 0.0001;

    u1 = -k1*w1; %задаем законы управления
    u2 = -k2*w2; %задаем законы управления
    u3 = -k3*w3; %задаем законы управления

    f = [(J2-J3)/J1*w2*w3 + u1;
         (J3-J1)/J2*w3*w1 + u2;
         (J1-J2)/J3*w1*w2 + u3]; % задаем динамику системы

    prog = sosprogram([w1 w2 w3]); % создаем программу

    [prog, beta] = sospolyvar(prog, 1); % создаем скалярную переменную beta

    Z = monomials([w1 w2 w3], 1:2); %создаем мономы из w1, w2, w3 со степенями 1 и 2

    [prog, V] = sossosvar(prog, Z); % добавляем в программу SOS-переменную V, которая будет иметь максимальную степень 2*2=4
    prog = sosineq(prog, V - epsilon*(w1^4 + w2^4 + w3^4)); % условие, гарантирующее радиальную неограниченность функции Ляпунова
    prog = sosineq(prog, -V + beta*(w1^4 + w2^4 + w3^4)); % верхняя граница для V 
    prog = soseq(prog, subs(V, [w1; w2; w3], [0; 0; 0])); % условие-равенство V(0,0,0)=0

    nablaV = [diff(V,w1); diff(V,w2); diff(V,w3)]; %находим градиент V(w1,w2,w3)

    prog = sosineq(prog, -nablaV'*f - gamma*V);  %добавляем в программу условие-неравенство для производной функции Ляпунова

    options.solver = 'sedumi';

    ok = false;
    Vn = []; beta_out = [];

    try
        prog = sossolve(prog, options); %решаем созданную программу
    catch
        return;
    end

    info = prog.solinfo.info;
    if info.numerr ~= 0 || info.pinf ~= 0 || info.dinf ~= 0 || info.feasratio < 0.9 % проверяем корректность решения
        return;
    end

    ok       = true;
    Vn       = sosgetsol(prog, V);
    beta_out = sosgetsol(prog, beta);
end
