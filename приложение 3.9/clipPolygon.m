%% отсечение выпуклого многоугольника полуплоскостью 
function Pout = clipPolygon(P, c, phi_k)
% Оставляет часть многоугольника P, где c'*x <= phi_k

    if isempty(P)
        Pout = P;
        return;
    end

    Pout = zeros(0,2);
    n = size(P,1);
    tol = 1e-10;

    for i = 1:n
        A = P(i,:);
        B = P(mod(i,n)+1,:);

        fa = c(1)*A(1) + c(2)*A(2) - phi_k;
        fb = c(1)*B(1) + c(2)*B(2) - phi_k;

        inA = fa <= tol;
        inB = fb <= tol;

        if inA
            Pout(end+1,:) = A;
        end

        if inA ~= inB
            t = fa/(fa - fb);
            I = A + t*(B - A);
            Pout(end+1,:) = I;
        end
    end

    % Удаление совпадающих соседних вершин
    if size(Pout,1) > 1
        keep = [true; sum(abs(diff(Pout,1,1)),2) > 1e-10];
        Pout = Pout(keep,:);
    end
end
