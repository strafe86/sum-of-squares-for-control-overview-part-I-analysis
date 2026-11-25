function s = fmt_poly(expr, tol, nd)
%FMT_POLY  Печатает полином expr, отбрасывая члены с |коэффициент| < tol
%          и оставляя ровно nd цифр после запятой.
    expr = expand(expr);
    [C, T] = coeffs(expr);
    s = '';
    for k = 1:numel(C)
        c = double(C(k));
        if abs(c) < tol, continue; end          % отсечение слагаемых
        cr = round(c * 10^nd) / 10^nd;          % округление до nd знаков
        if cr == 0, continue; end
        if strcmp(char(T(k)), '1')
            tk = '';
        else
            tk = char(T(k));
        end
        if isempty(s)
            if cr < 0
                s = sprintf('-%.*f', nd, abs(cr));
            else
                s = sprintf('%.*f', nd, cr);
            end
        else
            if cr < 0
                s = [s sprintf(' - %.*f', nd, abs(cr))];   
            else
                s = [s sprintf(' + %.*f', nd, cr)];        
            end
        end
        if ~isempty(tk)
            s = [s '*' tk];
        end
    end
    if isempty(s), s = '0'; end
end