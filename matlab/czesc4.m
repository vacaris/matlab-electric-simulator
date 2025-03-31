function czesc4()

% Parametry
R1 = 0.1; R2 = 10; C = 0.5; L1 = 3; L2 = 5; M = 0.8;
h = 1/256;
tol = 1e-6;

D1 = (L1 / M) - (M / L2);
D2 = (M / L1) - (L2 / M);

A = [
    -R1 / (M * D1), R2 / (L2 * D1), -1 / (M * D1);
    -R1 / (L1 * D2), R2 / (M * D2), -1 / (L1 * D2);
    1 / C,           0,              0
];
B = [1 / (M * D1); 1 / (L1 * D2); 0];
t = 0 : h : 30;
y_0 = [0; 0; 0];

function P = fun_F(f)
    y = solve_ieuler(@(t_i,y_i) A * y_i + (100 * sin(2 * pi * f * t_i)) * B, t, h, y_0);
    P = parab_intg(R1 * y(1,:) .^ 2 + R2 * y(2,:) .^ 2, h) - 406;
end

% wyznaczenie wartosci delta f
f = 0.03;
delta_f = 1;
err = 1;
deriv_f = (fun_F(f + delta_f) - fun_F(f)) / delta_f;
while err > 0.01
    delta_f = 0.5 * delta_f;
    deriv_f0 = deriv_f;
    deriv_f = (fun_F(f + delta_f) - fun_F(f)) / delta_f;
    err = abs(deriv_f - deriv_f0) / abs(deriv_f0);
end
fprintf('Znaleziona wartosci delta f: %e\n\n', delta_f);

f = 0.01 : 0.01 : 1;
n = length(f);
P = zeros(1, n);

for k = 1 : n
    P(k) = fun_F(f(k));
end

j = find(diff(sign(P(2 : n))) ~= 0);

for k = j
    a = f(k+1); b = f(k+2);
    fprintf('Przedzial poszukiwan: (%g,%g)\n', a, b);

    % Metoda Bisekcji
    [f_bisect, P_bisect, iter_bisect] = bisect_method(@fun_F, a, b, tol);
    fprintf('  Metoda Bisekcji\n    f = %.6f Hz\n', f_bisect);
    fprintf('    wartosc funkcji F: %e\n', P_bisect);
    fprintf('    liczba iteracji: %d\n', iter_bisect);
    fprintf('    liczba obliczen mocy: %d\n', iter_bisect + 1);

    % Metoda Siecznych
    [f_secant, P_secant, iter_secant] = secant_method(@fun_F, a, a + h/4, tol);
    fprintf('  Metoda Siecznych: f = %.6f Hz\n', f_secant);
    fprintf('    wartosc funkcji F: %e\n', P_secant);
    fprintf('    liczba iteracji: %d\n', iter_secant);
    fprintf('    liczba obliczen mocy: %d\n', iter_secant + 2);

    % Metoda Quasi-Newtona
    [f_newton, P_newton, iter_newton] = quasi_newton_method(@fun_F, a, delta_f, tol);
    fprintf('  Metoda Quasi-Newtona: f = %.6f Hz\n', f_newton);
    fprintf('    wartosc funkcji F: %e\n', P_newton);
    fprintf('    liczba iteracji: %d\n', iter_newton);
    fprintf('    liczba obliczen mocy: %d\n', 2 * iter_newton + 1);
end

end


function y = solve_ieuler(f, t, h, y_0)
    n = length(t);
    y = zeros(length(y_0),n);
    y(:,1) = y_0;
    hh = h / 2;

    for k = 1:(n-1)
        t_0 = t(k);
        y_0 = y_0 + h * f(t_0 + hh, y_0 + hh * f(t_0, y_0));
        y(:,k+1) = y_0;
    end
end

function res_intg = parab_intg(y, h)
    n = length(y);
    sum1 = sum(y(2 : 2 : (n - 1)));
    sum2 = sum(y(3 : 2 : (n - 1)));
    res_intg = h / 3 * (y(1) + y(n) + 4 * sum1 + 2 * sum2);
end

function [f, P, iter] = bisect_method(calculate_F, f_a, f_b, tol)
    iter = 0;
    F_a = calculate_F(f_a);
    while abs(f_b - f_a) > tol
        f_mid = (f_a + f_b) / 2;
        F_mid = calculate_F(f_mid);
        if F_a * F_mid < 0
            f_b = f_mid;
        else
            f_a = f_mid;
        end
        iter = iter + 1;
    end
    f = f_mid;
    P = F_mid;
end

function [f, P, iter] = secant_method(calculate_F, f0, f1, tol)
    iter = 0;
    F0 = calculate_F(f0);
    F1 = calculate_F(f1);
    while abs(f1 - f0) > tol
        f_new = f1 - F1 * (f1 - f0) / (F1 - F0);
        f0 = f1;
        f1 = f_new;
        F0 = F1;
        F1 = calculate_F(f_new);
        iter = iter + 1;
    end
    f = f1;
    P = F1;
end

function [f, P, iter] = quasi_newton_method(calculate_F, f0, delta_f, tol)
    iter = 0;
    while true
        F = calculate_F(f0);
        F_derivative = (calculate_F(f0 + delta_f) - F) / delta_f;
        f_new = f0 - F / F_derivative;
        if abs(f_new - f0) < tol
            break;
        end
        f0 = f_new;
        iter = iter + 1;
    end
    f = f0;
    P = calculate_F(f);
end