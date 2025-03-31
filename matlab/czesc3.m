function czesc3()

% Parametry
R1 = 0.1; R2 = 10; C = 0.5; L1 = 3; L2 = 5;
h1 = 1/512; % bardzo krótki krok
h2 = 1/32; % bardzo dlugi krok

% M(u1)
u_values = [20, 50, 100, 150, 200, 250, 280, 300]';
M_values = [0.46, 0.64, 0.78, 0.68, 0.44, 0.23, 0.18, 0.18]';

p = poly_interp(u_values, M_values);
fp = @(u) polyval(p, min(u, 300));

% wymuszenia
e_1 = @(t) ones(length(t), 1);
e_rect = @(t) 120 * (mod(t, 3) < 1.5); 
e_sin = @(t) sin(t);
e_sin1 = @(t) 240 * sin(t); 
e_sin5 = @(t) 210 * sin(2 * pi * 5 * t); 
e_sin50 = @(t) 120 * sin(2 * pi * 50 * t);



function res = deriv_fun(t_i, y_i, e, fun_M)
    res_e = e(t_i);
    M = fun_M(abs(res_e - y_i(3) - y_i(1) * R1));
    D1 = (L1 / M) - (M / L2);
    D2 = (M / L1) - (L2 / M);
    A = [
        -R1 / (M * D1), R2 / (L2 * D1), -1 / (M * D1);
        -R1 / (L1 * D2), R2 / (M * D2), -1 / (L1 * D2);
        1 / C,           0,              0
    ];
    B = [1 / (M * D1); 1 / (L1 * D2); 0];
    res = A * y_i + res_e * B;
end

t1 = 0 : h1 : 30;
t2 = 0 : h2 : 30;
y_0 = [0; 0; 0];

function subintg = simulate_power_nonlin(e, t, h)
    y = solve_ieuler(@(t0,y0) deriv_fun(t0, y0, e, fp), t1, h1, y_0);
    subintg = R1 * y(1,:) .^ 2 + R2 * y(2,:) .^ 2;
    power_rect = rect_intg(subintg, h);
    power_parab = parab_intg(subintg, h);
    fprintf('      met. prostokatów: %.10f\n', power_rect);
    fprintf('      met. parabol:     %.10f\n', power_parab);
end

function subintg = simulate_power_lin(e, t, h)
    M = 0.8;
    D1 = (L1 / M) - (M / L2);
    D2 = (M / L1) - (L2 / M);

    A = [
        -R1 / (M * D1), R2 / (L2 * D1), -1 / (M * D1);
        -R1 / (L1 * D2), R2 / (M * D2), -1 / (L1 * D2);
        1 / C,           0,              0
    ];
    B = [1 / (M * D1); 1 / (L1 * D2); 0];
    y = solve_ieuler(@(t_0,y_0) A * y_0 + e(t_0) * B, t1, h1, y_0);
    subintg = R1 * y(1,:) .^ 2 + R2 * y(2,:) .^ 2;
    power_rect = rect_intg(subintg, h);
    power_parab = parab_intg(subintg, h);
    fprintf('      met. prostokatów: %.10f\n', power_rect);
    fprintf('      met. parabol:     %.10f\n', power_parab);
end

function simulate_obwod_lin(e, title_text)
    fprintf('  liniowa indukcyjnosc wzajemna\n');
    fprintf('    krótki krok\n');
    subintg1 = simulate_power_lin(e, t1, h1);
    fprintf('    dlugi krok\n');
    simulate_power_lin(e, t2, h2);

    figure;
    plot(t1, subintg1); grid on;
    xlabel('Czas [s]'); ylabel('Napiêcie [V]');
    title(sprintf('Lin. - %s', title_text));
end

function simulate_obwod_nonlin(e, title_text)
    fprintf('  nieliniowa indukcyjnosc wzajemna\n');
    fprintf('    krótki krok\n');
    subintg1 = simulate_power_nonlin(e, t1, h1);
    fprintf('    dlugi krok\n');
    simulate_power_nonlin(e, t2, h2);

    figure;
    plot(t1, subintg1); grid on;
    xlabel('Czas [s]'); ylabel('Napiêcie [V]');
    title(sprintf('Nielin. - %s', title_text));
end

function simulate_obwod(e, title_text)
    fprintf('%s\n', title_text);
    simulate_obwod_lin(e, title_text);
    simulate_obwod_nonlin(e, title_text);
end

simulate_obwod(e_rect, 'Prostokatne wymuszenie');
simulate_obwod(e_1, 'Stale wymuszenie');
simulate_obwod(e_sin, 'Sinusoidalne wymuszenie (sin(t))');
simulate_obwod(e_sin1, 'Sinusoidalne wymuszenie (240*sin(t))');
simulate_obwod(e_sin5, 'Sinusoidalne wymuszenie (210*sin(2*pi*5*t))');
simulate_obwod(e_sin50, 'Sinusoidalne wymuszenie (120*sin(2*pi*50*t))');

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

function p = poly_interp(x, y)
    n = length(x) - 1;

    for st = 1 : n
        k = (n + 1) : -1 : (st + 1);
        y(k) = (y(k) - y(k - 1)) ./ (x(k) - x(k - st));
    end

    p = y(n + 1);

    for k = n : -1 : 1
        p = [p, y(k)] + [0, -x(k) * p];
    end
end

function res_intg = rect_intg(y, h)
    res_intg = h * sum(y(1 : (length(y) - 1)));
end

function res_intg = parab_intg(y, h)
    n = length(y);
    sum1 = sum(y(2 : 2 : (n - 1)));
    sum2 = sum(y(3 : 2 : (n - 1)));
    res_intg = h / 3 * (y(1) + y(n) + 4 * sum1 + 2 * sum2);
end