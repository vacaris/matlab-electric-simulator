function czesc2()

% Parametry
R1 = 0.1; R2 = 10; C = 0.5; L1 = 3; L2 = 5;
h = 1/128;

% M(u1)
u_values = [20, 50, 100, 150, 200, 250, 280, 300]';
M_values = [0.46, 0.64, 0.78, 0.68, 0.44, 0.23, 0.18, 0.18]';

% wielomiany interpolacyjne / aproksymacyjne
p3 = poly_approx(u_values, M_values, 3);
p5 = poly_approx(u_values, M_values, 5);
p = poly_interp(u_values, M_values);
p_deriv = polyder(p);
a = [polyval(p_deriv, u_values(1)), polyval(p_deriv, u_values(end))];
s = cubic_spline(u_values, M_values, a);

% funkcje aproksymujące / interpolujące dane
fp3 = @(u) polyval(p3, min(u, 300));
fp5 = @(u) polyval(p5, min(u, 300));
fp = @(u) polyval(p, min(u, 300));
fs = @(u) spline_val(s, u_values, min(u, 300));

% wykresy wszystkich funkcji interpolacyjnych i aproksymacyjnych
u = 0 : 0.01 : 300;
figure;
subplot(2, 2, 1);
plot(u, fp3(u), u_values, M_values, '.'); grid on;
xlabel('u [V]'); ylabel('M [H]');
title('Aproksymacja st. 3');
subplot(2, 2, 2);
plot(u, fp5(u), u_values, M_values, '.'); grid on;
xlabel('u [V]'); ylabel('M [H]');
title('Aproksymacja st. 5');
subplot(2, 2, 3);
plot(u, fp(u), u_values, M_values, '.'); grid on;
xlabel('u [V]'); ylabel('M [H]');
title('Interpolacja wielomianowa');
subplot(2, 2, 4);
plot(u, fs(u), u_values, M_values, '.'); grid on;
xlabel('u [V]'); ylabel('M [H]');
title('Interpolacja funkcjami sklejanymi st. 3');
sgtitle('Wykresy interpolacji / aproksymacji');

% wymuszenia
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

t = 0 : h : 30;
y_0 = [0; 0; 0];

function simulate_obwod(e, title_text)
    y1 = solve_ieuler(@(t0,y0) deriv_fun(t0, y0, e, fp3), t, h, y_0);
    y2 = solve_ieuler(@(t0,y0) deriv_fun(t0, y0, e, fp5), t, h, y_0);
    y3 = solve_ieuler(@(t0,y0) deriv_fun(t0, y0, e, fp), t, h, y_0);
    y4 = solve_ieuler(@(t0,y0) deriv_fun(t0, y0, e, fs), t, h, y_0);

    % wykresy natężenia i1
    figure;
    subplot(2, 2, 1);
    plot(t, y1(1,:)); grid on;
    xlabel('t [s]'); ylabel('i1 [A]');
    title('Aproksymacja st. 3');
    subplot(2, 2, 2);
    plot(t, y2(1,:)); grid on;
    xlabel('t [s]'); ylabel('i1 [A]');
    title('Aproksymacja st. 5');
    subplot(2, 2, 3);
    plot(t, y3(1,:)); grid on;
    xlabel('t [s]'); ylabel('i1 [A]');
    title('Interpolacja wielomianowa');
    subplot(2, 2, 4);
    plot(t, y4(1,:)); grid on;
    xlabel('t [s]'); ylabel('i1 [A]');
    title('Interpolacja funkcjami sklejanymi st. 3');
    sgtitle(title_text);
    
    % wykresy natężenia i2
    figure;
    subplot(2, 2, 1);
    plot(t, y1(2,:)); grid on;
    xlabel('t [s]'); ylabel('i2 [A]');
    title('Aproksymacja st. 3');
    subplot(2, 2, 2);
    plot(t, y2(2,:)); grid on;
    xlabel('t [s]'); ylabel('i2 [A]');
    title('Aproksymacja st. 5');
    subplot(2, 2, 3);
    plot(t, y3(2,:)); grid on;
    xlabel('t [s]'); ylabel('i2 [A]');
    title('Interpolacja wielomianowa');
    subplot(2, 2, 4);
    plot(t, y4(2,:)); grid on;
    xlabel('t [s]'); ylabel('i2 [A]');
    title('Interpolacja funkcjami sklejanymi st. 3');
    sgtitle(title_text);

    % wykresy napięcia uR2
    figure;
    subplot(2, 2, 1);
    plot(t, R2 * y1(2,:)); grid on;
    xlabel('t [s]'); ylabel('uR2 [V]');
    title('Aproksymacja st. 3');
    subplot(2, 2, 2);
    plot(t, R2 * y2(2,:)); grid on;
    xlabel('t [s]'); ylabel('uR2 [V]');
    title('Aproksymacja st. 5');
    subplot(2, 2, 3);
    plot(t, R2 * y3(2,:)); grid on;
    xlabel('t [s]'); ylabel('uR2 [V]');
    title('Interpolacja wielomianowa');
    subplot(2, 2, 4);
    plot(t, R2 * y4(2,:)); grid on;
    xlabel('t [s]'); ylabel('uR2 [V]');
    title('Interpolacja funkcjami sklejanymi st. 3');
    sgtitle(title_text);

    % wykresy napięcia uC
    figure;
    subplot(2, 2, 1);
    plot(t, y1(3,:)); grid on;
    xlabel('t [s]'); ylabel('uC [V]');
    title('Aproksymacja st. 3');
    subplot(2, 2, 2);
    plot(t, y2(3,:)); grid on;
    xlabel('t [s]'); ylabel('uC [V]');
    title('Aproksymacja st. 5');
    subplot(2, 2, 3);
    plot(t, y3(3,:)); grid on;
    xlabel('t [s]'); ylabel('uC [V]');
    title('Interpolacja wielomianowa');
    subplot(2, 2, 4);
    plot(t, y4(3,:)); grid on;
    xlabel('t [s]'); ylabel('uC [V]');
    title('Interpolacja funkcjami sklejanymi st. 3');
    sgtitle(title_text);
end

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

function p = poly_approx(x, y, n)
    [Q,R] = qr(bsxfun(@power, x, n : -1 : 0));
    n = n + 1;
    b = Q' * y;
    p = R(1:n, 1:n) \ b(1:n);
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

function s = cubic_spline(x, y, a)    
    n = length(x) - 1;
    h = diff(x);
    k = 2 : n;
    d = 6 * [(y(2) - y(1)) / h(1) - a(1);...
             (y(k+1) - y(k)) ./ h(k) - (y(k) - y(k-1)) ./ h(k-1);...
             (a(2) - (y(n+1) - y(n)) / h(n))];

    A = zeros(n+1, n+1);
    A(1,1) = 2 * h(1);
    A(1,2) = h(1);
    A(n+1,n) = h(n);
    A(n+1,n+1) = 2 * h(n);

    for k = 2 : n
        A(k,k-1) = h(k-1);
        A(k,k) = 2 * (h(k-1) + h(k));
        A(k,k+1) = h(k);
    end

    M = A \ d;
    s = zeros(n, 4);

    for k = 1 : n
        s(k,:) = [y(k),...
                  (y(k+1) - y(k)) / h(k) - h(k) * (2*M(k) + M(k+1)) / 6,...
                  M(k) / 2,...
                  (M(k+1) - M(k)) / (6 * h(k))];
    end
end

function yy = spline_val(s, x, xx)
    for i = 1 : (length(x) - 1)
        idx = xx >= x(i) & xx <= x(i+1);
        dx = xx(idx) - x(i);
        yy(idx) = s(i,1) + s(i,2) * dx + s(i,3) * dx.^2 + s(i,4) * dx.^3;
    end

    idx = xx < x(1);
    dx = xx(idx) - x(1);
    yy(idx) = s(1,1) + s(1,2) * dx + s(1,3) * dx.^2 + s(1,4) * dx.^3;
end