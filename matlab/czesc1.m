function czesc1

% Parametry
R1 = 0.1; R2 = 10; C = 0.5; L1 = 3; L2 = 5; M = 0.8; 

% D1 i D2
D1 = (L1 / M) - (M / L2);
D2 = (M / L1) - (L2 / M);

% Macierze A i B
A = [
    -R1 / (M * D1), R2 / (L2 * D1), -1 / (M * D1);
    -R1 / (L1 * D2), R2 / (M * D2), -1 / (L1 * D2);
     1 / C,           0,              0
];
B = [1 / (M * D1); 1 / (L1 * D2); 0];

% wymuszenia
e_rect = @(t) 120 * (mod(t, 3) < 1.5); 
e_sin1 = @(t) 240 * sin(t); 
e_sin5 = @(t) 210 * sin(2 * pi * 5 * t); 
e_sin50 = @(t) 120 * sin(2 * pi * 50 * t); 

% symulacje
simulate_obwod(A, B, e_rect, 1/128, 'Prostokatne wymuszenie');
simulate_obwod(A, B, e_sin1, 1/128, 'Sinusoidalne wymuszenie (240*sin(t))');
simulate_obwod(A, B, e_sin5, 1/64, 'Sinusoidalne wymuszenie (5 Hz)');
simulate_obwod(A, B, e_sin50, 1/32, 'Sinusoidalne wymuszenie (50 Hz)');

end



function y = solve_euler(f, t, h, y_0)
    n = length(t);
    y = zeros(length(y_0),n);
    y(:,1) = y_0;

    for k = 1:(n-1)
        y_0 = y_0 + h * f(t(k),y_0);
        y(:,k+1) = y_0;
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

function simulate_obwod(A, B, e, h, title_text)
    f = @(t_i,y_i) A * y_i + e(t_i) * B;
    y_0 = [0; 0; 0];
    t = 0 : h : 30;
    y_euler = solve_euler(f, t, h, y_0);
    y_ieuler = solve_ieuler(f, t, h, y_0);
    res_e = e(t);

    figure;

    subplot(2, 2, 1);
    plot(t, y_euler(1, :), 'b--', t, y_euler(2, :), 'r-');
    hold on;
    title('Euler');
    legend('i1', 'i2');
    xlabel('Czas [s]'); ylabel('Prąd [A]');
    grid on;

    subplot(2, 2, 2);
    plot(t, y_euler(3, :), 'g--', t, res_e, 'k-');
    title('Euler');
    legend('uc', 'e');
    xlabel('Czas [s]'); ylabel('Napięcie [V]');
    grid on;

    subplot(2, 2, 3);
    plot(t, y_ieuler(1, :), 'b--', t, y_ieuler(2, :), 'r-');
    hold on;
    title('Ulepszony Euler');
    legend('i1', 'i2');
    xlabel('Czas [s]'); ylabel('Prąd [A]');
    grid on;

    subplot(2, 2, 4);
    plot(t, y_ieuler(3, :), 'g--', t, res_e, 'k-');
    title('Ulepszony Euler');
    legend('uc', 'e');
    xlabel('Czas [s]'); ylabel('Napięcie [V]');
    grid on;

    sgtitle(title_text); 
end

