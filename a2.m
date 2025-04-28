clear; clc; close all;
format long;

K0 = 11;
K1 = 0.2;
x0 = 0;
x_slut = 0.5;
h = 0.5/64;

u0 = [0.1; tan(46*pi/180)];

[x, y, interpfel] = Rung2(@fprim, x0, x_slut, u0, h, K0, K1);

fprintf('vid x = %.5f är y = %.10f\n', x(end), y(end,1));

y_svar = y(end,1);

[x_half, y_half] = Rung2(@fprim, x0, x_slut, u0, h/2, K0, K1);

y_ny = y_half(end,1);

absfel = abs(y_ny - y_svar);
relfel = absfel / y_svar;

fprintf('\nDet absoluta felet blir %.10e\nDet relativa felet blir %.10e\n', absfel, relfel);
fprintf('\nInterpolationsfelet blir %.10e', interpfel)

plot(x, y(:,1), 'b-', 'LineWidth', 2);
xlabel('x'); ylabel('u(x)');
title('Lösning med RK4 och interpolation');
grid on;