clear; clc; close all;
format long;

% Parametrar
K1 = 0.2;               
x0 = 0;                 
x_slut = 0.5;           
h = 0.5/64;             
lutning = -0.51;        
tolerans = 1e-7;        

u0 = [0.1; tan(46*pi/180)];

% Använd K0_ber för att bestämma optimalt K0
[K0_nu, skillnad_nu, interpfel] = K0_ber(x0, x_slut, u0, h, K1, lutning, tolerans);

relfel_K0 = abs(skillnad_nu) / abs(K0_nu);

fprintf('Då lutningen i slutet är -0.51 blir K0 = %.10f \n', K0_nu)
fprintf('\nSkillnad i det beräknade lutningen och -0.51 blir %.10e\n', abs(skillnad_nu))
fprintf('\nDet absoluta felet för detta K0 blir %.10e\n', abs(skillnad_nu))
fprintf('Det relativa felet för detta K0 blir %.10e\n', relfel_K0)
fprintf('\nInterpolationsfelet blir %.10e\n', interpfel)

[x, y] = Rung2(@fprim, x0, x_slut, u0, h, K0_nu, K1);

figure;
plot(x, y(:,1))
grid on;
xlabel('x'), ylabel('y(x)');
title('Lösning med korrigerat K0');
