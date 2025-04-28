clear; clc; close all;
format long;

K1 = 0.2;                 
x0 = 0;
x_slut = 0.5;
h = 0.5/64;               
desiredSlope = -0.51;     
target_y_max = 0.255;     
tol_K = 1e-7;             
tol_u0 = 1e-7;            
maxIter = 100;            

u0_deriv_1 = tan(46*pi/180);      
u0_deriv_2 = u0_deriv_1 + 0.05;     

[f1, opt_K0_1, x_sol_1, y_sol_1, y_max_1, trunk_err_1, rel_trunk_err_1, interp_err_1, iterK1] = ...
    solve_inner(u0_deriv_1, h, x0, x_slut, K1, desiredSlope, target_y_max, tol_K, maxIter);
[f2, opt_K0_2, x_sol_2, y_sol_2, y_max_2, trunk_err_2, rel_trunk_err_2, interp_err_2, iterK2] = ...
    solve_inner(u0_deriv_2, h, x0, x_slut, K1, desiredSlope, target_y_max, tol_K, maxIter);

iter_outer = 0;
while abs(f2) > tol_u0 && iter_outer < maxIter
    delta_u0 = - f2*(u0_deriv_2 - u0_deriv_1) / (f2 - f1);
    u0_deriv_next = u0_deriv_2 + delta_u0;
    
    u0_deriv_1 = u0_deriv_2;
    f1 = f2;
    u0_deriv_2 = u0_deriv_next;
    
    [f2, opt_K0_2, x_sol_2, y_sol_2, y_max_2, trunk_err_2, rel_trunk_err_2, interp_err_2, iterK2] = ...
         solve_inner(u0_deriv_2, h, x0, x_slut, K1, desiredSlope, target_y_max, tol_K, maxIter);
    
    iter_outer = iter_outer + 1;
end

fprintf('Optimal initiala derivatan lutningen blir %.10f\n', u0_deriv_2);
fprintf('Optimal parameter K0 blir %.10f\n', opt_K0_2);
fprintf('\nBeräknad maxhöjd är %.10f (target = %.10f)\n', y_max_2, target_y_max);
fprintf('Absolut fel i maxhöjd blir %.10e\n', abs(y_max_2 - target_y_max));
fprintf('\nTrunkeringsfel (absolut) blir %.10e\n', trunk_err_2);
fprintf('Trunkeringsfel (relativt) blir %.10e\n', rel_trunk_err_2);
fprintf('\nInterpolationsfel blir %.10e\n', interp_err_2);

figure;
plot(x_sol_2, y_sol_2(:,1), 'b-', 'LineWidth', 2);
xlabel('x'); ylabel('y(x)');
title('Lösning med optimerade parametrar, d2.m');
grid on;

% --- Funktioner ---

function [opt_K0, x_out, y_out, interp_error, iterK] = solve_K0(u0_deriv, x0, x_slut, h, K1, desiredSlope, tol_K, maxIter)
    % Bestäm K0 med hjälp av K0_ber istället för att använda egen iteration.
    u0 = [0.1; u0_deriv];
    [opt_K0, ~, interp_error] = K0_ber(x0, x_slut, u0, h, K1, desiredSlope, tol_K);
    [x_out, y_out, ~] = Rung2(@fprim, x0, x_slut, u0, h, opt_K0, K1);
    iterK = NaN;
end

function [x_max, y_max] = maxHeight(x, y)
    x_fine = linspace(x(1), x(end), 1000);
    y_fine = zeros(1, length(x_fine));
    for idx = 1:length(x_fine)
        xq = x_fine(idx);
        idx_nodes = choose_four_nodes(x, xq);
        y_fine(idx) = cubicLagrangeScalar(x(idx_nodes), y(idx_nodes,1), xq);
    end
    [y_max, ind] = max(y_fine);
    x_max = x_fine(ind);
end

function [f_val, opt_K0, x_sol, y_sol, y_max, trunk_err, rel_trunk_err, interp_err, iterK] = solve_inner(u0_deriv, h, x0, x_slut, K1, desiredSlope, target_y_max, tol_K, maxIter)
    [opt_K0, x_sol, y_sol, interp_err, iterK] = solve_K0(u0_deriv, x0, x_slut, h, K1, desiredSlope, tol_K, maxIter);
    [x_max, y_max] = maxHeight(x_sol, y_sol);
    
    % Beräkna trunkeringsfel med steglängden h och h/2
    h_half = h/2;
    u0 = [0.1; u0_deriv];
    [x_half, y_half, interp_err_half] = Rung2(@fprim, x0, x_slut, u0, h_half, opt_K0, K1);
    [~, y_max_half] = maxHeight(x_half, y_half);
    
    trunk_err = abs(y_max_half - y_max);
    rel_trunk_err = trunk_err / abs(y_max_half);
    
    f_val = y_max - target_y_max;
    interp_err = max(interp_err, interp_err_half);
end

% Hjälpfunktion: Kubisk Lagrange-interpolering för en enskild skalfunktion.
function val = cubicLagrangeScalar(x_nodes, y_nodes, xq)
    n = length(x_nodes);
    val = 0;
    for i = 1:n
        L = 1;
        for j = 1:n
            if j ~= i
                L = L * (xq - x_nodes(j)) / (x_nodes(i) - x_nodes(j));
            end
        end
        val = val + y_nodes(i) * L;
    end
end

% Hjälpfunktion: Välj fyra noder från sorterat x för interpolation vid xq.
function idx = choose_four_nodes(x, xq)
    n = length(x);
    if n <= 4
        idx = 1:n;
        return;
    end
    i = find(x <= xq, 1, 'last');
    if isempty(i)
        i = 1;
    end
    if i < 2
        idx = 1:4;
    elseif i > n-2
        idx = n-3:n;
    else
        idx = (i-1):(i+2);
    end
end
