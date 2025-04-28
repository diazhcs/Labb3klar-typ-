clear; clc; close all;
format long;

% Parametrar
K1 = 0.2;
x_start = 0;
x_slut = 0.5;
h = 0.5/64;
lutning = -0.51;
tolerans = 1e-7;

u0 = [0.1; tan(46*pi/180)];

% Bestäm K0 med hjälp av K0_ber
[K0_nuvarande, diff_val, interp_fel_h] = K0_ber(x_start, x_slut, u0, h, K1, lutning, tolerans);

% Beräkna lösningen med Rung1 med det optimerade K0
[x, y, interp_fel_h] = Rung1(@fprim, x_start, x_slut, u0, h, K0_nuvarande, K1);

% Hitta den lokala maximumhöjden
y_prim_h = y(:,2);
kors_idx = find(diff(sign(y_prim_h)) < 0, 1);
if isempty(kors_idx)
    error('Inget maximum hittades.');
end
i = kors_idx;
x_min = x(i); 
x_max = x(i+1);

x_fine = linspace(x_min, x_max, 50);
y_fine = zeros(1, length(x_fine));
for idx = 1:length(x_fine)
    xq = x_fine(idx);
    nodes_idx = choose_four_nodes(x, xq);
    y_fine(idx) = cubicLagrangeScalar(x(nodes_idx), y(nodes_idx,1), xq);
end
[y_max_h, idx_max] = max(y_fine);
x_max_interp = x_fine(idx_max);

% Lösning med halverad steglängd för trunkeringsfel
h_halv = h/2;
[x_halv, y_halv, interp_fel_halv] = Rung1(@fprim, x_start, x_slut, u0, h_halv, K0_nuvarande, K1);
y_prim_halv = y_halv(:,2);
kors_idx_halv = find(diff(sign(y_prim_halv)) < 0, 1);
if isempty(kors_idx_halv)
    error('Inget maximum hittades med h/2.');
end
i_halv = kors_idx_halv;
x_min_h = x_halv(i_halv); 
x_max_h = x_halv(i_halv+1);

x_fine_halv = linspace(x_min_h, x_max_h, 50);
y_fine_halv = zeros(1, length(x_fine_halv));
for idx = 1:length(x_fine_halv)
    xq = x_fine_halv(idx);
    nodes_idx = choose_four_nodes(x_halv, xq);
    y_fine_halv(idx) = cubicLagrangeScalar(x_halv(nodes_idx), y_halv(nodes_idx,1), xq);
end
[y_max_halv, idx_max_halv] = max(y_fine_halv);
x_max_h_interp = x_fine_halv(idx_max_halv);

abs_trunkeringsfel = abs(y_max_halv - y_max_h);
rel_trunkeringsfel = abs_trunkeringsfel / y_max_halv;
interpolationsfel = max(interp_fel_h, interp_fel_halv);

fprintf('Högsta höjden (y_max):    %.10f\n', y_max_h);
fprintf('x-position för maxhöjd:  %.10f\n', x_max_interp);
fprintf('\nInterpolationsfel (fin):%.10e\n', interp_fel_h);
fprintf('\nTrunkeringsfel (absolut): %.10e\n', abs_trunkeringsfel);
fprintf('Trunkeringsfel (relativt):%.10e\n', rel_trunkeringsfel);
fprintf('\nInterpolationsfel (max):  %.10e\n', interpolationsfel);

figure;
plot(x, y(:,1), 'b-', 'LineWidth', 2);
hold on;
plot(x_max_interp, y_max_h, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('x');
ylabel('y(x)');
title('Lösning med korrigerat K_{0} och beräknad maxhöjd');
grid on;
hold off;

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

% Hjälpfunktion: Välj fyra lämpliga noder från ett sorterat fält x för interpolation vid xq.
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
