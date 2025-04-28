function [x, y, interpfel] = Rung2(fprim, x0, x_slut, u0, h, K0, K1)
    % RK4-solver med fast steglängd
    N = floor((x_slut - x0)/h);
    x = x0 + (0:N)*h;
    y = zeros(N+1, length(u0));
    y(1,:) = u0';
    
    for i = 1:N
        k1 = fprim(x(i), y(i,:)', K0, K1);
        k2 = fprim(x(i)+h/2, y(i,:)' + (h/2)*k1, K0, K1);
        k3 = fprim(x(i)+h/2, y(i,:)' + (h/2)*k2, K0, K1);
        k4 = fprim(x(i)+h, y(i,:)' + h*k3, K0, K1);
        y(i+1,:) = (y(i,:)' + (h/6)*(k1 + 2*k2 + 2*k3 + k4))';
    end

    % --- Handgjord kubisk interpolation (med Lagrange-polynom) ---
    % Interpolera med alla punkter
    x_interp = x0 : h/2 : x_slut;
    y_interp = zeros(length(x_interp), size(y,2));

    for idx = 1:length(x_interp)
        xq = x_interp(idx);
        nodes_idx = choose_four_nodes(x, xq);
        x_nodes = x(nodes_idx);
        for j = 1:size(y,2)
            y_interp(idx, j) = cubicLagrangeScalar(x_nodes, y(nodes_idx, j), xq);
        end
    end

    % Interpolera med varannan punkt (använd bara varje 2:a nod)
    x_coarse = x(1:2:end);
    y_coarse = y(1:2:end, :);
    y_interp_varannan = zeros(length(x_interp), size(y,2));
    
    for idx = 1:length(x_interp)
        xq = x_interp(idx);
        nodes_idx = choose_four_nodes(x_coarse, xq);
        x_nodes = x_coarse(nodes_idx);
        for j = 1:size(y,2)
            y_interp_varannan(idx, j) = cubicLagrangeScalar(x_nodes, y_coarse(nodes_idx, j), xq);
        end
    end

    % Beräkna interpolationsfelet som största absoluta skillnaden mellan de två metoderna
    felmatris = abs(y_interp - y_interp_varannan);
    interpfel = max(felmatris(:));
end

function val = cubicLagrangeScalar(x_nodes, y_nodes, xq)
    % Handgjord kubisk Lagrange-interpolering för en enskild skalfunktion.
    % x_nodes: Vektor med fyra noder, y_nodes: motsvarande funktionsvärden.
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

function idx = choose_four_nodes(x, xq)
    % Välj fyra index från den sorterade vektorn x som täcker xq på ett bra sätt.
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

