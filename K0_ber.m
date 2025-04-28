function [K0_opt, skillnad_nu, interpfel] = K0_ber(x0, x_slut, u0, h, K1, lutning, tolerans)
    % K0_ber beräknar det optimala värdet på K0 så att den slutliga
    % lutningen blir lika med det önskade värdet (lutning) inom angiven tolerans.
    %
    % Parametrar:
    %   x0      - Startvärde för x
    %   x_slut  - Slutvärde för x
    %   u0      - Vektor med startvärden [y(0); y'(0)]
    %   h       - Steglängd
    %   K1      - Parameter K1
    %   lutning - Önskad slutlig lutning (t.ex. -0.51)
    %   tolerans- Konvergenstolerans för skillnaden i lutning
    %
    % Returnerar:
    %   K0_opt      - Det optimerade värdet på K0
    %   skillnad_nu - Slutlig differens mellan beräknad och önskad lutning
    %   interpfel   - Det maximala interpolationsfelet (returnerat från Rung2)
    
    % Startgissningar
    K0_innan = 0;
    [~, y] = Rung2(@fprim, x0, x_slut, u0, h, K0_innan, K1);
    lutningSlut_innan = y(end, 2);
    skillnad_innan = lutningSlut_innan - lutning;
    
    % Nästa gissning
    K0_nu = K0_innan + 0.1;
    [~, y] = Rung2(@fprim, x0, x_slut, u0, h, K0_nu, K1);
    lutningSlut_nu = y(end, 2);
    skillnad_nu = lutningSlut_nu - lutning;
    
    % Newton-iteration
    while abs(skillnad_nu) > tolerans
        deltaK0 = -skillnad_nu * (K0_nu - K0_innan) / (skillnad_nu - skillnad_innan);
        K0_next = K0_nu + deltaK0;
        
        K0_innan   = K0_nu;
        skillnad_innan = skillnad_nu;
        K0_nu     = K0_next;
        
        [~, y, interpfel] = Rung2(@fprim, x0, x_slut, u0, h, K0_nu, K1);
        lutningSlut_nu = y(end, 2);
        skillnad_nu = lutningSlut_nu - lutning;
    end
    
    K0_opt = K0_nu;
end
