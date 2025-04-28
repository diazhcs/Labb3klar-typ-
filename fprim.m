function uprim = fprim(x,u,K0,K1)
    uprim = zeros(2,1);
    uprim(1) = u(2);
    uprim(2) = -(K0 - K1*x)*u(1)*(1 + (u(2))^2)^(3/2);
end