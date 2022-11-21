function rhs = rhs_1dHeatFT(t, k, a, u)
    %u = fft(u);
    rhs = -a^2* k.^2.* u;
end

