function rhs = rhs_1dHeatFT(t, k, a, u, dt, h)
rhs = -a .* k.^2.* u + fft(h(:, int32(t/dt) + 1));
end

