function rhs=a_rhs(t, a, dummy, phi, phixx, alpha, dt, h);
rhs = phi.'*phixx*a*alpha  + phi.' * h(:, int32(t/dt) + 1);
end