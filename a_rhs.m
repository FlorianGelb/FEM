function rhs=a_rhs(t, a, dummy, phi, phixx, alpha);
rhs = phi.'*phixx*a*alpha;
end