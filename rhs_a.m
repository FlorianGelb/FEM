function rhs=a_rhs(t, a, dummy, phi, phixx);
rhs = (1/2)*(phi.')*phixx*a+i*(phi.')*((abs(phi*a).^2).*(phi*a));
end