T = 0.1;
L = 1;
n = 50;
dx = L/n;
nt = 10000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 1
u0 = 0*sin(2*pi/L * X) + 1;

m = main(["FT"], L, u0, T, n, nt, alpha);
m = m.solve();
m = m.pod(0.95);
m.plot()