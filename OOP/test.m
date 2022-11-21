T = 1;
L = 1;
n = 50;
dx = L/n;
nt = 10000;
alpha = 1;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
u0 = sin(2*pi/L * X);

m = main(["FT"], L, u0, T, n, nt, alpha);
m = m.solve();
m = m.pod(0.95);
m.plot()