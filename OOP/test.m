T = 5;
L = 1;
n = 50;
dx = L/(n-1);
nt = 100000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.05;
u0 = 0*sin(2*pi/L * X) + 1;


m = main(["FEM"], L, u0, T, n, nt, alpha);
m = m.solve();
m = m.pod(0.999955);
m.plot();

