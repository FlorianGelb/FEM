T = 1;
L = 1;
n = 70;
dx = L/(n);
nt = 10000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.01;
u0 = 0*(-4*sin(2*pi/L * X)) + 10000;
h = 0*randi([-1000 1000], [n,nt]);

%h = eye(n, nt);


m = main(["FEM"], L, u0, h, T, n, nt, alpha);
m = m.solve();
m = m.rom(0.90, 10);
m.compare(true);
m.plot();