T = 5;
L = 1;
n = 100;
dx = L/(n);
nt = 25000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.01;
u0 = (-4*sin(2*pi/L * X))+1;
h = randi([-1000 1000], [n,nt]);

%h = eye(n, nt);


m = main(["FEM"], L, u0, h, T, n, nt, alpha);
m = m.solve();
m = m.rom(1, 50);
%m.compare();
m.plot();