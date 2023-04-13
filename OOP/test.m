T = 1;
L = 1;
n = 70;
dx = L/(n-1);
nt = 10000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.1;
u0 = 0*(-4*sin(2*pi/L * X));
h = randi([-10000 10000], [n,nt]);
%h = eye(n, nt);
figure(20)
imagesc(h)



m = main(["FEM"], L, u0, h, T, n, nt, alpha);
m = m.solve();
m = m.rom(0.95, 15);
m.plot();