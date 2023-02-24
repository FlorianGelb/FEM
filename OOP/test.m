T = 1;
L = 1;
n = 50;
dx = L/(n-1);
nt = 4000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.01;
u0 = -4*sin(2*pi/L * X);
h = randi([-10 10], [n,nt])/10;
%h = eye(n, nt);
figure(200)
imagesc(h);


m = main(["FEM" "FT"], L, u0, h, T, n, nt, alpha);
m = m.solve();
m = m.pod(0.95);
m.plot();