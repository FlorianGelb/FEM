T = 0.1;
L = 1;
n = 100;
dx = L/(n-1);
nt = 150;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.01;
u0 = sin(2*pi/L * X);
h = randi([100 500], [n,nt]);
h(:, 100:nt) = -500*eye(n, nt-100+1);
figure(200)
imagesc(h);

m = main(["FEM" "FT"], L, u0, h, T, n, nt, alpha);
m = m.solve();
m = m.pod(0.9995);
m.plot();