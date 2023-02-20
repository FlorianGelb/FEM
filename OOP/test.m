T = 0.1;
L = 1;
n = 100;
dx = L/(n-1);
nt = 150;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.01;
u0 = 0*sin(2*pi/L * X);
h = zeros(n,nt);
h(:, 1:nt) = 1;%5*eye(n, nt-100+1);


m = main(["FEM" "FT"], L, u0, h, T, n, nt, alpha);
m = m.solve();
m = m.pod(0.9995);
m.plot();