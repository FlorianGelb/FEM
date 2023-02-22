T = 1;
L = 1;
n = 50;
dx = L/(n-1);
nt = 1000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.25;
u0 = sin(2*pi/L * X);
h = randi([1 5], [n,nt]);
%h = eye(n, nt);
h(1, :) = 0;
h(n, :) = 0;
figure(200)
imagesc(h);


m = main(["FEM" "FT"], L, u0, h, T, n, nt, alpha);
m = m.solve();
m = m.pod(0.95);
m.plot();