clear all; clc; close all;

L = 30; n = 512;
x2 = linspace(-L/2, L/2, n+1); x=x2(1:n);
k = (2*pi/L) * [0:n/2-1  -n/2:-1].';
u = sech(x); ut = fft(u);
t = linspace(0, 2*pi, 41);
[t, utsol] = ode45("nls_rhs", t, ut, [], k);
for j = 1:length(t)
usol(j,:) = ifft(utsol(j,:));    
end

surfl(x, t, abs(usol));
X = usol.';
[u,s, v] = svd(X, "econ");

phi = u(:, 1:2); 
a0 = [phi(:,1).' * (2*(sech(x))).'
    phi(:,2).' * (2*sech(x)).'];

mode1_xx = ifft(-(k.^2).*fft(u(:, 1)));
mode2_xx = ifft(-(k.^2).*fft(u(:, 2)));
phixx = [mode1_xx mode2_xx];
[t, asol]=ode45("rhs_a", t, a0, [], phi, phixx);

for j = 1:length(t)
usol(j,:) = asol(1)*phi(:, 1) + asol(2)*phi(:, 2); 
end
figure(6);
surfl(abs(usol));