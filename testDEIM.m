clear all; close all; clc;

L = 40;
n = 256;
x2 = linspace(-L/2, L/2, n+1);
x = x2(1:n);
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
t = linspace(0, 2*pi,61);

N = 2;
u0 = N * sech(x).';
ut = fft(u0);
[t, utsol] = ode45("ch_pod_sol_rhs", t, ut, [], k);
for j = 1:length(t)
    usol(j, :) = ifft(utsol(j,:));
end
