clc; clear all; close all;
T = 1;
L = 1;
n = 50;
dx = L/n;
nt = 10000;
alpha = 0.5;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
u_0 = sin(2*pi/L * X);
%u_0 = 0*X;
% u_0((L/2 - L/10)/dx:(L/2 + L/10)/dx) = 1;

S = []
[t, sol] = ode45(@(t,sol)rhs_1dHeatFT(t,k, alpha, sol),t,fft(u_0));
%[t, sol] = ode45("rhs_1dHeatFT", t, k, [], alpha, u_0);

usol = [];
for j = 1:length(t)
usol(j, :) = ifft(sol(j,:));    
end
imagesc(usol.');
colorbar
colormap turbo


function rhs = f_rhs(t, k, a, u)
size(u)
rhs = -a^2* k.^2.* u;
end
