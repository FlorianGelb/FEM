clear;clc;close all
tic
N = 50;
dt = 0.0001
X = 1;
T = 1;
K = zeros(N);
M = zeros(N);
delta_x = X / (N-1);
syms t(x) t2(x) t3(x) f(x) p(x)
t3(x) = triangularPulse(0, delta_x, 2 * delta_x, x);
t(x) = triangularPulse(0, delta_x, 2 * delta_x, x)^2;
t2(x)= triangularPulse(0, delta_x, 2 * delta_x, x) * triangularPulse(delta_x, 2 * delta_x, 3* delta_x, x);
ii = integral(matlabFunction(t(x)), 0, X);
ij = integral(matlabFunction(t2(x)), 0, X);
f(x) = 100*exp(-0.5*(x-5)) + exp(0.5*(x-X+15));
F = zeros(N, 1);
p(x) = 0
H(1) = f(0);
H(end) = f(X);
alpha = 0.0005;

for i = 1:N-2
  F(i+1) =   integral(matlabFunction(f(x) * triangularPulse((i-1)*delta_x, (i)*delta_x, (i+1)*delta_x, x)), 0, X);
end

F(1) = f(0);
F(end) = f(X);


for i = 1:N
    K(i,i) = -2/delta_x;
    M(i, i) = ii;
    if i-1 > 0
        K(i, i-1) = 1/delta_x;
        M(i, i-1) =  ij;
    end
    if i + 1 < N+1
        K(i, i+1) = 1/delta_x;
        M(i, i+1) =   ij;
    end
end

K = alpha * K;
Mmod = M;
Mmod(1,:) = [1, zeros(1, N-1)];
Mmod(N,:) = [zeros(1, N-1), 1];
C = zeros(N, uint8(1+(T/dt)));
c_o = inv(Mmod)*F;

C(:, 1) = c_o;
j = 2;
im = inv(M);
for i=dt:dt:T
    c_n = alpha*dt * im * K * c_o + c_o;
    c_o = c_n;
    C(:, j) = c_o;
    j= j +1;
end


dx = 0.001;
S = []
for t = 1:1:T/dt
c = C(:, t);
interp_c = interp1([0:delta_x:X], c, [0:dx:X-dx]);
S(:,t) = interp_c;
end

[U, Z, V] = svd(S,"econ");
energies = diag(Z) / trace(Z);
e_t = 0.90;
e = 0;
trunc = 0;
for i = 1:size(energies, 1)
e = e+energies(i);
if e >= e_t
trunc = i;
    break
end
end

U2 = U(:, 1:trunc);

phi = [];
%for i = 1:trunc
%phi(:, i) = Nderiv2(Nderiv2(U2(:, i), dx), dx);
%end

k = (2*pi/X) * [0:X/dx/2-1 -X/dx/2:-1].';
for i = 1:trunc
    phi(:,i) = ifft(-(k.^2).*fft(U2(:, i)));
end

a = U2.'*S(:, 1);
a0 = a;
A = []
A(1, :) = U2*a;
j = 2;
Y = U2.'*phi;

%for i=dt:dt:T
%    a_n = dt * alpha*Y*a + a;
 %   a = a_n;
 %   A(j, :) = U2*a;
  %  j= j +1;
%end
A = A.';

[t, asol] = ode45("a_rhs", [0:dt:T], a0,[], U2, phi, alpha);
usol = zeros(size(U2, 1), length(t)).';
for j = 1:length(t)
    for i = 1:trunc;
        usol(j,:) = usol(j, :) + asol(j,i)*U2(:,i).';
    end
end
A = usol.';
figure(1)
imagesc(C)
colormap turbo
colorbar
figure(2)
imagesc(S);
colormap turbo
colorbar
figure(3)
imagesc(A)
colormap turbo
colorbar
figure(4)
plot(energies);
toc


function dy = Nderiv2(y,h)
% Compute the second derivative of input vector y, with spacing h
n = length(y);
for i=1:n;
    switch i
        case 1
            % use Forward difference equation for y''
            dy(i) =  (y(i+1) - y(i)) / h;
        case n
            % use backward difference equation for y''
             dy(i) = (y(i) - y(i-1))/h;
        otherwise
            % use central difference equation for y''
            dy(i) = (y(i+1) - y(i-1)) / 2*h;
            
end
end
end

