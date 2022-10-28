tic
N = 100;
dt = 0.00001
X = 1;
T = 0.1;
K = zeros(N);
M = zeros(N);
delta_x = X / (N-1);
syms t(x) t2(x) t3(x) f(x) p(x)
t3(x) = triangularPulse(0, delta_x, 2 * delta_x, x);
t(x) = triangularPulse(0, delta_x, 2 * delta_x, x)^2;
t2(x)= triangularPulse(0, delta_x, 2 * delta_x, x) * triangularPulse(delta_x, 2 * delta_x, 3* delta_x, x);
ii = integral(matlabFunction(t(x)), 0, X);
ij = integral(matlabFunction(t2(x)), 0, X);
f(x) = exp(-0.5*(x-5)) + exp(0.5*(x-N+15));
F = zeros(N, 1);
p(x) = 0
H = zeros(N, 1);
H(1) = f(0);
H(end) = f(X);
alpha = 0.1

for i = 1:N-2
  F(i+1) =   integral(matlabFunction(f(x) * triangularPulse((i-1)*delta_x, (i)*delta_x, (i+1)*delta_x, x)), 0, X);
end

F(1) = f(1);
F(end) = f(N);


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
C = zeros(N, uint8(1+(1/dt)));
c_o = inv(Mmod)*F;

C(:, 1) = c_o;
j = 2;
im = inv(M);
for i=dt:dt:T
    c_n = dt * im * K * c_o + c_o;
    c_o = c_n;
    C(:, j) = c_o;
    j= j +1;
end
p(x) = 0;
c = C(:, 15);
%for t = 1:dt:1
%for i = 0:N-1
%  hold on
%  p = p + c(i+1) * triangularPulse((i-1)*delta_x, (i)*delta_x, (i+1)*delta_x, x);
% 
%end
%end


colorbar

imagesc(C)
colormap turbo
colorbar
toc
