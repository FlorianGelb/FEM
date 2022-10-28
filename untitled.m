tic
N = 9;
X = 1;
K = zeros(N);
M = zeros(N);
delta_x = X / (N-1);
a = 1;
syms p(x) G(x) x
G(x) = cos(x);
F = zeros(N, 0);
p(x) = 0
syms t(x) t2(x) t3(x)
t3(x) = triangularPulse(0, delta_x, 2 * delta_x, x);
t(x) = triangularPulse(0, delta_x, 2 * delta_x, x)^2;w
t2(x)= triangularPulse(0, delta_x, 2 * delta_x, x) * triangularPulse(delta_x, 2 * delta_x, 3* delta_x, x);
F(1) = 1;
F(N) = 5;
for i = 1:N-2
  F(i+1) =   integral(matlabFunction(G(x) * triangularPulse((i-1)*delta_x, (i)*delta_x, (i+1)*delta_x, x)), 0, X);
%  hold on
%  fplot(triangularPulse((i-1)*delta_x, (i)*delta_x, (i+1)*delta_x, x), [0, X])
%  hold off
end
%stop
for i = 1:N
    K(i,i) = -2/delta_x;
    M(i, i) = a * integral(matlabFunction(t(x)), 0, X);
    if i-1 > 0
        K(i, i-1) = 1/delta_x;
        M(i, i-1) =  a * integral(matlabFunction(t2(x)), 0, X);
    end
    if i + 1 < N+1
        K(i, i+1) = 1/delta_x;
        M(i, i+1) =  a * integral(matlabFunction(t2(x)), 0, X);
    end
end
M(1,:) = [1/2, zeros(1, N-1)]
M(N,:) = [zeros(1, N-1), 1/2]
K(1,:) = [1/2, zeros(1, N-1)]
K(N,:) = [zeros(1, N-1), 1/2]
%[Z, D, V] = svd(K+M);
%c = (transpose(V) * inv(D) * transpose(Z)) * transpose(F)
c = inv(transpose((K+M)) * (K+M))*transpose(K+M) * transpose(F)
%c = pinv(K + M) * transpose(F);
j = 1;
for i = 0:N-1
%F(i+1) =   integral(matlabFunction(G(x) * triangularPulse((i-1)*delta_x, (i)*delta_x, (i+1)*delta_x, x)), 0, X);
  hold on
  p = p + c(i+1) * triangularPulse((i-1)*delta_x, (i)*delta_x, (i+1)*delta_x, x);
  fplot(c(i+1) * triangularPulse((i-1)*delta_x, (i)*delta_x, (i+1)*delta_x, x), [0, X]);
end
%for i = 0:delta_x:X-delta_x
% disp(i);
% p = p + c(j) * triangularPulse(i-delta_x, i, i + delta_x, x);
% hold on
% fplot(c(j) * triangularPulse(i-delta_x, i, i + delta_x, x), [0, X]);
% j = j + 1;
%end
syms u(x) U(x)
ode = diff(diff(u,x), x) + a * u(x) == G(x);
U(x) = dsolve(ode, u(0) == 1, u(X) == 5);
fplot(p, [0, X])
fplot(U(x), [0, X])
integral(matlabFunction(abs(U(x)) - abs(p(x))), 0, X)
toc
