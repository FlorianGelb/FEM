T = 1;
L = 1;
n = 100;
dx = L/(n);
nt = 10000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.1;
u0 = (sin(2*pi/L * X));
h = 0*randi([-10 10], [n,nt]);
%h = eye(n, nt);



m = main(["FEM"], L, u0, h, T, n, nt, alpha);
m = m.solve();
m = m.rom(100, 1);
m.compare(false);
m.plot()

E_2 = {};
for i = 1:n
    if mod(i, 10) == 0
        disp(i);
    end
    m = main(["FEM"], L, u0, h, T, n, nt, alpha);
    m = m.solve();
    m = m.rom(i);
    e2 = m.compare(false);
    E_2{end+1} = e2;
end

E_POD = [];
E_BT = [];
E_MT = [];
E_HNA = [];
for j = 1:length(E_2)
    E_POD(end+1) = E_2{j}("Proper Orthorgonal Decomposition");
    E_BT(end+1) = E_2{j}("Balanced Reduction");
    E_MT(end+1) = E_2{j}("Modal Truncation");
    E_HNA(end+1) = E_2{j}("Hankel Norm Approximation")
end
figure;
plot(E_POD)
set(gca, 'FontSize', 14);
title("$$||\epsilon(r)||_{L2}$$ Proper Orthogonal Decomposition", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||Y - Ŷ||_2");
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_POD_SIN.png");

figure;
plot(E_MT)
set(gca, 'FontSize', 14);
title("$$||\epsilon(r)||_{L2}$$ Modal Truncation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||Y - Ŷ||_2");
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_MT_SIN.png");

figure;
plot(E_BT)
set(gca, 'FontSize', 14);
title("$$||\epsilon(r)||_{L2}$$ Balanced Truncation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||Y - Ŷ||_2");
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_BT_SIN.png");

figure;
plot(E_HNA)
set(gca, 'FontSize', 14);
title("$$||\epsilon(r)||_{L2}$$ Hankel Norm Approximation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||Y - Ŷ||_2");
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_HNA_SIN.png");

%m.plot();