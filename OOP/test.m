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



%m = main(["FEM"], L, u0, h, T, n, nt, alpha);
%m = m.solve();
%m = m.rom(100, 1);
%m.compare(false);
%m.plot()

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
plot(E_POD(1:end-1))
set(gca, 'FontSize', 14);
title("$$||G_{\epsilon(r)}||_{H_{\infty}}$$ Proper Orthogonal Decomposition", ' ', "interpreter", "latex")
%title("$$||\epsilon(r)||_{L2}$$ Proper Orthogonal Decomposition", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||G - G_r||_{H_{\infty}}");
%ylabel("||Y - Ŷ||_2");
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/freq/H_POD_SIN.png");

figure;
plot(E_MT)
set(gca, 'FontSize', 14);
title("$$||G_{\epsilon(r)}||_{H_{\infty}}$$ Modal Truncation", ' ', "interpreter", "latex")
%title("$$||\epsilon(r)||_{L2}$$ Modal Truncation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||G - G_r||_{H_{\infty}}");
%ylabel("||Y - Ŷ||_2");
%exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/freq/H_MT.png");

figure;
plot(E_BT)
set(gca, 'FontSize', 14);
title("$$||G_{\epsilon(r)}||_{H_{\infty}}$$ Balanced Truncation", ' ', "interpreter", "latex")
%title("$$||\epsilon(r)||_{L2}$$ Balanced Truncation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||G - G_r||_{H_{\infty}}");
%ylabel("||Y - Ŷ||_2");
%exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/freq/H_BT.png");

figure;
plot(E_HNA(1:int8(n/3)+1))
set(gca, 'FontSize', 14);
title("$$||G_{\epsilon(r)}||_{H_{\infty}}$$ Hankel Norm Approximation", ' ', "interpreter", "latex")
%title("$$||\epsilon(r)||_{L2}$$ Hankel Norm Approximation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||G - G_r||_{H_{\infty}}");
%ylabel("||Y - Ŷ||_2");
%exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/freq/H_HNA.png");

E_POD_1 = load("E_POD.mat").E_POD;
figure;

boxplot([log10(E_POD_1); log10(E_POD); log10(E_BT); log10(E_MT); log10(E_HNA)].', "Labels", {'POD', 'POD-SIN', 'BT', 'MT', 'HNA'})
set(gca, 'FontSize', 14);
ylabel("$\log_{10}(\epsilon)$", "interpreter", "latex");
title("Comparison of error", ' ', "interpreter", "latex")
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/freq/H_BOX.png");
%m.plot();