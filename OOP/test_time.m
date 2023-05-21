delete(gcp('nocreate'));
T = 1;
L = 1;
n = 100;
dx = L/(n);
nt = 10000;
k = (2*pi/L) * [0:n/2-1 -n/2:-1].';
X = linspace(0, L, n);
t = linspace(0, T, nt);
alpha = 0.1;
u0 = 0*(sin(2*pi/L * X)) + 1;
h = 0*randi([-10 10], [n,nt]);


delete(gcp('nocreate'))
E_2 = {};
m = main(["FEM"], L, u0, h, T, n, nt, alpha);
m = m.solve();

for i = 1:100
    if mod(i, 10) == 0
        disp(i);
    end
    m = m.flush_rom();
    m = m.rom(i, true);
    e2 = m.compare("time");
    E_2{end+1} = e2;
 end
E_POD = [];
E_BT = [];
E_MT = [];
E_HNA = [];
for j = 1:length(E_2)
    temp = E_2{j}("Proper Orthogonal Decomposition");
    E_POD(j, :) = temp{1};

    temp = E_2{j}("Balanced Reduction");
    E_BT(j, :) = temp{1};
    
    temp = E_2{j}("Modal Truncation");
    E_MT(j, :) = temp{1};
    
    temp = E_2{j}("Hankel Norm Approximation");
    E_HNA(j, :) = temp{1};

end

figure
boxplot(E_POD.');
set(gca, 'FontSize', 18);
ylabel("Time in s", "interpreter", "latex");
xlabel("r");
title("Time POD", ' ', "interpreter", "latex")
xticks([0 10 20 30 40 50 60 70 80 90 100])
xticklabels({0 10 20 30 40 50 60 70 80 90 100})
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/time/POD.png")


figure
boxplot(E_BT.');
set(gca, 'FontSize', 18);
ylabel("Time in s", "interpreter", "latex");
xlabel("r");
title("Time BT", ' ', "interpreter", "latex")
xticks([0 10 20 30 40 50 60 70 80 90 100])
xticklabels({0 10 20 30 40 50 60 70 80 90 100})
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/time/BT.png")


figure
boxplot(E_MT.');
set(gca, 'FontSize', 18);
ylabel("Time in s", "interpreter", "latex");
xlabel("r");
title("Time MT", ' ', "interpreter", "latex")
xticks([0 10 20 30 40 50 60 70 80 90 100])
xticklabels({0 10 20 30 40 50 60 70 80 90 100})
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/time/MT.png")

figure
boxplot(E_HNA.');
set(gca, 'FontSize', 18);
ylabel("Time in s", "interpreter", "latex");
xlabel("r");
title("Time HNA", ' ', "interpreter", "latex")
xticks([0 10 20 30 40 50 60 70 80 90 100])
xticklabels({0 10 20 30 40 50 60 70 80 90 100})
exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/time/HNA.png")
