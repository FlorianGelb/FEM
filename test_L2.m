tic
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
%h = eye(n, nt);


E_2 = {};
m = main(["FEM"], L, u0, h, T, n, nt, alpha);
m = m.solve();


p = parpool("threads");

F1 = parfeval(@loop, 1, 1, 30, m);
F2 = parfeval(@loop, 1, 31, 51, m);
F3 = parfeval(@loop, 1, 52, 69, m);
F4 = parfeval(@loop, 1, 70, 84, m);
F5 = parfeval(@loop, 1, 85, 94, m);
F6 = parfeval(@loop, 1, 95, 100, m);

R1 = fetchOutputs([F1, F2, F3, F4, F5, F6], 'UniformOutput', false);
delete(gcp('nocreate'))
E_2 = cat(2, R1{:});
E_POD = [];
E_BT = [];
E_MT = [];
E_HNA = [];
for j = 1:length(E_2)
    E_POD(j, :) = E_2{j}("Proper Orthogonal Decomposition");

    E_BT(j, :) = E_2{j}("Balanced Reduction");
    
    E_MT(j, :) = E_2{j}("Modal Truncation");
    

    E_HNA(j, :) = E_2{j}("Hankel Norm Approximation");

     if E_2{j}("Hankel Norm Approximation") < 1e+20
         
     else
         E_HNA(end+1) = NaN;
     end
end

figure;
plot(E_POD(1:end-1))
set(gca, 'FontSize', 14);
title("$$||\epsilon(r)||_{L2}$$ Proper Orthogonal Decomposition", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||Y - Ŷ||_2");
%exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_POD_SIN.png");
%exportgraphics(gcf, "/home/f/Documents/Studienarbeit/images/L2_POD_SIN.png");


figure;
plot(E_MT)
set(gca, 'FontSize', 14);
title("$$||\epsilon(r)||_{L2}$$ Modal Truncation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||Y - Ŷ||_2");
%exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_MT.png");
%exportgraphics(gcf, "/home/f/Documents/Studienarbeit/images/L2_MT_SIN.png");

figure;
plot(E_BT)
set(gca, 'FontSize', 14);
title("$$||\epsilon(r)||_{L2}$$ Balanced Truncation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||Y - Ŷ||_2");
%exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_BT.png");
%exportgraphics(gcf, "/home/f/Documents/Studienarbeit/images/L2_BT_SIN.png");

figure;
plot(E_HNA);
set(gca, 'FontSize', 14);
title("$$||\epsilon(r)||_{L2}$$ Hankel Norm Approximation", ' ', "interpreter", "latex")
xlabel("r");
ylabel("||Y - Ŷ||_2");
%exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_HNA.png");
%exportgraphics(gcf, "/home/f/Documents/Studienarbeit/images/L2_HNA_SIN.png");

figure;
boxplot([log10(E_POD), log10(E_BT), log10(E_MT), log10(E_HNA)], "Labels", {'POD', 'BT', 'MT', 'HNA'})
set(gca, 'FontSize', 14);
ylabel("$\log_{10}(\epsilon)$", "interpreter", "latex");
title("Comparison of $L_{2}$ error", ' ', "interpreter", "latex")
%exportgraphics(gcf, "C:/Users/Florian/Documents/Studienarbeit/images/L2_BOX.png");
%exportgraphics(gcf, "/home/f/Documents/Studienarbeit/images/L2_BOX_SIN.png");
% %m.plot();
toc

function E_2 = loop(a, b, m)
    E_2 = {}
    for i = a:b
        if mod(i, 10) == 0
            disp(i);
        end
        m = m.flush_rom();
        m = m.rom(i, false);
        e2 = m.compare("L2");
        E_2{end+1} = e2
    end
end