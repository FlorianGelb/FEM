%-----Unsteady 1-D Heat Transfer Code-----
%-----All Units are in S.I.---------------
clear all
clc
%Initialization of parameters
L = 0.3; % Thickness of wall
Nx = 200; % Discretisation
x(1) = 0;
x(Nx) = L;
time = 86400; %Time domain 
Nt = 4*86400; % Time steps
dx = L/Nx;
dt = time/Nt;
Tb = 27;
h2 = 8;      % Coefficient of heat transfer
h1 = 28;     % Coefficient of heat transfer 
den = 2000;  % Material Density
c = 1;        
k = 2;       % Thermal conductivity of wall
t = [0:dt:time];
Tw = 35 + 5*cos(2*pi*(t-16*3600)/(24*3600)); %Temperature variation at one side of the wall
%---Initialisation of matrices------
T = zeros(Nx, Nt);
ap1 = zeros(Nx, Nt);
aw = zeros(Nx, Nt);
ae = zeros(Nx,Nt);
apo = zeros(Nx, Nt);
for i = 1:Nx
 T(i,1) = 30;
end
%------Calculation of Matrix Coefficients and Solution of Linear System------------
for j = 1:Nt
    for i = 2:Nx-1
        ap1(i,j+1) = (den*c*dx/dt + 2*k/dx);
        aw(i,j+1) = -k/dx;
        ae(i,j+1) = -k/dx;
        apo(i,j+1) = (den*c*dx/dt)*T(i,j);
        x(i) = x(i-1) + dx;
    end
    ap1(1,j+1) = (den*c*dx/dt +(1/(1/h1 + (dx/2*k))) + k/dx);
    ae(1,j+1) = -k/dx;
    apo(1,j+1) = (den*c*dx/dt)*T(1,j)+Tw(j+1)*(1/(1/h1 + dx/(2*k)));
    ap1(Nx,j+1) = (den*c*dx/dt + 1/(1/h2 + dx/(2*k)) + k/dx);
    aw(Nx,j+1) = -k/dx;
    apo(Nx,j+1) = (den*c*dx/dt)*T(Nx,j)+Tb*(1/(1/h2 + dx/(2*k)));
    T(:,j+1) = tridiag(ap1(:,j+1), aw(:,j+1), ae(:,j+1), apo(:,j+1));   
        
end

%-----------Results Plotting------------
figure(1)
plot(t, T(Nx/2,:), 'red', t, Tw, 'blue', 'LineWidth',2.0)
grid on
figure(2)
plot(x, T(:,3600), 'red', x, T(:, 7200), 'blue', x, T(:,10800), 'green', 'LineWidth',2.0)
grid on
figure(3)
plot(x, T(:,10),'red', x, T(:,50),'blue', x, T(:,90),'green','LineWidth',2.0)
grid on

