del = -15.5e6; %Detuning
L = -del/3; %Linewidth of Laser
s0 = 11; %Saturation Intensity
lam = 852.3e-9; %Laser Wavelength
I = 800; %Intensity of trapping laser
R = L*sqrt(1+I/s0); %Rabi Frequency
s = R^2/(del^2+(L/2)^2); %Saturation Parameter
k = 2*pi/lam; %Wave Number of Laser
Current = 1; %Current through magnetic coils
%Grad = 0.1; %Test Gradient
Grad = 7.5*Current/100; %Gradient produced by coils
mu = 9.27e-24; %Bohr Magneton
hbar = 1.055e-34; %hbar
Na = 6.022e23; %Avogadro's Number
M = 2.207e-25; %Mass of Cesium atom

%Damping Coefficient
alpha = -2*hbar*k^2*(s/(1+s)^2)*(del*L/(del^2+L^2/4));

%Spring Constant
kappa = (mu/hbar)*(Grad/k)*alpha;

%Damping Constant
gamma = alpha/(2*M);

%Oscillation Frequency

omega = sqrt(kappa/M);

%Restoring time

tau = 2*gamma/omega^2;