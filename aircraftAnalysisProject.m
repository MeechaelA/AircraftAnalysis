hold off, close all, clear, clc, format long
% Hydra X6 attack bomber

% Test flight 1 - Max velocity data
% Fuel load - 3000 lb, half max
% take-off weight - 15,000 lb
% alt - 20,000 ft
% max vel achieved - Mach 0.68
% range at max vel - 300 mi
% L/D at flight conditions - 6.594

% Test flight 2 - Max endurance
% Fuel load: 3,000 lb, half max
% Take-off weight: 15,000 lb
% Altitude: 20,000 ft
% Max endurance: 11.06 hr
% Flight velocity for max endurance: Mach 0.3233
% L/D at flight conditions: 15.331

% Plot the following:
% (1) Tr and Ta vs velocity at sea level
% (2) Pr and Pa vs velocity at sea level, 15,000 ft, and 20,000 ft
% (3) Rate of climb vs velocity at sea level

% Table with info requested by Fury (range, endurance, ceiling)
% 1) What’s their maximum range, and at what velocity and altitude?
% 2) What’s their maximum endurance, and at what velocity and altitude?
% 3) The helicarrier has an operating altitude of 15,000 ft. If these aircrafts engage the
% helicarrier, what’s the max duration could they operate?
% 4) What’s their best absolute and service ceiling? And at what velocities?
% 5) What’s their max rate of climb at sea level?
%% General Calculations
W = 18000; % [lb] THIS IS FULLY FUELED WEIGHT
span = 50.87; % Tip2tip span [ft]
S = 345.12; % Planform Wing Area [ft^2]
AR = span^2/S; % Aspect Ratio [ft^2/ft^2]
e = 1.78*(1-(0.045*AR^(0.68)))-0.64;

% atmospheric conditions
rhoSL  = 2.377*10^-3;  % [slugs/ft^3]
rho15k = 1.4962*10^-3; % [slugs/ft^3]
rho20k = 1.2673*10^-3; % [slugs/ft^3]

Temp20k = 447.43; % [degR]

% repeated functions
%% Test Flight 1 Calcs (Max Vel)
MTest1 = 0.68;
a20k = sqrt(1.4*1716*Temp20k);
V20kTest1 = MTest1*a20k;

%% Test Flight 2 Calcs (Max Endurance)

MTest2 = 0.3233;
V20kTest2 = MTest2*a20k;
Cl20k = W/(0.5*rho20k*V20kTest2^2*S);
Cd20k = Cl20k/15.331;
Cdi20k = (Cl20k^2)/(pi*e*AR);
Cd0 = Cd20k - Cdi20k;

%% Plots
% (1) Tr and Ta vs velocity at sea level

% 20k - use max velocity at 20k to find max velocity at SL (needed to find Ta)
figure(1)
Tr20k = @(v) W./((W./(((.5)*rho20k*(v.^2)) * S))./(Cd0 + (((W./(((.5)*rho20k*(v.^2)) * S)).^2)/(pi*e*AR))));
fplot(Tr20k, [150 900]), title('Thrust Required and Thrust Available vs Velocity at 20,000 [ft]')
xlabel('Velocity [ft/s]'), ylabel('Thrust [lb_f]'), xlim([0 900]), ylim([0 6000]), grid on, hold on
Ta20k = Tr20k(V20kTest1);
fplot(@(V) Ta20k+(0.*V), [0 900]); plot(V20kTest1, Ta20k, 'r*')
legend('Thrust Required', 'Thrust Available', 'Maximum Velocity at 20,000 [ft]','Location', 'southeast'), hold off
TaSL = Ta20k*(rhoSL/rho20k);
TaSLTable = @(V) TaSL + (0.*V);
% SL
figure(2)
TrSL = @(v) W./((W./(((.5)*rhoSL*(v.^2)) * S))./(Cd0 + (((W./(((.5)*rhoSL*(v.^2)) * S)).^2)/(pi*e*AR))));
fplot(TrSL, [150 900]), title('Thrust Required and Thrust Available vs Velocity at Sea Level')
xlabel('Velocity [ft/s]'), ylabel('Thrust [lb_f]'), xlim([0 900]), ylim([0 6000]), grid on, hold on
fplot(@(V) TaSL+(0.*V), [0 900])
SLdataPoint = [722.78:.001:722.79; TaSLTable(722.78:.001:722.79); TrSL(722.78:.001:722.79);]'
plot(722.789, TrSL(722.789), 'r*')
legend('Thrust Required', 'Thrust Available', 'Maximum Velocity at Sea Level','Location', 'southeast'), hold off

% (2) Pr and Pa vs velocity at sea level, 15,000 ft, and 20,000 ft
figure(3)
PrSL = @(V) TrSL(V)*V;
PaSL = @(V) TaSL*V;
fplot(PrSL, [60 850]), hold on, grid on
fplot(PaSL, [0 850])
%% Questions and table data

% 1) What's their maximum range, and at what velocity and altitude?

% Cl12Cd = (((1/3)*Cd0*pi*e*AR)^(1/4))/((4/3)*Cd0);
% R = (2*sqrt((2)/(rho8km*S))*(1/ct)*(Cl12Cd)*((Wi.^(1/2))-(Wf.^(1/2))))/1000;