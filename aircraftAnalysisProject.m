clear, clc, format long
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
% 1) What�s their maximum range, and at what velocity and altitude?
% 2) What�s their maximum endurance, and at what velocity and altitude?
% 3) The helicarrier has an operating altitude of 15,000 ft. If these aircrafts engage the
% helicarrier, what�s the max duration could they operate?
% 4) What�s their best absolute and service ceiling? And at what velocities?
% 5) What�s their max rate of climb at sea level?
%% General Calculations
W = 18000; % [lb] THIS IS FULLY FUELED WEIGHT
span = 50.87; % Tip2tip span [ft]
S = 345.12; % Planform Wing Area [ft^2]
AR = span^2/S; % Aspect Ratio [ft^2/ft^2]
e = 1.78*(1-(0.045*AR^(0.68)))-0.64;

% atmospheric conditions
rhoSL  = 2.3769*10^-3;  % [slugs/ft^3]
rho15k = 1.4962*10^-3; % [slugs/ft^3]
rho20k = 1.2673*10^-3; % [slugs/ft^3]

T20k = 447.43; % [degR]
%% Test Flight 2 Calcs

MTest2 = 0.3233;
a20k = sqrt(1.4*1716*T20k);
V20kTest2 = MTest2*a20k;
Cl20k = W/(0.5*rho20k*V20kTest2^2*S);
Cd20k = Cl20k/15.331;
Cdi20k = (Cl20k^2)/(pi*e*AR);
Cd0 = Cd20k - Cdi20k;

%% Plots
% (1) Tr and Ta vs velocity at sea level
% SL
% ClSL = W/(.5*rhoSL*V^2*S)
% CdSL = Cd0 + ((ClSL^2)/(pi*e*AR));

TrSL = @(V) W/((W/(.5*rhoSL*V^2*S))/(Cd0 + (((W/(.5*rhoSL*V^2*S))^2)/(pi*e*AR))));
fplot(TrSL, [0 400], 'b-'), xlim([0 350])
% (2) Pr and Pa vs velocity at sea level, 15,000 ft, and 20,000 ft

%% Questions and table data

% 1) What's their maximum range, and at what velocity and altitude?

% Cl12Cd = (((1/3)*Cd0*pi*e*AR)^(1/4))/((4/3)*Cd0);
% R = (2*sqrt((2)/(rho8km*S))*(1/ct)*(Cl12Cd)*((Wi.^(1/2))-(Wf.^(1/2))))/1000;