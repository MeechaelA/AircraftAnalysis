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
%%
W = 18000; % [lb] THIS IS FULLY FUELED WEIGHT
e;
AR;
S; % [ft^2]

% atmospheric conditions
rhoSL  = 2.3769*10^-3;  % [slugs/ft^3]
rho15k = 1.4962*10^-3; % [slugs/ft^3]
rho20k = 1.2673*10^-3; % [slugs/ft^3]

%% Plots
% (1) Tr and Ta vs velocity at sea level
% SL
CdSL = Cd0+(Cl^2/pi*e*AR)
ClSL = W/(0.5*rhoSL*V^2*S);
TrSL = @(V) W/(Cl/Cd);
% (2) Pr and Pa vs velocity at sea level, 15,000 ft, and 20,000 ft

%% Questions and table data

% 1) What's their maximum range, and at what velocity and altitude?

% Cl12Cd = (((1/3)*Cd0*pi*e*AR)^(1/4))/((4/3)*Cd0);
% R = (2*sqrt((2)/(rho8km*S))*(1/ct)*(Cl12Cd)*((Wi.^(1/2))-(Wf.^(1/2))))/1000;