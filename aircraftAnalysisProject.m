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
Tr20k = @(V) W./((W./(((.5)*rho20k*(V.^2)) * S))./(Cd0 + (((W./(((.5)*rho20k*(V.^2)) * S)).^2)/(pi*e*AR))));
Ta20k = Tr20k(V20kTest1);
TaSL = Ta20k*(rhoSL/rho20k);
TaSLTable = @(V) TaSL + (0.*V);
% SL
figure(1)
TrSL = @(V) W./((W./(((.5)*rhoSL*(V.^2)) * S))./(Cd0 + (((W./(((.5)*rhoSL*(V.^2)) * S)).^2)/(pi*e*AR))));
fplot(TrSL, [150 900]), title('Thrust Required and Thrust Available vs Velocity at Sea Level')
xlabel('Velocity [ft/s]'), ylabel('Thrust [lb_f]'), xlim([0 900]), ylim([0 6000]), grid on, hold on
fplot(@(V) TaSL+(0.*V), [0 900])
ThrustSLdataPoints = [722.78:.001:722.79; TaSLTable(722.78:.001:722.79); TrSL(722.78:.001:722.79);]';
plot(722.789, TrSL(722.789), 'r*')
legend('Thrust Required', 'Thrust Available', 'Maximum Velocity at Sea Level','Location', 'southeast'), hold off

% (2) Pr and Pa vs velocity at sea level, 15,000 ft, and 20,000 ft
% SL
figure(2)
PrSL = @(V) TrSL(V).*V;
PaSL = @(V) TaSL.*V;
fplot(@(V) PrSL(V).*0.001818, [60 775]), hold on, grid on
fplot(@(V) PaSL(V).*0.001818, [0 775]), xlabel('Velocity [ft/s]'), ylabel('Power [hp]')
title('Power Requred and Power Available vs Velocity at Sea Level')
legend('Power Required', 'Power Available', 'Location', 'southeast'), hold off
% 15k
figure(3)
Tr15k = @(V) W./((W./(((.5).*rho15k.*(V.^2)) .* S))./(Cd0 + (((W./(((.5).*rho15k.*(V.^2)) .* S)).^2)./(pi.*e.*AR))));
Ta15k = Ta20k.*(rho15k./rho20k);
Pr15k = @(V) Tr15k(V).*V;
Pa15k = @(V) Ta15k.*V;
fplot(@(V) Pr15k(V).*0.001818, [110 775]), hold on, grid on
fplot(@(V) Pa15k(V).*0.001818, [0 775]), xlabel('Velocity [ft/s]'), ylabel('Power [hp]')
title('Power Requred and Power Available vs Velocity at 15,000 [ft]')
legend('Power Required', 'Power Available', 'Location', 'southeast'), hold off
% 20k
figure(4)
Tr20k = @(V) W./((W./(((.5).*rho20k.*(V.^2)) .* S))./(Cd0 + (((W./(((.5).*rho20k.*(V.^2)) .* S)).^2)/(pi.*e.*AR))));
Pr20k = @(V) Tr20k(V).*V;
Pa20k = @(V) Ta20k.*V;
fplot(@(V) Pr20k(V).*0.001818, [125 775]), hold on, grid on
fplot(@(V) Pa20k(V).*0.001818, [0 775]), xlabel('Velocity [ft/s]'), ylabel('Power [hp]')
title('Power Requred and Power Available vs Velocity at 20,000 [ft]')
legend('Power Required', 'Power Available', 'Location', 'southeast'), hold off

% (3) Rate of climb vs velocity at sea level
% RC = (Pa - Pr)/W
figure(5)
RCSL = @(V) (PaSL(V)-PrSL(V))/W;
fplot(RCSL, [0 800]), ylim([0 100]), grid on, hold on
xlabel('Velocity [ft/s]'), ylabel('Rate of Climb [ft/s]')
title('Rate of Climb at Sea Level')
plot(4.4088*10^2, 0.6323*10^2, 'r*')
legend('Rate of Climb at Sea Level', 'Maximum Rate of Climb at Sea Level (63.23 [ft/s])')
RCSLdataPoints = [400:.01:450; RCSL(400:.01:450)]';
%% Questions and table data

% 1) What's their maximum range, and at what velocity and altitude?

% Cl12Cd = (((1/3)*Cd0*pi*e*AR)^(1/4))/((4/3)*Cd0);
% R = (2*sqrt((2)/(rho8km*S))*(1/ct)*(Cl12Cd)*((Wi.^(1/2))-(Wf.^(1/2))))/1000;