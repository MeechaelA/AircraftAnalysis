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
% Tr and Ta vs velocity at sea level
% Pr and Pa vs velocity at sea level, 15,000 ft, and 20,000 ft
% Rate of climb vs velocity at sea level
% Table with info requested by Fury (range, endurance, ceiling)
% 1) What's their maximum range, and at what velocity and altitude?
% 2) What's their maximum endurance, and at what velocity and altitude?
% 3) The helicarrier has an operating altitude of 15,000 ft. If these aircrafts engage the
% helicarrier, what�s the max duration could they operate?
% 4) What's their best absolute and service ceiling? And at what velocities?
% 5) What's their max rate of climb at sea level?


Tr = W./(L./D);
Ta = q*S*CD0 + q*S*((CL^2)/(pi*e*AR))
% Pr and Pa vs velocity at sea level, 15,000 ft, and 20,000 ft
% SL

Pr = @(v) (W./((W ./ (.5.*rhoSL.*S.*(v.^2)))./ (Cd0 + (((W ./ (.5.*rhoSL.*S.*(v.^2))).^2)./(pi.*e.*AR))))) .* v;
Pa = @(v) Tsl.*v;

% 15,000 ft

% 20,000 ft