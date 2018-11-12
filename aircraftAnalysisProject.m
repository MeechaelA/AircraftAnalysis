% MAE 200 Aircraft Analysis Project
%
% Brady Burnsides
%
% Michael Angeles
%
% Giselle Padilla
%
% Kevin Lang
%
% Stephen Penney
hold off, close all,  clear,  clc, format long

%% General Calculations

W = 18000; % [lb] THIS IS FULLY FUELED WEIGHT
span = 50.87; % Tip2tip span [ft]
S = 345.12; % Planform Wing Area [ft^2]
AR = span^2/S; % Aspect Ratio [ft^2/ft^2]
e = 1.78*(1-(0.045*AR^(0.68)))-0.64;
% Atmospheric Conditions
rhoSL  = 2.377*10^-3;  % [slugs/ft^3]
rho15k = 1.4962*10^-3; % [slugs/ft^3]
rho20k = 1.2673*10^-3; % [slugs/ft^3]
Temp20k = 447.43; % [degR]
rho34k = 7.6696*10^-4; % [slugs/ft^3]
rho35k = 7.3820*10^-4; % [slugs/ft^3]

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

%% Thrust and Power Calculations
Tr20k = @(V) W./((W./(((.5)*rho20k*(V.^2)) * S))./(Cd0 + (((W./(((.5)*rho20k*(V.^2)) * S)).^2)/(pi*e*AR))));
Ta20k = Tr20k(V20kTest1);
Pr20k = @(V) Tr20k(V).*V;
Pa20k = @(V) Ta20k.*V;
RC20k = @(V)(Pa20k(V)-Pr20k(V))/W;

TaSL = Ta20k*(rhoSL/rho20k);
TrSL = @(V) W./((W./(((.5)*rhoSL*(V.^2)) * S))./(Cd0 + (((W./(((.5)*rhoSL*(V.^2)) * S)).^2)/(pi*e*AR))));
PrSL = @(V) TrSL(V).*V;
PaSL = @(V) TaSL.*V;
RCSL = @(V) (PaSL(V)-PrSL(V))/W;

Tr15k = @(V) W./((W./(((.5).*rho15k.*(V.^2)) .* S))./(Cd0 + (((W./(((.5).*rho15k.*(V.^2)) .* S)).^2)./(pi.*e.*AR))));
Ta15k = Ta20k.*(rho15k./rho20k);
Pr15k = @(V) Tr15k(V).*V;
Pa15k = @(V) Ta15k.*V;

Tr34k = @(V) W./((W./(((.5).*rho34k.*(V.^2)) .* S))./(Cd0 + (((W./(((.5).*rho34k.*(V.^2)) .* S)).^2)./(pi.*e.*AR))));
Ta34k = Ta20k.*(rho34k./rho20k);
Pr34k = @(V) Tr34k(V).*V;
Pa34k = @(V) Ta34k.*V;
RC34k = @(V)(Pa34k(V)-Pr34k(V))/W;

Tr35k = @(V) W./((W./(((.5).*rho35k.*(V.^2)) .* S))./(Cd0 + (((W./(((.5).*rho35k.*(V.^2)) .* S)).^2)./(pi.*e.*AR))));
Ta35k = Ta20k.*(rho35k./rho20k);
Pr35k = @(V) Tr35k(V).*V;
Pa35k = @(V) Ta35k.*V;
RC35k = @(V)(Pa35k(V)-Pr35k(V))/W;

%%

TSFC = (3000/11.06)/Ta20k;
ct = TSFC/3600;

%% Specified Plots

%%%%%%%%
% (1) Tr and Ta vs velocity at sea level

% Figure 1 - Tr and Ta [lbf] vs Velocity [ft/s] at Sea Level
figure(1), hold on, grid on
title('Thrust Required and Thrust Available [lb_f] vs Velocity [ft/s] at Sea Level')
xlabel('Velocity [ft/s]'), ylabel('Thrust [lb_f]'), xlim([0 900]), ylim([0 6000])
fplot(TrSL, [150 900]), plot([0 900], [TaSL TaSL])
legend('Thrust Required', 'Thrust Available','Location', 'southeast')
hold off

%%%%%%%%
% (2) Pr and Pa vs velocity at sea level, 15,000 ft, and 20,000 ft

% SL
% Figure 2 - Pr and Pa [hp] vs Velocity [ft/s] at Sea Level
figure(2), hold on, grid on
title('Power Requred and Power Available [hp] vs Velocity [ft/s] at Sea Level')
xlabel('Velocity [ft/s]'), ylabel('Power [hp]')
fplot(@(V) PrSL(V).*0.001818, [60 775]), fplot(@(V) PaSL(V).*0.001818, [0 775])
legend('Power Required', 'Power Available', 'Location', 'southeast')
hold off

% 15k
% Figure 3 - Pr and Pa [hp] vs Velocity at 15,000 [ft]
figure(3), hold on, grid on
title('Power Requred and Power Available [hp] vs Velocity [ft/s] at 15,000 [ft]')
xlabel('Velocity [ft/s]'), ylabel('Power [hp]')
fplot(@(V) Pr15k(V).*0.001818, [110 775]), fplot(@(V) Pa15k(V).*0.001818, [0 775])
legend('Power Required', 'Power Available', 'Location', 'southeast')
hold off

% 20k
% Figure 4 - Pr and Pa [hp] vs Velocity at 20,000 [ft]
figure(4), hold on, grid on
title('Power Requred and Power Available [hp] vs Velocity [ft/s] at 20,000 [ft]')
xlabel('Velocity [ft/s]'), ylabel('Power [hp]')
fplot(@(V) Pr20k(V).*0.001818, [125 775]), fplot(@(V) Pa20k(V).*0.001818, [0 775]),
legend('Power Required', 'Power Available', 'Location', 'southeast')
hold off

%%%%%%%%
% (3) Rate of climb vs velocity at sea level

% Figure 5 - R/C [ft/s] vs Velocity [ft/s] at SL
figure(5), hold on, grid on
title('Rate of Climb [ft/s] vs Velocity [ft/s] at Sea Level')
xlabel('Velocity [ft/s]'), ylabel('Rate of Climb [ft/s]')
fplot(RCSL, [0 800]), ylim([0 100])

negativeRCSL = @(V) -1*RCSL(V);
[velAtMaxRCSL, negMaxRCSL] = fminbnd(negativeRCSL, 0, 800); maxRCSL = -1*negMaxRCSL;
plot(velAtMaxRCSL, maxRCSL, 'r.', 'MarkerSize', 15)

legend('Rate of Climb at Sea Level', 'Maximum Rate of Climb at Sea Level (63.29 [ft/s])')
hold off
%% Questions

%%%%%%%%
% 1) What's their maximum range, and at what velocity and altitude?
% Maximum range occurs at the highest altitudes - so take abs ceiling of
% 35k
Cl12Cd35k = @(V) (((W./(((.5).*rho35k.*(V.^2)) .* S))).^(0.5)./(Cd0 + (((W./(((.5).*rho35k.*(V.^2)) .* S)).^2)./(pi.*e.*AR))));
V = 0:0.001:1000;
[maxCl12Cd35k, I] = max(Cl12Cd35k(V));
VelMaxRange = V(I);

figure(6)
fplot(Tr35k, [150 1000]), hold on, grid on
minTV35k = min(Tr35k(150:.01:1000)./(150:.01:1000));
fplot(@(V) minTV35k*V, [0 1000])
fplot(@(V) minTV35k*V, [0 1000])
plot(VelMaxRange,Tr35k(VelMaxRange),'g*')
xlabel('Velocity [ft/s]')
ylabel('Thrust Required')
% ANSWER
R = 2*sqrt(2/(rho35k*S))*(1/ct)*maxCl12Cd35k*(W^(0.5)-12000^(0.5));
%%%%%%%%

% 2)
W0 = 18000;
W1 = 13000;
rho = (2.3769*10^-3):-0.0000001:(7.382*10^-4);
VMaxE = 1:0.061020259:1001;
ClE = W./(rho.*VMaxE.^2.*S);
CdE = Cd0 + ((ClE.^2)./(pi*e*AR));
LoD = (0.5.*rho.*VMaxE.^2.*S.*ClE);
E = (1/ct).*(LoD).*log(W0./W1);
plot(VMaxE,Tr35k(VMaxE));
MaxEndurance = max(E);

%%%%%%%%
% 3) The helicarrier has an operating altitude of 15,000 ft. If these aircrafts engage the
% helicarrier, what�s the max duration could they operate?
% need to find max endurance at 15,000 ft

% Figure 7
figure(7), hold on, grid on
title('Thrust Required and Thrust Available [lb_f] vs Velocity [ft/s] at 15,000 [ft]')
xlabel('Velocity [ft/s]'), ylabel('Thrust [lb_f]'), xlim([0 900]), ylim([0 6000])
fplot(Tr15k, [150 900]), plot([0 900], [Ta15k Ta15k])
legend('Thrust Required', 'Thrust Available', 'Location', 'southeast')
hold off

%%%%%%%%
% 4) What�s their best absolute and service ceiling? And at what velocities?

% Need to obtain another max Rate of Climb value in order to linearly plot max 
% R/C against altitude and obtain absolute and service (100 ft/min)
% R/Cs
negativeRC20k = @(V) -1*RC20k(V);
[velAtMaxRC20k, negMaxRC20k] = fminbnd(negativeRC20k, 0, 800); maxRC20k = -1*negMaxRC20k;

% Figure 8 - Plots two previously determined maxR/Cs and then linearly
% interpolates the two to find the ceilings
figure(8), hold on, grid on
title('Maximum Rate of Climb [ft/s] vs Altitude [ft/10^3]')
ylabel('Altitude [ft/10^3]'), xlabel('Maximum Rate of Climb [ft/s]'), ylim([0 45])
plot([maxRCSL maxRC20k], [0 20], 'r.', 'MarkerSize', 15);

% Determine the linear function that passes through these two points
coeff = polyfit([maxRCSL maxRC20k], [0 20], 1);
m = coeff(1); absCeiling = coeff(2); maxRCvAlt = @(maxRC) m.*maxRC+absCeiling;
fplot(maxRCvAlt, [0 70]);

% Service ceiling represents a max RC of 100 ft/min or 1.667 ft/s
serviceRC = 1.667; % [ft/s]
serviceCeil = (m*1.66667)+absCeiling;

% Plot Ceilings
plot(serviceRC, serviceCeil, 'b.', 'MarkerSize', 15)
plot(0, absCeiling, 'm.', 'MarkerSize', 15)

legend('Previously calculated max R/Cs', 'Interpolated line from those points', 'Service Ceiling (R/C = 100 [ft/min])', 'Absolute Ceiling')
hold off


% Find the velocities associated with the absolute and service ceilings
% Absolute Ceiling (~35,000 [ft])
negativeRC35k = @(V) -1*RC35k(V);
[velAtMaxRC35k, negMaxRC35k] = fminbnd(negativeRC35k, 0, 800); maxRC35k = -1*negMaxRC35k;
% Service Ceiling (~34,000 [ft])
negativeRC34k = @(V) -1*RC34k(V);
[velAtMaxRC34k, negMaxRC34k] = fminbnd(negativeRC34k, 0, 800); maxRC34k = -1*negMaxRC34k;

%%%%%%%%
% 5) What�s their max rate of climb at sea level?
velAtMaxRCSL; maxRCSL;

%%
Velocity = [VelMaxRange VMaxE velAtMaxRC35k velAtMaxRC34k]; % Various Velocities
Q1 = [R Velocity(1,1) absCeiling];
% Q2 = [E Velocity(1,2) rho)];
% Q3 = [Duration]
Q4 = table([absCeiling Velocity(1,3); serviceCeil Velocity(1,4)]);      %Table for Question 4 Answers
%Q5 = table(RCSL);                                                       %Table for Q5 Answers


Titles = table({'Question1'; 'Question2'; 'Question3'; 'Question4'; 'Question5'}); %Creates Title Names for .xlsx file


%Q1,Q2,Q3,Q4,Q5 Title Writes
writetable(Titles(1,1),'MAE200AnalysisAnswers.xlsx','WriteVariableNames',0,'Sheet',1,'Range','A1')
writetable(Titles(2,1),'MAE200AnalysisAnswers.xlsx','WriteVariableNames',0,'Sheet',1,'Range','A6')
writetable(Titles(3,1),'MAE200AnalysisAnswers.xlsx','WriteVariableNames',0,'Sheet',1,'Range','A11')
writetable(Titles(4,1),'MAE200AnalysisAnswers.xlsx','WriteVariableNames',0,'Sheet',1,'Range','A16')
writetable(Titles(5,1),'MAE200AnalysisAnswers.xlsx','WriteVariableNames',0,'Sheet',1,'Range','A21')



%Q1,Q2,Q3,Q4,Q5 Answer Writes
% writetable(Q1,'MAE200AnalysisAnswers.xlsx','Sheet',1,'Range','B1')
% writetable(Q2,'MAE200AnalysisAnswers.xlsx','Sheet',1,'Range','B6')
% writetable(Q3,'MAE200AnalysisAnswers.xlsx','Sheet',1,'Range','B12')
% writetable(Q4,'MAE200AnalysisAnswers.xlsx','Sheet',1,'Range','B18')
% writetable(Q5,'MAE200AnalysisAnswers.xlsx','Sheet',1,'Range','B24')