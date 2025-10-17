% Progetto Ugello
clc, close all, clear all
%% DATI DI PROGETTO
P00 = .9; %atm
T00 = 300; % K
R = 287; % J/kg/K
k = 1.4; % Specific heat ratio for air
RHO00 = P00*101325/R/T00; % kg/m^3

Pchock_P00 = (2/(1+k))^(k/(k-1));
P3_P00_Progetto = .07/P00;
portata_progetto = 3.6; % kg/s

%% CALCOLI
% calcolo area ristretta e grafico rhoc sezione ristretta
rho_c_choking = sqrt(7*(P00*101325)*RHO00*((Pchock_P00)^(2/k)-(Pchock_P00)^((1+k)/k)));
Area_choking = portata_progetto/rho_c_choking; % m^2
raggio_chocking = 1000*sqrt(Area_choking/pi); % mm
x = linspace(0, 1, 1000); % rapporto pressione pressione gola/pressione totale
y = sqrt(7.*(P00*101325).*RHO00.*((x).^(2/k)-(x).^((1+k)/k))); % rho_c della gola
plot(x,y);
hold on
% plot rho_c della sezione di chocking
plot(x, (portata_progetto/Area_choking).*ones(length(x)));

% calcolo della area di valle
rho_c_valle = sqrt(7*(P00*101325)*RHO00*((P3_P00_Progetto)^(2/k)-(P3_P00_Progetto)^((1+k)/k)));
Area_valle = portata_progetto/rho_c_valle; % m^2
raggio_valle = 1000*sqrt(Area_valle/pi); % mm
hold on
% Plot rho_c della sezione di valle
plot(x, (portata_progetto/Area_valle).*ones(length(x)), 'r--');

% calcolo della area di monte 
T1_monte = T00 - .5; % K
c1_monte = sqrt(2*(k*R/(k-1))*(T00-T1_monte)); % m/s
P1_monte = P00 *(T00 / T1_monte)^(k/(1-k)); % atm
rho_1_monte = (101325*P1_monte)/R/T1_monte; % kg/m^3
Area_monte = portata_progetto/rho_1_monte/c1_monte; % m^2
raggio_monte = 1000*sqrt(Area_monte/pi); % mm
% Plot rho_c della sezione di monte
hold on
plot(x, (portata_progetto/Area_monte).*ones(length(x)), 'g--');









