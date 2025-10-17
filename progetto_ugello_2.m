% Progetto Ugello
clc, close all, clear all

%% DATI DI PROGETTO
P00 = .9;            % atm (pressione totale a monte)
T00 = 300;           % K   (temperatura totale a monte)
R   = 287;           % J/kg/K (gas costante specifica aria)
k   = 1.4;           % (-) rapporto calori specifici per aria
RHO00 = P00*101325/R/T00;     % kg/m^3 densità totale corrispondente

Pchock_P00       = (2/(1+k))^(k/(k-1));  % rapporto P* / P0 a choking isentropico
P3_P00_Progetto  = .07/P00;              % rapporto P3 / P0 di progetto
portata_progetto = 3.6;                  % kg/s (portata di progetto)

%% CALCOLI
% calcolo area ristretta e grafico rhoc sezione ristretta
rho_c_choking  = sqrt(7*(P00*101325)*RHO00*((Pchock_P00)^(2/k)-(Pchock_P00)^((1+k)/k)));
Area_choking   = portata_progetto/rho_c_choking;       % m^2
raggio_chocking = 1000*sqrt(Area_choking/pi);          % mm

% dominio per x = P_gola / P0 e funzione rho_c(x) per la gola
x = linspace(0, 1, 1000);                              % (-) rapporto pressione gola/pressione totale
y = sqrt(7.*(P00*101325).*RHO00.*((x).^(2/k)-(x).^((1+k)/k))); % rho_c alla gola

% grafico base
plot(x, y); hold on

% linea orizzontale corrispondente a m_dot / A_gola
plot(x, (portata_progetto/Area_choking).*ones(length(x)));

% calcolo della area di valle
rho_c_valle = sqrt(7*(P00*101325)*RHO00*((P3_P00_Progetto)^(2/k)-(P3_P00_Progetto)^((1+k)/k)));
Area_valle  = portata_progetto/rho_c_valle;            % m^2
raggio_valle = 1000*sqrt(Area_valle/pi);               % mm

% linea orizzontale corrispondente a m_dot / A_valle
hold on
plot(x, (portata_progetto/Area_valle).*ones(length(x)), 'r--');

% calcolo della area di monte (condizione leggermente accelerata)
T1_monte   = T00 - .5;                                 % K
c1_monte   = sqrt(2*(k*R/(k-1))*(T00-T1_monte));       % m/s (da energia)
P1_monte   = P00 *(T00 / T1_monte)^(k/(1-k));          % atm (da isentropica)
rho_1_monte = (101325*P1_monte)/R/T1_monte;            % kg/m^3
Area_monte  = portata_progetto/rho_1_monte/c1_monte;   % m^2
raggio_monte = 1000*sqrt(Area_monte/pi);               % mm

% linea orizzontale corrispondente a m_dot / A_monte
hold on
plot(x, (portata_progetto/Area_monte).*ones(length(x)), 'g--');

%% MIGLIORIE DI LETTURA (solo etichette/estetica, nessun impatto sui calcoli)
xlabel('x = P_{gola} / P_0 (-)');
ylabel('\rho_c  [kg/(m^2 s)]');
title('Progetto ugello — \rho_c(x) e vincoli di portata/area');
legend('\rho_c(x) alla gola', '\dot{m}/A_{gola}', '\dot{m}/A_{valle}', '\dot{m}/A_{monte}', 'Location','best');
grid on; box on;
xlim([0 1]);
