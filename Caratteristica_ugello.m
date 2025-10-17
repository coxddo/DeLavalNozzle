%% plot_pressione_portata.m
% P–Q da CSV (CFD) + curva teorica isentropica, choking e intersezioni con target.
% - Lettura colonne robusta (gestisce nomi diversi)
% - Etichette DP rialzate
% - Assi: X a 3 decimali (atm), Y a 2 decimali (kg/s)
% - Nessuna funzione di toolbox
% CALCOLO CURVA CARATTERISTICA UGELLO PORTATA - PRESSIONE

clear; clc; close all;

%% === INPUT ===
fname = 'DesignPoints.csv';
T = readtable(fname);

%% === Lettura colonne robusta ===
% Etichette DP
if ismember('DesignPoints', T.Properties.VariableNames)
    dp = string(T.DesignPoints);
else
    dp = "DP" + (1:height(T)).';
end

% Portata (prova vari nomi)
cand_q = ["minus-mout-plot-op","minus_mout_plot_op","portata","m_dot","flow"];
q = [];
for c = cand_q
    if ismember(c, T.Properties.VariableNames)
        q = T.(c); break;
    end
end
if isempty(q), error('Colonna portata non trovata. Nomi provati: %s', strjoin(cand_q, ', ')); end
q = q(:);

% Pressione in atm (prova vari nomi)
cand_p = ["p_out_atm","p_atm","P_atm","pout_atm"];
p_atm = [];
for c = cand_p
    if ismember(c, T.Properties.VariableNames)
        p_atm = T.(c); break;
    end
end
if isempty(p_atm), error('Colonna pressione [atm] non trovata. Nomi provati: %s', strjoin(cand_p, ', ')); end
p_atm = p_atm(:);

%% === Tabella e ordinamento ===
S = table(p_atm, q, dp, 'VariableNames', {'p_atm','q','dp'});
S = sortrows(S, 'p_atm', 'ascend');

%% === Grafico CFD P–Q ===
figure('Name','Pressione - Portata (CFD + Teorica)','NumberTitle','off');
plot(S.p_atm, S.q, 'o-','LineWidth',1.5,'MarkerSize',6); hold on; grid on; box on;
xlabel('Pressione [atm]'); ylabel('Portata [kg/s]');
title('Pressione - Portata (CFD + Teorica)');
ax = gca; ax.XAxis.TickLabelFormat = '%.3f'; ax.YAxis.TickLabelFormat = '%.2f';

% Etichette DP rialzate
spanX = max(S.p_atm) - min(S.p_atm);
spanY = max(S.q)     - min(S.q);
dx = 0.00 * (spanX + eps);
dy = 0.03 * (spanY + eps);
for k = 1:height(S)
    text(S.p_atm(k)+dx, S.q(k)+dy, S.dp(k), ...
        'FontSize',9,'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','Interpreter','none');
end

%% === Riepilogo console ===
fprintf('Riepilogo (ordinati per pressione):\n');
for k = 1:height(S)
    fprintf('  %s: p = %.3f atm, q = %.2f kg/s\n', S.dp(k), S.p_atm(k), S.q(k));
end

%% === Teorica isentropica (ugello) ===
P0_atm = 0.9;                % pressione totale monte [atm]
T0     = 300;                % K
k      = 1.4;                % gamma
R      = 287;                % J/(kg*K)
Area_valle   = pi*(0.1)^2;   % m^2 (r=0.1 m)
RHO10 = P0_atm*101325/R/T0;  % kg/m^3

P0_pa  = P0_atm*101325;                         % Pa
P_pa   = linspace(0, P0_pa, 600);      % Pa (evita 0 esatto)
Pr     = P_pa./P0_pa;                           % P/P0
P_iso  = P_pa/101325;                           % atm (asse X teorica)

% calcolo rho_c teorica
rho_c1 = Pr.^(2/k)-Pr.^((1+k)/k);
rho_c2 = 2*k/(k-1).*rho_c1;
rho_c = sqrt(rho_c2.*P0_pa*RHO10);
m_teo  = Area_valle .*rho_c;

% Assicura vettori colonna (evita errori di concatenazione)
P_iso = P_iso(:);
m_teo = m_teo(:);

% Choking
Area_chocking = pi*(73.38/1000)^2; % m^2
Pr_crit   = (2/(k+1))^(k/(k-1));
Pcrit_atm = Pr_crit * P0_atm;
rho_c1ch = Pr_crit.^(2/k)-Pr_crit.^((1+k)/k);
rho_c2ch = 2*k/(k-1).*rho_c1ch;
rho_c_ch = sqrt(rho_c2ch.*P0_pa*RHO10);
m_chocking = Area_chocking*rho_c_ch;

% Plot teorica + linea P*
plot(P_iso, m_teo, 'r--', 'LineWidth', 1.5);
if Pcrit_atm>=min(P_iso) && Pcrit_atm<=max(P_iso)
    xline(Pcrit_atm, ':', 'P^* (critica)', 'LabelVerticalAlignment','middle');
end

legend('CFD','Teorica isentropica','Location','best');

%% === Intersezioni con target ===
target = m_chocking;                           % kg/s

y  = m_teo - target; s = sign(y);
ix = find(s(1:end-1).*s(2:end) <= 0 & ~isnan(s(1:end-1)) & ~isnan(s(2:end)));

P_int_atm = [];
for i = ix(:).'
    x1 = P_iso(i); x2 = P_iso(i+1);
    y1 = m_teo(i); y2 = m_teo(i+1);
    if y2~=y1
        Pk = x1 + (target - y1) * (x2 - x1) / (y2 - y1);
        P_int_atm(end+1) = Pk; %#ok<AGROW>
    end
end

% Linea orizzontale target + markers
yline(target, ':', sprintf('\\it m = %.1f kg/s', target), 'LabelVerticalAlignment','bottom');
for i = 1:numel(P_int_atm)
    plot(P_int_atm(i), target, 'ko', 'MarkerFaceColor','k');
end

% Stampa intersezioni
if isempty(P_int_atm)
    fprintf('Nessuna intersezione: target %.2f kg/s fuori range [%.2f, %.2f]\n', ...
        target, min(m_teo), max(m_teo));
else
    fprintf('Intersezioni per m_dot = %.2f kg/s:\n', target);
    for i = 1:numel(P_int_atm)
        fprintf('  P = %.3f atm\n', P_int_atm(i));
    end
end

%% === Padding verticale (senza errori di dimensione) ===
y_all  = [S.q(:); m_teo(:)];
ymin_a = min(y_all);
ymax_a = max(y_all);
yspan  = max(ymax_a - ymin_a, eps);
ylim([ymin_a - 0.05*yspan, ymax_a + 0.10*yspan]);
