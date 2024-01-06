
% Code écrit par Solal BITTOUN, Lilian DELORY et Maxime LELONG - MASTER SAR
% Dans le cadre du TP/DM du cours d'Estimation et Identification.
% Code qui permet d'étudier le décalage entre les deux interfaces.

%% Generate measures
clear 
clc;
% Define a trajectory
TR.T = 2 ; %s
TR.N = 5; % nb periods
TR.Q = [pi/2; 0];
TR.C = [1 0 0 0; 1 0 0 0];
TR.S = [0 1 0 2; 0 1 0 2]; 

%TR.Q = [0; 0];
%TR.C = [0 0 0 0; 0 0 0 0];
%TR.S = [0 0 0 0; 0 0 0 0]; 
    
% Run script
option = 'only_acquisition_hardware' ; % 'full_robot' or  'only_acquisition_hardware'
[q, tau] = myrobot(TR,option);

params.fe_q = 1e3;
params.fe_tau = 2.5e3;

signals.q = q;
signals.tau_1k= resample(tau,params.fe_q,params.fe_tau);

signals.tau_1k_reduit(:, 1) = signals.tau_1k(1:length(signals.q(:,1)), 1);

signals.tau_1k_reduit(:, 2) = signals.tau_1k(1:length(signals.q(:,2)), 2);

signals.tau_1k = signals.tau_1k_reduit;

%% Plot results
figure()
subplot(211); plot(q); xlabel('Time Sample'); ylabel('q');
legend('q_1', 'q_2'); title("Mesures avec interface position") ;
subplot(212); plot(signals.tau_1k); xlabel('Time Sample'); ylabel('\tau');
legend('\tau_1', '\tau_2') ; title("Mesures avec interface couple ");


%% FFT
figure;
[freq_q1, fft_q1] = compute_fft(signals.q(:,1), params.fe_q);
[freq_tau2, fft_tau2] = compute_fft(signals.tau_1k(:,2), params.fe_q);


subplot(2,1,1)
plot_fft_q1 = plot(freq_q1, fft_q1, LineWidth=1);
plot_fft_orig.Color(4) = 0.30;
title("Signal 1");
xlabel('Fréquence (Hz)')
ylabel('Amplitude (dB)')

subplot(2,1,2)
plot_fft_ = plot(freq_tau2, fft_tau2, LineWidth=1);
plot_fft_orig.Color(4) = 0.30;
title("Signal 2");
xlabel('Fréquence (Hz)')
ylabel('Amplitude (dB)')
sgtitle("FFT des signaux utilisés pour l'étude du décalage");


%%
% Spécifications du filtre passe bas
order_aa = 8; % Ordre du filtre

% Création du filtre Butterworth passe-bas
[b_aa, a_aa] = butter(order_aa, 20/(1000/2), 'low');

% Application du filtre passe bas.
signals.filtered_tau = filtfilt(b_aa, a_aa, signals.tau_1k);
signals.filtered_q = filtfilt(b_aa, a_aa, signals.q);

% %% On selectionne les signaux en tronquant les extrémités
% % Sélectionner des échantillons au centre du signal de tau
% signals.tau_1k_reduit1(:, 1) = signals.filtered_tau(100:length(signals.filtered_tau(:, 1))-100, 1);
% signals.tau_1k_reduit1(:, 2) = signals.filtered_tau(100:length(signals.filtered_tau(:, 2))-100, 2);
% signals.filtered_tau = signals.tau_1k_reduit1;
% 
% 
% % Sélectionner des échantillons au centre du signal de q 
% signals.q_reduit(:, 1) = signals.filtered_q(100:length(signals.filtered_q(:, 1))-100, 1);
% signals.q_reduit(:, 2) = signals.filtered_q(100:length(signals.filtered_q(:, 2))-100, 2);
% signals.filtered_q = signals.q_reduit;

%% Plot results after filtering
figure()
subplot(211); plot(signals.filtered_q); xlabel('Time Sample'); ylabel('q');
legend('q_1', 'q_2'); title("Mesures avec interface position") ;
subplot(212); plot(signals.filtered_tau); xlabel('Time Sample'); ylabel('\tau');
legend('\tau_1', '\tau_2') ; title("Mesures avec interface couple ");

%% Plot results of same signal on same figure
figure()
subplot(211); plot(signals.filtered_q(:,1)); xlabel('Time Sample'); ylabel('q');
hold on;
subplot(211); plot(signals.filtered_tau(:,1));
legend('q_1', 'tau_1'); title("Signal 1") ;

subplot(212); plot(signals.filtered_q(:,2)); xlabel('Time Sample'); ylabel('q');
hold on;
subplot(212); plot(signals.filtered_tau(:,2));
legend('q_2', '\tau_2') ; title("Signal 2");

%% Decalage 
figure()
[correlation, lag] = xcorr(signals.filtered_tau(:,1),signals.filtered_q(:,1));
plot(correlation);
title('Intercorrelation');
[~, idx] = max(correlation);
time_shift_samples = lag(idx);
disp("Décalage temporelle calculée (signal 1) - (en nombre d'échantillons) : " + time_shift_samples);

figure()
[correlation, lag] = xcorr(signals.filtered_tau(:,2),signals.filtered_q(:,2));
plot(correlation);
title('Intercorrelation');
[~, idx] = max(correlation);
time_shift_samples = lag(idx);
disp("Décalage temporelle calculée (signal 2) - (en nombre d'échantillons) : " + time_shift_samples);

%% Avec finddelay 
d1 = finddelay(signals.filtered_q(:,1),signals.filtered_tau(:,1));
d2 = finddelay(signals.filtered_q(:,2),signals.filtered_tau(:,2));

disp("Décalage temporelle calculée (signal 1) - (en nombre d'échantillons) : " + d1);
disp("Décalage temporelle calculée (signal 2) - (en nombre d'échantillons) : " + d2);

%% t - alignment

signals.tau_aligned = signals.filtered_tau(5:end,:);
signals.q_aligned = signals.filtered_q(1:end-4,:);

%% Plot results on same figure
figure()
subplot(211); plot(signals.q_aligned(:,1)); xlabel('Time Sample'); ylabel('q');
hold on;
subplot(211); plot(signals.tau_aligned(:,1));
legend('q_1', 'tau_1'); title("Signal 1") ;

subplot(212); plot(signals.q_aligned(:,2)); xlabel('Time Sample'); ylabel('q');
hold on;
subplot(212); plot(signals.tau_aligned(:,2));
legend('q_2', '\tau_2') ; title("Signal 2");