
% Code écrit par Solal BITTOUN, Lilian DELORY et Maxime LELONG - MASTER SAR
% Dans le cadre du TP/DM du cours d'Estimation et Identification.
% Code permettant d'étudier le bruit blanc associé aux signaux.

%% Cleanup
clear variables; clc; close all;

%% Generate measures

% Define a trajectory
TR.T = 40 ; %s
TR.N = 4; % nb periods
% TR.Q = [pi/2; 0];
% TR.C = [0 1 0 0; 1 0 0 0];
% TR.S = [0 0 0 0; 0 1 0 2]; 

TR.Q = [0; pi];
TR.C = [0 0 0 0; 0 0 0 0];
TR.S = [0 0 0 0; 0 0 0 0]; 
    
% Run script
option = 'full_robot' ; % 'full_robot' or  'only_acquisition_hardware'
[q, tau] = myrobot(TR,option);

params.fe_q = 1e3;
params.fe_tau = 2.5e3;
signals.q = q;
signals.tau = tau;

%% White noise analysis before processing
% With no input to the system => only measurement noise (assumed gaussian)
fprintf("------------------- BRUITS DE MESURES COUPLE -----------------------------\n")
tau1 = signals.tau(:,1);
mu_1 = mean(tau1);
V_1 = var(tau1);
fprintf("Bruit de mesures de tau 1 :\n");
fprintf("   -> N(%f, %f)\n", mu_1, V_1);
tau2 = signals.tau(:,2);
mu_2 = mean(tau2);
V_2 = var(tau2);
fprintf("Bruit de mesures de tau 2 :\n");
fprintf("   -> N(%f, %f)\n", mu_2, V_2);

fprintf("------------------- BRUITS DE MESURES POSITION -----------------------------\n")

q1 = signals.q(:,1);
mu_1q = mean(q1);
V_1q = var(q1);
fprintf("Bruit de mesures de q_1 :\n");
fprintf("   -> N(%f, %f)\n", mu_1q, V_1q);

q2 = signals.q(:,2);
mu_2q = mean(q2);
V_2q = var(q2);
fprintf("Bruit de mesures de q_1 :\n");
fprintf("   -> N(%f, %f)\n", mu_2q, V_2q);

%% Plot results
figure()
subplot(211); plot(signals.q); xlabel('Time Sample'); ylabel('q');
legend('q_1', 'q_2'); title("Mesures des positions") ;
subplot(212); plot(signals.tau); xlabel('Time Sample'); ylabel('\tau');
legend('\tau_1', '\tau_2') ; title("Mesures des couples ");

%% Filter 50Hz
filt_50hz_1k = designfilt('bandstopiir', ...       % Response type
       'PassbandFrequency1',48, ...    % Frequency constraints
       'StopbandFrequency1',49.5, ...
       'StopbandFrequency2',50.5, ...
       'PassbandFrequency2',52, ...
       'PassbandRipple1',0.1, ...         % Magnitude constraints
       'StopbandAttenuation',15, ...
       'PassbandRipple2',0.1, ...
       'DesignMethod','ellip', ...      % Design method
       'SampleRate',params.fe_q) ;              % Sample rate
filt_50hz_2_5K = designfilt('bandstopiir', ...       % Response type
       'PassbandFrequency1',48, ...    % Frequency constraints
       'StopbandFrequency1',49.5, ...
       'StopbandFrequency2',50.5, ...
       'PassbandFrequency2',52, ...
       'PassbandRipple1',0.1, ...         % Magnitude constraints
       'StopbandAttenuation',15, ...
       'PassbandRipple2',0.1, ...
       'DesignMethod','ellip', ...      % Design method
       'SampleRate',params.fe_tau) ;              % Sample rate

signals.q_no_50 = filtfilt(filt_50hz_1k, q);
signals.tau_no_50 = filtfilt(filt_50hz_2_5K, tau);

%% Plot ffts Positions
figure;
[freq_q1, fft_q1] = compute_fft(signals.q(:,1), params.fe_q);
[~, fft_q1_filt] = compute_fft(signals.q_no_50(:,1), params.fe_q);

[freq_q2, fft_q2] = compute_fft(signals.q(:,2), params.fe_q);
[~, fft_q2_filt] = compute_fft(signals.q_no_50(:,2), params.fe_q);

subplot(2,1,1)
plot_fft_orig = plot(freq_q1, fft_q1, LineWidth=2); hold on;
plot_fft_orig.Color(4) = 0.30;
plot(freq_q1, fft_q1_filt,LineWidth=2);
xlabel('Fréquence (Hz)')
ylabel('Amplitude (dB)')
legend("q_1" ,"q_1 filtré");

subplot(2,1,2)
plot_fft_orig = plot(freq_q2, fft_q2, LineWidth=2); hold on;
plot_fft_orig.Color(4) = 0.30;
plot(freq_q2, fft_q2_filt,LineWidth=2);
xlabel('Fréquence (Hz)')
ylabel('Amplitude (dB)')
legend("q_2" ,"q_2 filtré");
sgtitle("Comparaison des FFT des signaux de positions avec et sans bruit de secteur");

%% Plot ffts Couples

figure;
[freq_tau1, fft_tau1] = compute_fft(tau(:,1), params.fe_tau);
[~, fft_tau1_filt] = compute_fft(signals.tau_no_50(:,1), params.fe_tau);

[freq_tau2, fft_tau2] = compute_fft(tau(:,2), params.fe_tau);
[~, fft_tau2_filt] = compute_fft(signals.tau_no_50(:,2), params.fe_tau);

subplot(2,1,1)
plot_fft_orig = plot(freq_tau1, fft_tau1, LineWidth=2); hold on;
plot_fft_orig.Color(4) = 0.30;
plot(freq_tau1, fft_tau1_filt,LineWidth=2);
axis([0, 1250, -300, 0])
xlabel('Fréquence (Hz)')
ylabel('Amplitude (dB)')
legend("\tau_1" ,"\tau_1 filtré");

subplot(2,1,2)
plot_fft_orig = plot(freq_tau2, fft_tau2, LineWidth=2); hold on;
plot_fft_orig.Color(4) = 0.30;
plot(freq_tau2, fft_tau2_filt,LineWidth=2);
axis([0, 1250, -300, 0])
xlabel('Fréquence (Hz)')
ylabel('Amplitude (dB)')
legend("\tau_2" ,"\tau_2 filtré");
sgtitle("Comparaison des FFT des signaux de couples avec et sans bruit de secteur");

%% Resample
signals.tau_1k = resample(signals.tau_no_50,params.fe_q,params.fe_tau);
signals.tau_1k = signals.tau_1k(1:end-5,:);

%% White noise analysis
% With no input to the system => only measurement noise (assumed gaussian)
fprintf("------------------- BRUITS DE MESURES COUPLES APRES FILTRAGE 50 HZ -----------------------------\n")
tau1 = signals.tau_no_50(:,1);
mu_1 = mean(tau1);
V_1 = var(tau1);
fprintf("Bruit de mesures de tau 1 :\n");
fprintf("   -> N(%f, %f)\n", mu_1, V_1);
tau2 = signals.tau_no_50(:,2);
mu_2 = mean(tau2);
V_2 = var(tau2);
fprintf("Bruit de mesures de tau 2 :\n");
fprintf("   -> N(%f, %f)\n", mu_2, V_2);

fprintf("------------ BRUITS DE MESURES COUPLES APRES FILTRAGE 50 HZ et REECHANTILLONNAGE -------\n")
tau1 = signals.tau_1k(:,1);
mu_1 = mean(tau1);
V_1 = var(tau1);
fprintf("Bruit de mesures de tau 1 :\n");
fprintf("   -> N(%f, %f)\n", mu_1, V_1);
tau2 = signals.tau_1k(:,2);
mu_2 = mean(tau2);
V_2 = var(tau2);
fprintf("Bruit de mesures de tau 2 :\n");
fprintf("   -> N(%f, %f)\n", mu_2, V_2);


