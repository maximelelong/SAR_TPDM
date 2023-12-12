%close all 
clear variables; clc; close;

% Define a trajectory
TR.T = 1 ; %s
TR.N = 3; % nb periods
% TR.Q = [pi/2; 0];
TR.Q = [0; 0];
TR.C = [0 0 0 0; 0 0 0 0];
TR.S = [0 0 0 0; 0 0 0 0]; 
% TR.C = [0 1 0 0; 1 0 0 0];
% TR.S = [0 0 0 0; 0 1 0 2]; 
    
% Run script
option = 'only_acquisition_hardware' ; % 'full_robot' or  'only_acquisition_hardware'
[q, tau] = myrobot(TR,option);

%% Plot results

% Plot results
a = figure() ;
subplot(211); plot(q); xlabel('Time Sample'); ylabel('q');
title("Mesures de position")
legend('q_1', 'q_2')
subplot(212); plot(tau); xlabel('Time Sample'); ylabel('\tau');
legend('\tau_1', '\tau_2')
title("Mesures de couple")


%% FFT analysis
fe_q = 1e3;
fe_tau = 2.5e3;

filt_50hz = designfilt('bandstopiir', ...       % Response type
       'PassbandFrequency1',48, ...    % Frequency constraints
       'StopbandFrequency1',49.5, ...
       'StopbandFrequency2',50.5, ...
       'PassbandFrequency2',52, ...
       'PassbandRipple1',0.1, ...         % Magnitude constraints
       'StopbandAttenuation',15, ...
       'PassbandRipple2',0.1, ...
       'DesignMethod','ellip', ...      % Design method
       'SampleRate',fe_q) ;  

filtered_q1 = filtfilt(filt_50hz, q(:,1));


figure()
[freq_q1, fft_q1] = compute_fft(q(:,1), fe_q);
[~, fft_q1_filt] = compute_fft(filtered_q1, fe_q);
plot(freq_q1, fft_q1, LineWidth=2); hold on;
plot(freq_q1, fft_q1_filt, LineWidth=2);

