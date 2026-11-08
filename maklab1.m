%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIGITAL PULSE RADAR PROJECT
% Part 1 : Radar Initialization & Signal Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%% ========================================================================
% RADAR PARAMETERS (According to Simulink Model)
% ========================================================================

c = 3e8;                      % Speed of light (m/s)

fc = 5e9;                     % Carrier frequency (Hz)

lambda = c/fc;

B = 3e6;                      % Sweep Bandwidth (Hz)

tau = 20e-6;                  % Pulse Width (s)

PRF = 6000;                   % Pulse Repetition Frequency (Hz)

PRI = 1/PRF;

fs = 12e6;                    % Sample Rate

N_pulses = 301;               % Same as Simulink

K = B/tau;                    % Chirp slope

%% ========================================================================
% TARGET PARAMETERS
% ========================================================================

targets = [...
    8000     20     1000;
    15000    80   100000];

N_targets = size(targets,1);

%% ========================================================================
% Derived Parameters
% ========================================================================

N_samples_PRI = round(PRI*fs);

N_samples_pulse = round(tau*fs);

range_resolution = c/(2*B);

max_range = c*PRI/2;

fprintf('\n');
fprintf('=====================================================\n');
fprintf('Digital Pulse Radar Simulation\n');
fprintf('=====================================================\n');

fprintf('Carrier Frequency      : %.2f GHz\n',fc/1e9);
fprintf('Bandwidth              : %.2f MHz\n',B/1e6);
fprintf('Sample Rate            : %.2f MHz\n',fs/1e6);
fprintf('Pulse Width            : %.2f us\n',tau*1e6);
fprintf('PRF                    : %.0f Hz\n',PRF);

fprintf('Range Resolution       : %.2f m\n',range_resolution);

fprintf('Maximum Unambiguous Range : %.2f km\n',max_range/1000);

fprintf('Number of Pulses       : %d\n',N_pulses);

fprintf('=====================================================\n');

%% ========================================================================
% LFM CHIRP GENERATION
% ========================================================================

t_pulse = (0:N_samples_pulse-1)/fs;

tx_pulse = exp(1j*pi*K*t_pulse.^2);

mf_kernel = conj(fliplr(tx_pulse));

%% ========================================================================
% Instantaneous Frequency
% ========================================================================

inst_freq = K*t_pulse;

%% ========================================================================
% Allocate Memory
% ========================================================================

rx_matrix = complex(zeros(N_pulses,N_samples_PRI));

%% ========================================================================
% Noise
% ========================================================================

SNR_dB = 20;

noise_sigma = 10^(-SNR_dB/20)/sqrt(2);

%% ========================================================================
% Generate Echoes
% ========================================================================

for pulse = 1:N_pulses

    slow_time = (pulse-1)*PRI;

    rx_signal = complex(zeros(1,N_samples_PRI));

    for target = 1:N_targets

        R0  = targets(target,1);

        vel = targets(target,2);

        rcs = sqrt(targets(target,3));

        current_range = R0 + vel*slow_time;

        delay = round((2*current_range/c)*fs);

        index = delay + (1:N_samples_pulse);

        if index(end)<=N_samples_PRI

            phase = exp(-1j*4*pi*fc*current_range/c);

            rx_signal(index)=rx_signal(index)+...
                rcs*tx_pulse.*phase;

        end

    end

    noise = noise_sigma*(randn(1,N_samples_PRI)+...
        1j*randn(1,N_samples_PRI));

    rx_matrix(pulse,:) = rx_signal + noise;

end

%% ========================================================================
% Axes
% ========================================================================

time_axis = t_pulse*1e6;

range_axis = (0:N_samples_PRI-1)*c/(2*fs)/1000;

%% ========================================================================
% COLORS
% ========================================================================

clrBlue = [0 0.4470 0.7410];

clrOrange = [0.8500 0.3250 0.0980];

clrGreen = [0.4660 0.6740 0.1880];

%% ========================================================================
% Figure 1
% TX LFM CHIRP
% ========================================================================

figure('Name','LFM Waveform',...
       'Position',[100 100 1100 420]);

subplot(1,2,1)

plot(time_axis,...
     real(tx_pulse),...
     'Color',clrBlue,...
     'LineWidth',1.5);

hold on

plot(time_axis,...
     imag(tx_pulse),...
     '--',...
     'Color',clrOrange,...
     'LineWidth',1.5);

xlabel('Time (\mus)')

ylabel('Amplitude')

title('Transmitted LFM Chirp')

legend('Real','Imaginary')

grid on

subplot(1,2,2)

plot(time_axis,...
     inst_freq/1e6,...
     'Color',clrGreen,...
     'LineWidth',1.6);

xlabel('Time (\mus)')

ylabel('Frequency (MHz)')

title('Instantaneous Frequency')

grid on

%% ========================================================================
% Figure 2
% Received Signal
% ========================================================================

figure('Name','Received Signal');

plot(range_axis,...
     real(rx_matrix(1,:)),...
     'Color',clrBlue,...
     'LineWidth',1.4);

xlabel('Range (km)')

ylabel('Amplitude')

title('Received Signal')

grid on

xlim([0 max(range_axis)])

disp('Part 1 Completed Successfully')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2 : Matched Filtering & Doppler Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================================================================
% MATCHED FILTER
% ========================================================================

mf_matrix = complex(zeros(size(rx_matrix)));

for pulse = 1:N_pulses

    mf_matrix(pulse,:) = conv(rx_matrix(pulse,:),...
                              mf_kernel,...
                              'same');

end

%% ========================================================================
% Matched Filter Magnitude
% ========================================================================

mf_mag = abs(mf_matrix);

%% ========================================================================
% Figure 3
% Matched Filter Output
% ========================================================================

figure('Name','Matched Filter Output');

mf_db = 20*log10(mf_mag(1,:) + eps);

mf_db = mf_db - max(mf_db);

plot(range_axis,...
     mf_db,...
     'Color',clrBlue,...
     'LineWidth',1.5);

grid on

xlabel('Range (km)')

ylabel('Magnitude (dB)')

title('Matched Filter Output')

ylim([-60 5])

xlim([0 max(range_axis)])

%% ========================================================================
% Average Range Profile
% ========================================================================

avg_profile = mean(mf_mag,1);

avg_db = 20*log10(avg_profile + eps);

avg_db = avg_db - max(avg_db);

%% ========================================================================
% Figure 4
% Range Profile
% ========================================================================

figure('Name','Range Profile');

plot(range_axis,...
     avg_db,...
     'Color',clrGreen,...
     'LineWidth',1.8);

grid on

xlabel('Range (km)')

ylabel('Magnitude (dB)')

title('Average Range Profile')

ylim([-60 5])

xlim([0 max(range_axis)])

%% ========================================================================
% Automatic Peak Detection
% ========================================================================

threshold = max(avg_profile)*0.30;

[pks,locs] = findpeaks(avg_profile,...
                       'MinPeakHeight',threshold,...
                       'MinPeakDistance',10);

estimated_ranges = range_axis(locs);

%% ========================================================================
% Figure 5
% Detected Targets
% ========================================================================

figure('Name','Detected Targets');

plot(range_axis,...
     avg_db,...
     'Color',clrBlue,...
     'LineWidth',1.5);

hold on

plot(estimated_ranges,...
     avg_db(locs),...
     'ro',...
     'MarkerFaceColor','r',...
     'MarkerSize',8);

for k = 1:length(estimated_ranges)

    text(estimated_ranges(k)+0.15,...
         avg_db(locs(k))+2,...
         sprintf('%.2f km',estimated_ranges(k)));

end

grid on

xlabel('Range (km)')

ylabel('Magnitude (dB)')

title('Detected Target Locations')

ylim([-60 5])

xlim([0 max(range_axis)])

%% ========================================================================
% Print Estimated Ranges
% ========================================================================

fprintf('\n');
fprintf('=========================================\n');
fprintf('Estimated Target Ranges\n');
fprintf('=========================================\n');

for k = 1:length(estimated_ranges)

    fprintf('Target %d : %.3f km\n',...
            k,...
            estimated_ranges(k));

end

fprintf('=========================================\n');

%% ========================================================================
% Find Closest Bins to True Targets
% ========================================================================

[~,bin1] = min(abs(range_axis*1000-targets(1,1)));

[~,bin2] = min(abs(range_axis*1000-targets(2,1)));

%% ========================================================================
% Slow-Time Signals
% ========================================================================

slow_signal1 = mf_matrix(:,bin1);

slow_signal2 = mf_matrix(:,bin2);

%% ========================================================================
% Figure 6
% Slow-Time Signal
% ========================================================================

figure('Name','Slow-Time Signal');

subplot(2,1,1)

plot(real(slow_signal1),...
     'Color',clrBlue,...
     'LineWidth',1.4);

hold on

plot(imag(slow_signal1),...
     '--',...
     'Color',clrOrange,...
     'LineWidth',1.4);

title('Target 1 Slow-Time Signal')

xlabel('Pulse Number')

ylabel('Amplitude')

legend('Real','Imaginary')

grid on

subplot(2,1,2)

plot(real(slow_signal2),...
     'Color',clrBlue,...
     'LineWidth',1.4);

hold on

plot(imag(slow_signal2),...
     '--',...
     'Color',clrOrange,...
     'LineWidth',1.4);

title('Target 2 Slow-Time Signal')

xlabel('Pulse Number')

ylabel('Amplitude')

legend('Real','Imaginary')

grid on

%% ========================================================================
% Doppler FFT
% ========================================================================

window = hamming(N_pulses);

doppler1 = fftshift(fft(slow_signal1.*window));

doppler2 = fftshift(fft(slow_signal2.*window));

doppler1_db = 20*log10(abs(doppler1)+eps);

doppler2_db = 20*log10(abs(doppler2)+eps);

doppler1_db = doppler1_db - max(doppler1_db);

doppler2_db = doppler2_db - max(doppler2_db);

velocity_axis = (-N_pulses/2:N_pulses/2-1)...
    *(PRF/N_pulses)*lambda/2;

%% ========================================================================
% Figure 7
% Doppler Spectrum
% ========================================================================

figure('Name','Doppler Spectrum');

plot(velocity_axis,...
     doppler1_db,...
     'Color',clrBlue,...
     'LineWidth',1.6);

hold on

plot(velocity_axis,...
     doppler2_db,...
     '--',...
     'Color',clrOrange,...
     'LineWidth',1.6);

xlabel('Velocity (m/s)')

ylabel('Magnitude (dB)')

title('Target Doppler Spectrum')

legend('Target 1','Target 2')

ylim([-60 5])

grid on

disp('Part 2 Completed Successfully');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3 : Range-Doppler Processing & Performance Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================================================================
% RANGE-DOPPLER MAP
% ========================================================================

RD = fftshift(fft(mf_matrix.*hamming(N_pulses),N_pulses,1),1);

RD_dB = 20*log10(abs(RD)+eps);

RD_dB = RD_dB - max(RD_dB(:));

%% ========================================================================
% Figure 8 : Range-Doppler Map
% ========================================================================

figure('Name','Range-Doppler Map',...
       'Position',[100 100 900 650]);

imagesc(range_axis,...
         velocity_axis,...
         RD_dB,...
         [-50 0]);

axis xy

xlabel('Range (km)')

ylabel('Velocity (m/s)')

title('Range-Doppler Map')

colormap(jet)

colorbar

hold on

plot(targets(1,1)/1000,...
     targets(1,2),...
     'w+',...
     'LineWidth',2,...
     'MarkerSize',12);

plot(targets(2,1)/1000,...
     targets(2,2),...
     'ws',...
     'MarkerFaceColor','w',...
     'MarkerSize',8);

legend('Target 1','Target 2');

%% ========================================================================
% Velocity Estimation
% ========================================================================

[~,idx1] = max(abs(doppler1));

[~,idx2] = max(abs(doppler2));

estimated_velocity1 = velocity_axis(idx1);

estimated_velocity2 = velocity_axis(idx2);

%% ========================================================================
% Print Estimated Velocities
% ========================================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('Estimated Target Velocities\n');
fprintf('============================================\n');

fprintf('Target 1 = %.2f m/s\n',estimated_velocity1);

fprintf('Target 2 = %.2f m/s\n',estimated_velocity2);

fprintf('============================================\n');

%% ========================================================================
% True vs Estimated Comparison
% ========================================================================

fprintf('\n');
fprintf('========================================================\n');
fprintf('TRUE VS ESTIMATED TARGET PARAMETERS\n');
fprintf('========================================================\n');

fprintf('Target 1\n');

fprintf(' True Range      : %.2f km\n',targets(1,1)/1000);

fprintf(' Estimated Range : %.2f km\n',estimated_ranges(1));

fprintf(' True Velocity   : %.2f m/s\n',targets(1,2));

fprintf(' Estimated Vel.  : %.2f m/s\n\n',estimated_velocity1);

fprintf('Target 2\n');

if length(estimated_ranges)>=2

fprintf(' True Range      : %.2f km\n',targets(2,1)/1000);

fprintf(' Estimated Range : %.2f km\n',estimated_ranges(2));

end

fprintf(' True Velocity   : %.2f m/s\n',targets(2,2));

fprintf(' Estimated Vel.  : %.2f m/s\n',estimated_velocity2);

fprintf('========================================================\n');

%% ========================================================================
% Velocity Ambiguity Demonstration
% ========================================================================

PRF_list = [2000 4000 6000 8000 10000];

figure('Name','Velocity Ambiguity',...
       'Position',[100 100 1500 350]);

for k=1:length(PRF_list)

    v_amb = lambda*PRF_list(k)/4;

    v = linspace(-v_amb,v_amb,2000);

    alias1 = mod(targets(1,2)+v_amb,...
                2*v_amb)-v_amb;

    alias2 = mod(targets(2,2)+v_amb,...
                2*v_amb)-v_amb;

    spectrum = exp(-(v-alias1).^2/2);

    spectrum = spectrum + ...
        0.8*exp(-(v-alias2).^2/2);

    subplot(1,length(PRF_list),k)

    plot(v,...
         spectrum/max(spectrum),...
         'LineWidth',1.5);

    hold on

    xline(alias1,'--b','LineWidth',1.2);

    xline(alias2,'--r','LineWidth',1.2);

    xlabel('Velocity')

    ylabel('Normalized')

    title(sprintf('PRF = %d Hz',PRF_list(k)))

    grid on

end

%% ========================================================================
% Radar Performance Summary
% ========================================================================

velocity_resolution = PRF/N_pulses*lambda/2;

maximum_velocity = lambda*PRF/4;

fprintf('\n');
fprintf('================ RADAR PERFORMANCE ==================\n');

fprintf('Carrier Frequency          : %.2f GHz\n',fc/1e9);

fprintf('Bandwidth                  : %.2f MHz\n',B/1e6);

fprintf('Pulse Width                : %.2f us\n',tau*1e6);

fprintf('Sampling Frequency         : %.2f MHz\n',fs/1e6);

fprintf('Number of Pulses           : %d\n',N_pulses);

fprintf('Range Resolution           : %.2f m\n',range_resolution);

fprintf('Maximum Range              : %.2f km\n',max_range/1000);

fprintf('Velocity Resolution        : %.2f m/s\n',...
        velocity_resolution);

fprintf('Maximum Unambiguous Speed  : %.2f m/s\n',...
        maximum_velocity);

fprintf('=====================================================\n');

%% ========================================================================
% Final Message
% ========================================================================

disp(' ');
disp('*****************************************************');
disp(' DIGITAL PULSE RADAR SIMULATION COMPLETED SUCCESSFULLY');
disp('*****************************************************');
disp(' ');
disp('Generated Outputs:');
disp(' ');
disp('1. LFM Chirp');
disp('2. Instantaneous Frequency');
disp('3. Received Signal');
disp('4. Matched Filter Output');
disp('5. Range Profile');
disp('6. Automatic Target Detection');
disp('7. Slow-Time Signal');
disp('8. Doppler Spectrum');
disp('9. Range-Doppler Map');
disp('10. Velocity Ambiguity');
disp('11. Estimated Target Range');
disp('12. Estimated Target Velocity');
disp('13. Performance Summary');
disp(' ');
disp('*****************************************************');
