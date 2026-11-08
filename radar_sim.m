%% ============================================================
% Pulse Radar Post Processing
% MATLAB R2024b
% =============================================================
clc;
close all;

%% Radar Parameters
c = 3e8;                  % Speed of light (m/s)
fc = 5.6e9;               % Carrier frequency
lambda = c/fc;
PRF = 2000;               % Pulse Repetition Frequency (Hz)
fs = 10e6;                % Sampling frequency (CHANGE if different)

%% Extract Data

rx_raw = squeeze(out_rx_raw(:,1,:));      % 2000 x 301
mf = squeeze(out_mf(:,1,:));              % 2000 x 301
mf_mag = squeeze(out_mf_mag(:,1,:));      % 2000 x 301

Ns = size(rx_raw,1);
Np = size(rx_raw,2);

%% =============================================================
% Figure 1 : Received Signal
%% =============================================================

figure;
plot(abs(rx_raw(:,1)),'LineWidth',1.4);
grid on;
xlabel('Sample');
ylabel('Amplitude');
title('Received Signal');

%% =============================================================
% Figure 2 : Matched Filter Output
%% =============================================================

figure;
plot(abs(mf(:,1)),'LineWidth',1.4);
grid on;
xlabel('Sample');
ylabel('Amplitude');
title('Matched Filter Output');

%% =============================================================
% Figure 3 : Range Profile
%% =============================================================

range_axis = (0:Ns-1)*c/(2*fs);

figure;
plot(range_axis,abs(mf(:,1)),'LineWidth',1.6);
grid on;
xlabel('Range (m)');
ylabel('Magnitude');
title('Range Profile');

%% =============================================================
% Figure 4 : Automatic Target Detection
%% =============================================================

threshold = 0.30*max(abs(mf(:,1)));

[pks,locs] = findpeaks(abs(mf(:,1)),...
    'MinPeakHeight',threshold,...
    'MinPeakDistance',20);

ranges = range_axis(locs);

figure;
plot(range_axis,abs(mf(:,1)),'b','LineWidth',1.5);
hold on;

plot(ranges,pks,'ro',...
    'MarkerFaceColor','r',...
    'MarkerSize',8);

grid on;

xlabel('Range (m)');
ylabel('Magnitude');
title('Detected Targets');

%% Display Estimated Range

disp('-------------------------------------');
disp('Estimated Target Ranges');

for k=1:length(ranges)

    fprintf('Target %d : %.2f meters\n',k,ranges(k));

end

disp('-------------------------------------');

%% =============================================================
% Figure 5 : Range Doppler Map
%% =============================================================

RD = fftshift(fft(mf,Np,2),2);

fd = (-Np/2:Np/2-1)*(PRF/Np);

velocity = fd*lambda/2;

figure;

imagesc(velocity,range_axis,...
    20*log10(abs(RD)+eps));

axis xy;

xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Range-Doppler Map');

colorbar;

colormap jet;

%% =============================================================
% Figure 6 : Doppler Spectrum
%% =============================================================

midBin = locs(1);

doppler = abs(RD(midBin,:));

figure;

plot(velocity,doppler,...
    'LineWidth',1.5);

grid on;

xlabel('Velocity (m/s)');
ylabel('Magnitude');
title('Target Doppler Spectrum');

%% =============================================================
% Estimate Velocity
%% =============================================================

[~,idx] = max(doppler);

estimated_velocity = velocity(idx);

fprintf('\nEstimated Velocity = %.2f m/s\n',estimated_velocity);

%% =============================================================
% Summary
%% =============================================================

fprintf('\n====================================\n');
fprintf('Number of Pulses      : %d\n',Np);
fprintf('Samples per Pulse     : %d\n',Ns);
fprintf('Range Resolution      : %.2f m\n',c/(2*fs));
fprintf('Estimated Velocity    : %.2f m/s\n',estimated_velocity);
fprintf('====================================\n');