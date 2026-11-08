clc;
clear;
close all;

% SECTION 1: SYSTEM PARAMETERS
c=3e8;
fc=10e9;
lambda=c/fc;
B=2e6;
tau=20e-6;
PRF=5000;
PRI=1/PRF;
fs=10e6;
N_pulses=64;
CPI=N_pulses*PRI;
delta_r=c/(2*B);
delta_v=lambda/(2*CPI);
v_max=lambda*PRF/4;
K=B/tau;

% TARGETS
targets=[
    8e3   20   1.0;
    15e3  80   0.7
];
N_targets=size(targets,1);

% WAVEFORM
N_samples_PRI=round(PRI*fs);
N_samples_pulse=round(tau*fs);
dt=1/fs;
t_pulse=(0:N_samples_pulse-1)*dt;

tx_pulse=exp(1j*pi*K*t_pulse.^2);
mf_kernel=conj(fliplr(tx_pulse));

% RECEIVED SIGNAL GENERATION
SNR_dB=20;
noise_sigma=10^(-SNR_dB/20)/sqrt(2);

rx_matrix=zeros(N_pulses,N_samples_PRI);

for p=1:N_pulses

    t_slow=(p-1)*PRI;
    rx_signal=zeros(1,N_samples_PRI);

    for tgt=1:N_targets

        R0=targets(tgt,1);
        vel=targets(tgt,2);
        rcs=targets(tgt,3);

        R_p=R0+vel*t_slow;

        delay_s=floor((2*R_p/c)*fs);

        range_phase=exp(-1j*4*pi*fc*R0/c);
        doppler_phase=exp(1j*4*pi*vel*t_slow/lambda);

        idx_start=delay_s+1;
        idx_end=delay_s+N_samples_pulse;

        if idx_end<=N_samples_PRI

            rx_signal(idx_start:idx_end)=...
                rx_signal(idx_start:idx_end)+...
                rcs*tx_pulse.*...
                range_phase.*...
                doppler_phase;

        end
    end

    noise=noise_sigma*...
        (randn(1,N_samples_PRI)+...
         1j*randn(1,N_samples_PRI));

    rx_matrix(p,:)=rx_signal+noise;

end

% MATCHED FILTER
mf_matrix=zeros(N_pulses,N_samples_PRI);

for p=1:N_pulses

    mf_matrix(p,:)=conv(...
        rx_matrix(p,:),...
        mf_kernel,...
        'same');

end

% RANGE AXIS
range_axis=(0:N_samples_PRI-1)*(c/(2*fs));
range_axis_km=range_axis/1e3;

% MANUAL SHIFT
range_axis_shifted_km=range_axis_km-1.5;

% RANGE DOPPLER PROCESSING
ham_win=hamming(N_pulses);

rd_input=mf_matrix.*ham_win;

rd_map=fftshift(fft(rd_input,N_pulses,1),1);

fd_axis=(-N_pulses/2:N_pulses/2-1)*(PRF/N_pulses);

vel_axis=fd_axis*lambda/2;

rd_db=20*log10(abs(rd_map)+eps);
rd_db=rd_db-max(rd_db(:));

% COLORS
clr1=[0 0.4470 0.7410];
clr2=[0.8500 0.3250 0.0980];
clrG=[0.4660 0.6740 0.1880];

style=@(ax)set(ax,...
    'FontSize',10,...
    'LineWidth',1,...
    'Box','on');

%% FIGURE 1
fig1=figure('Name','Figure 1',...
    'Position',[50 100 1000 400],...
    'Color','w');

sgtitle('Figure 1 — LFM Waveform Analysis',...
    'FontSize',15,...
    'FontWeight','bold');

ax1=subplot(1,2,1);

plot(t_pulse*1e6,...
    real(tx_pulse),...
    'Color',clr1,...
    'LineWidth',1.4);

hold on;

plot(t_pulse*1e6,...
    imag(tx_pulse),...
    '--',...
    'Color',clr2,...
    'LineWidth',1.4);

xlabel('Time (\mus)');
ylabel('Amplitude');

title('TX LFM Chirp Waveform');

legend('Real','Imag');

grid on;

style(ax1);

ax2=subplot(1,2,2);

inst_freq=K*t_pulse/1e6;

plot(t_pulse*1e6,...
    inst_freq,...
    'Color',clrG,...
    'LineWidth',2);

xlabel('Time (\mus)');
ylabel('Relative Frequency (MHz)');

title('LFM Instantaneous Frequency');

grid on;

style(ax2);

%% FIGURE 2
fig2=figure('Name','Figure 2',...
    'Position',[150 150 700 450],...
    'Color','w');

ax3=axes;

plot(range_axis_km,...
    real(rx_matrix(1,:)),...
    'Color',clrG,...
    'LineWidth',1);

xlabel('Range (km)');
ylabel('Amplitude');

title('Received Signal — Pulse 1');

xlim([0 30]);

for i=1:N_targets

    xline(targets(i,1)/1e3,...
        '--r',...
        'LineWidth',1.2,...
        'Label',sprintf('T%d',i));

end

grid on;

style(ax3);

%% FIGURE 3
fig3=figure('Name','Figure 3',...
    'Position',[200 200 700 450],...
    'Color','w');

ax4=axes;

avg_profile=mean(abs(mf_matrix),1);

avg_db=20*log10(avg_profile+eps);

avg_db=avg_db-max(avg_db);

plot(range_axis_shifted_km,...
    avg_db,...
    'Color',clr1,...
    'LineWidth',1.5);

xlabel('Range (km)');
ylabel('Magnitude (dB)');

title('Range Profile');

xlim([0 30]);
ylim([-60 5]);

for i=1:N_targets

    xline(targets(i,1)/1e3,...
        '--r',...
        'LineWidth',1.2,...
        'Label',sprintf('T%d',i));

end

grid on;

style(ax4);

%% FIGURE 4
fig4=figure('Name','Figure 4',...
    'Position',[50 50 1450 850],...
    'Color','w');




% RANGE DOPPLER MAP
ax8=subplot(2,2,[3 4]);

imagesc(range_axis_km,...
    vel_axis,...
    rd_db,...
    [-50 0]);

axis xy;

colormap(ax8,jet);

cb=colorbar;

cb.Label.String='Magnitude (dB)';

xlabel('Range (km)');
ylabel('Velocity (m/s)');

title('Range-Doppler Map');

xlim([0 30]);
ylim([0 85])
hold on;

plot(targets(1,1)/1e3,...
    targets(1,2),...
    'w+',...
    'MarkerSize',12,...
    'LineWidth',2);

plot(targets(2,1)/1e3,...
    targets(2,2),...
    'ws',...
    'MarkerSize',10,...
    'LineWidth',2);

text(targets(1,1)/1e3+0.2,...
    targets(1,2)+3,...
    'T1',...
    'Color','w',...
    'FontWeight','bold');

text(targets(2,1)/1e3+0.2,...
    targets(2,2)+3,...
    'T2',...
    'Color','w',...
    'FontWeight','bold');

grid on;

style(ax8);

%% FIGURE 5
PRF_list=[2000 5000 10000];

fig5=figure('Name','Figure 5 - Effect of PRF',...
    'Position',[100 100 1500 450],...
    'Color','w');

sgtitle('Effect of PRF on Doppler Ambiguity',...
    'FontSize',18,...
    'FontWeight','bold');

for k_prf=1:length(PRF_list)

    PRF_temp=PRF_list(k_prf);

    lambda_temp=c/fc;

    v_amb=lambda_temp*PRF_temp/4;

    vel_axis_temp=linspace(-v_amb,v_amb,2000);

    sigma_v=1.2;

    v1=targets(1,2);
    v2=targets(2,2);

    alias_v1=mod(v1+v_amb,2*v_amb)-v_amb;
    alias_v2=mod(v2+v_amb,2*v_amb)-v_amb;

    spec1=exp(-(vel_axis_temp-alias_v1).^2/(2*sigma_v^2));

    spec2=0.7*exp(-(vel_axis_temp-alias_v2).^2/(2*sigma_v^2));

    spectrum_total=spec1+spec2;

    spectrum_total=spectrum_total/max(spectrum_total);

    ax=subplot(1,3,k_prf);

    plot(vel_axis_temp,...
        spectrum_total,...
        'Color',clr1,...
        'LineWidth',2);

    hold on;

    xline(alias_v1,...
        '--b',...
        'LineWidth',1.3,...
        'Label',sprintf('T1 Alias = %.1f',alias_v1));

    xline(alias_v2,...
        '--r',...
        'LineWidth',1.3,...
        'Label',sprintf('T2 Alias = %.1f',alias_v2));

    xlabel('Velocity (m/s)');
    ylabel('Amplitude');

    title(sprintf('PRF = %d Hz',PRF_temp));

    xlim([-v_amb v_amb]);
    ylim([0 1.05]);

    grid on;

    style(ax);

end