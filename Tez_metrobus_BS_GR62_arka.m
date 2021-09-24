clear all; close all; clc;
tic%zaman başladı

filename = 'C:\Users\Aleyna\Documents\MATLAB\BS_GR62.xlsx';
[filepath,name,ext] = fileparts(filename)

%Excel'den datayı varible olarak matlab'e aktarma
%filename = 'BS_YRS1.xlsx'; %path görebilmek için filename üzerinden aynı dosyayı uzantısı ile çektik
data =xlsread(filename);
time =xlsread(filename,'A:A');
acc_signal_x =xlsread(filename,'B:B');
acc_signal_y =xlsread(filename,'C:C');
acc_signal_z =xlsread(filename,'D:D');
acc_signal_total =xlsread(filename,'E:E');

l_raw=length(time) %length of the data: zaman vektöründen çekilen veri sayısı
l_intp=2*l_raw-1; %(2n-1) will give new interpolated length
dt_raw= time(2)-time(1);

 %interpolation of acceleration vectors:  
for i = 1:l_raw
l_array(i)=i   ; %length's array form for pchip interpolation
end
 tq = 1:0.5:l_raw; %see help pchip
 xqp = pchip(l_array,acc_signal_x,tq);
 acc_intp_x=transpose (xqp);
 yqp = pchip(l_array,acc_signal_y,tq);
 acc_intp_y=transpose (yqp);
 zqp = pchip(l_array,acc_signal_z,tq);
 acc_intp_z=transpose (zqp);
 tqp = pchip(l_array,acc_signal_total,tq);
 acc_intp_total=transpose (tqp);
  %interpolation of time vectors:  
time_qp = pchip(l_array,time,tq);
time_inp=transpose (time_qp);
%% Sampling Period and Frequency of Interpolated Data

%Sampling Period of Interpolated Time Series
dt=time_inp(2)-time_inp(1);%0.005
L=length(time_inp);%size of interpolated data
Fs=1/dt; %198.86Hz
f_nyquist=Fs/2; %According to Nyquist maximum observable frequency of the system
f_practice=Fs/5; %According to Low Freq Ptactical Approah maximum observable frequency of the system
 %% 8th Order Butterworth HP & LP Filtering
 
 [b,a]=butter(8,0.03,'high'); % designing butterworth highpass filter: high-pass cutoff frequency 0.03 * 150 = 4.5 hz
 c_off=floor(f_nyquist)/150; % Floor for rounding down the Estimation: 150 Hz is for cutoff frequency ratioing for Butterworth filter package in matlab
 [bb,aa]=butter(8,c_off,'low');% low-pass cutoff frequency ie: (99)/150 * 150 = 99 hz
  
 
 %X-Direction Filtering & Integrating from Acc to Vel to Displacement
 % acceleration x-direction signal highpass filtered at 1Hz:
 accf_x=filter(b,a,acc_intp_x);
% acc-X signal low pass filtered at ~100Hz:
 accfl_x=filter(bb,aa,acc_intp_x);
% acc signal low&high pass filtered between: 4.5 to ~100Hz:
 acc_hp_lp_filtered_x=filter(bb,aa,accf_x); %(4.5Hz-100Hz).
%Integrating acceleration signal to gain displacement:
 veli_x=1000*cumsum(acc_hp_lp_filtered_x)*dt; % filtered acceleration integrated & conv to mm
 velf_x=filter(b,a,veli_x); %integrated velocity highpass filtered
 dispi_x=cumsum(velf_x)*dt; % filtered velocity integrated to displacement
 disp_x=filter(b,a,dispi_x); %integrated displacement is highpass filtered
 %%
 %%FOR AXISES, FIGURES ARE CODED AS SUCH: AXIS-X: figure (11"1"),figure (11"2"), AXIS-Y: FIG(2"1") FIG(2"2") AXIS-Z: FIG(3"1") FIG(3"2"),AXIS-TOTAL: FIG(4"1") FIG(4"2")..etc
 %% 
 %Plotting Time varying- Acceleration, Velocity and Displacement Signals:
 %Time varying- Acceleration in m/s^2 and 'g':
 figure (11)
 plot(time_inp, accfl_x,'m')%magenta color
 title('Filtered Acceleration-X Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Acceleration-X (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 11. Time-varying Filtered Input Acceleration-X Signal(m/s^2) \n ')
 
 %Time Varying Filtered & Integrated Velocity-X:
 figure (12)
 plot(time_inp, accfl_x/(9.81),'m')%where g=9.81 m/s^2, magenta color
 title('Filtered Acceleration-X Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Acceleration-X (g)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 12. Time-varying Filtered Input Acceleration-X Signal(g) \n ')

 %Time varying- Velocity
 figure (13)
 plot(time_inp, velf_x,'g')%green color
 title('Filtered Velocity-X','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Velocity -X(mm/sec)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 13. Time-varying Filtered Input Velocity-X(mm/sec) \n ')

 
  %Time varying- Displacement
 figure (14)
 plot(time_inp, disp_x)%original blue color
 title('Displacement-X Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Displacement-X (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 14. Time-varying Input Displacement-X (mm) \n ')

 
 %% FFTs of acceleration and displacement signals
 
% FFT of Low-pass filtered Acceleration signal with f_cutoff= 100Hz.
 NFFT = 2^nextpow2(L); % Next power of 2 from length of y: it rounds up to the number such that will be as a power of 2
 accflfft_x = fft(accfl_x,NFFT)/L; % FFT of Low-pass filtered Acceleration signal
 f = Fs/2*linspace(0,1,NFFT/2); % FFT freq axis, sampling frequency is Fs=99
 
 
 %fft of Acceleration plot
 figure (15)
 plot(f,2*abs(accflfft_x(1:NFFT/2)),'m')%magenta color. For Single-Sided fft: 2*(abs value).
 xlim([0 100])
 title_head = sprintf('Single-Sided Low-pass Filtered \n Acceleration-X Amplitude Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|Y(f)|, (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 15. Single-Sided low-pass Filtered Acceleration-X Amplitude Spectrum (m/s^2, Hz)\n ')
 
 %the same fft of acc plot in 'dB'
 figure (16)
 plot(f,20*log10(2*abs(accflfft_x(1:NFFT/2))),'r')%red color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 title_head = sprintf('Single-Sided Low-pass Filtered \n Acceleration-X Amplitude Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('20*log10(|Y(f)|), (dB)','fontsize',11,'Fontname','Timesnewroman');% dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 16. Single-Sided low-pass Filtered Acceleration-x Amplitude Spectrum (dB, Hz) \n ')
 
 
 
 %FFT of Position signal:
 dispfft_x = fft(disp_x,NFFT)/L; % in (mm)s. FFT of displacement signal
 %fft of position plot:
 figure (17)
 plot(f,2*abs(dispfft_x(1:NFFT/2)))%For Single-Sided fft: 2*(abs value)
 xlim([0 100])
 title_head = sprintf('Single-Sided \n Displacement Amplitude-X Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|U(f)|, (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 17. Single-Sided Displacement Amplitude-X Spectrum (mm, Hz) \n ')

%ZOOMED fft of position plot:
 figure (18)
 plot(f,2*abs(dispfft_x(1:NFFT/2)))%For Single-Sided fft: 2*(abs value)
 axis([3.28 11 0 0.005])
 title_head = sprintf('ZOOMED Single-Sided \n Displacement Amplitude-X Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('ZOOMED Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|U(f)|, (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 18. Zoomed Single-Sided Displacement Amplitude-X Spectrum (mm, Hz) \n ')
 %% 
 %the same FFT of position plot in 'dB'
 figure (19)
 plot(f,20*log10(2*abs(dispfft_x(1:NFFT/2))), 'c')%cyan color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 title('Single-Sided Displacement Amplitude-X Spectrum','fontsize',11,'Fontname','Timesnewroman')
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('20*log10(|U(f)|), (dB)','fontsize',11,'Fontname','Timesnewroman');% dB=20log10(A1/A2) & dB=10log10(Power1/Power2),
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 19. Single-Sided Displacement Amplitude-X Spectrum (dB, Hz) \n ')
 
 %% PSDs of acceleration FFT and displacement FFTs:
 
 
%PSD of the acceleration FFT:
 xdft_acc_x = fft(accfl_x); %FFT of the acc is selected for PSD evaluation
 xdft_acc_x = xdft_acc_x(1:(L+1)/2);%where L=length(accfl) not length(acc_signal) not length(accfft_x)
 psdx_acc_x = (1/(Fs*L)) * abs(xdft_acc_x).^2; %Power is directly proportional with the square of the absolute value of the position amplitude
 psdx_acc_x(2:end-1) = 2*psdx_acc_x(2:end-1); %dB/Hz
 freq_acc_x = 0:Fs/L:Fs/2;
 % PSD of FFT of Acceleration Plot
 figure (110)
 plot(freq_acc_x,psdx_acc_x,'m')%magenta color
 xlim([0 100])
 grid on
 title_head = sprintf('Periodogram Using FFT of the \n Low-pass Filtered Acceleration-X Signal');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency ( ((m/s^2)^2)/Hz )','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 fprintf(' \n Figure 110. Periodogram Using FFT of the low-pass filtered Acceleration Signal (((m/s^2)^2)/Hz, Hz) \n ')
 
 
 %the same PSD of Acceleration Plot in 'dB'
 figure (111)
 plot(freq_acc_x,10*log10(psdx_acc_x),'r')%magenta color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 %xlim([0 100])
 grid on
 title_head = sprintf('PSD Using FFT of the \n Low-pass Filtered Acceleration-X Signal');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 fprintf(' \n Figure 111. PSD Using FFT of the low-pass filtered Acceleration-X Signal (dB/Hz, Hz) \n ')

 
 % PSD of position signal calculation for position:
 xdft_disp_x = fft(disp_x/1000); %in meters , FFT of the acc is selected for PSD evaluation
 xdft_disp_x = xdft_disp_x(1:(L+1)/2);%where L=length(accfl) not length(acc_signal) not length(accfft)
 psdx_disp_x = (1/(Fs*L)) * abs(xdft_disp_x).^2; %Power is directly proportional with the square of the absolute value of the position amplitude
 psdx_disp_x(2:end-1) = 2*psdx_disp_x(2:end-1); %dB/Hz
 freq_disp_x = 0:Fs/L:Fs/2;
 % PSD of FFT of the position signal Plot
 figure (112)
 plot(freq_disp_x,(psdx_disp_x))%original blue color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 ylim([0 0.3e-7])
 title('Periodogram Using FFT of the Displacement-X','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency ((m^2)/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 112. Periodogram Using FFT of Displacement-X ((m^2)/Hz, Hz) \n ')

 
 %the same PSD of Acceleration Plot in 'dB'
 figure (114)
 plot(freq_disp_x,10*log10(psdx_disp_x),'c')%magenta color
 %xlim([0 100])
 title('Periodogram Using FFT of the Displacement-X','fontsize',11,'Fontname','Timesnewroman')
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 114. PSD Using FFT of Displacement-X (dB/Hz, Hz) \n ')
%% CONSLUSIVE WELCH PSD OVER ACCELERATION
%VWelch Power Spectral Estimator for Acceleration and displacement signals:
% Create a Welch spectral estimator for acc and disp.
 h_welch= spectrum.welch;
% Welch PSD of Raw and Filtered Acceleration:
 %One-sided PSD of welch of the raw Acceleration-X signal:
 Hpsd_raw_acc_x=psd(h_welch,acc_intp_x,'Fs',Fs);
 %One-sided PSD of welch of the low-pass filtered Acceleration-X signal:
 Hpsd_filtered_acc_x=psd(h_welch,accflfft_x,'Fs',Fs);
 % PSD of Welch spectral of the raw Acceleration-X Signal Plot:
 figure (115)
 plot(Hpsd_raw_acc_x)
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Raw Acceleration-X Signal \n Welch Power/Freq, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Raw Metrobus Acceleration-X Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 115. Welch Power Spectrum of Raw Metrobus Acceleration-X Signal(dB/Hz, Hz) \n ')

  figure (116)
 plot(Hpsd_filtered_acc_x)
 xlim([0 150])
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Metrobus Low-pass Filtered Acc-X Signal \n Welch Power/Frequency, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Metrobus Low-pass Filtered Acc-X Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 116. Welch Power Spectrum of Metrobus Filtered Acceleration-X Signal (dB/Hz, Hz) \n ')

%% FINAL CONCLUSION OVER WELCH PSD OF DISPLACEMENT

%One-sided PSD of welch of Displacement:
 Hpsd_disp_x=psd(h_welch,disp_x,'Fs',Fs);% Calculate and plot the one-sided PSD.
 % PSD of Welch spectral of Displacement Plot:
 figure (117)
 plot(Hpsd_disp_x)
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Metrobus Displacement-X Welch Power/Frequency, \n (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('One-sided Welch PSD of Metrobus Displacement-X','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 117. Welch Power Spectrum of Metrobus Displacement-X(dB/Hz, Hz) \n ')
 
 %% FOR Y AXIS:
 
 
 %Y-Direction Filtering & Integrating from Acc to Vel to Displacement
 % acceleration Y-direction signal highpass filtered at 1Hz:
 accf_y=filter(b,a,acc_intp_y);
% acc signal low pass filtered at ~100Hz:
 accfl_y=filter(bb,aa,acc_intp_y);
% acc signal low&high pass filtered between: 4.5 to ~100Hz:
 acc_hp_lp_filtered_y=filter(bb,aa,accf_y); %(4.5Hz-100Hz).
%Integrating acceleration signal to gain displacement:
 veli_y=1000*cumsum(acc_hp_lp_filtered_y)*dt; % filtered acceleration integrated & conv to mm
 velf_y=filter(b,a,veli_y); %integrated velocity highpass filtered
 dispi_y=cumsum(velf_y)*dt; % filtered velocity integrated to displacement
 disp_y=filter(b,a,dispi_y); %integrated displacement is highpass filtered
 
 %Plotting Time varying- Acceleration, Velocity and Displacement Signals:
 %Time varying- Acceleration in m/s^2 and 'g':
 figure (21)
 plot(time_inp, accfl_y,'m')%magenta color
 title('Filtered Acceleration-Y Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Acceleration-Y (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 21. Time-varying Filtered Input Acceleration-Y Signal(m/s^2) \n ')
 
 %Time Varying Filtered & Integrated Velocity-Y:
 figure (22)
 plot(time_inp, accfl_y/(9.81),'m')%where g=9.81 m/s^2, magenta color
 title('Filtered Acceleration-Y Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Acceleration-Y (g)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 22. Time-varying Filtered Input Acceleration-Y Signal(g) \n ')

 %Time varying- Velocity
 figure (23)
 plot(time_inp, velf_y,'g')%green color
 title('Filtered Velocity-Y','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Velocity -Y(mm/sec)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 23. Time-varying Filtered Input Velocity-Y(mm/sec) \n ')

 
  %Time varying- Displacement
 figure (24)
 plot(time_inp, disp_y)%original blue color
 title('Displacement-Y Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Displacement-Y (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 24. Time-varying Input Displacement-Y (mm) \n ')

 
 %% FFTs of acceleration and displacement signals
 
% FFT of Low-pass filtered Acceleration signal with f_cutoff= 100Hz.
 NFFT = 2^nextpow2(L); % Next power of 2 from length of y: it rounds up to the number such that will be as a power of 2
 accflfft_y = fft(accfl_y,NFFT)/L; % FFT of Low-pass filtered Acceleration signal
 f = Fs/2*linspace(0,1,NFFT/2); % FFT freq axis, sampling frequency is Fs=99
 
 %fft of Acceleration plot
 figure (25)
 plot(f,2*abs(accflfft_y(1:NFFT/2)),'m')%magenta color. For Single-Sided fft: 2*(abs value).
 xlim([0 100])
 title_head = sprintf('Single-Sided Low-pass Filtered \n Acceleration-Y Amplitude Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|Y(f)|, (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 25. Single-Sided low-pass Filtered Acceleration-Y Amplitude Spectrum (m/s^2, Hz)\n ')
 
 %the same fft of acc plot in 'dB'
 figure (26)
 plot(f,20*log10(2*abs(accflfft_y(1:NFFT/2))),'r')%red color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 title_head = sprintf('Single-Sided Low-pass Filtered \n Acceleration-Y Amplitude Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('20*log10(|Y(f)|), (dB)','fontsize',11,'Fontname','Timesnewroman');% dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 26. Single-Sided low-pass Filtered Acceleration-Y Amplitude Spectrum (dB, Hz) \n ')
 
 
 
 %FFT of Position signal:
 dispfft_y = fft(disp_y,NFFT)/L; % in (mm)s. FFT of displacement signal
 %fft of position plot:
 figure (27)
 plot(f,2*abs(dispfft_y(1:NFFT/2)))%For Single-Sided fft: 2*(abs value)
 xlim([0 100])
 title_head = sprintf('Single-Sided \n Displacement Amplitude-Y Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|U(f)|, (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 27. Single-Sided Displacement Amplitude-Y Spectrum (mm, Hz) \n ')

%ZOOMED fft of position plot:
 figure (28)
 plot(f,2*abs(dispfft_y(1:NFFT/2)))%For Single-Sided fft: 2*(abs value)
 axis([3.28 11 0 0.005])
 title_head = sprintf('ZOOMED Single-Sided \n Displacement Amplitude-Y Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('ZOOMED Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|U(f)|, (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 28. Zoomed Single-Sided Displacement Amplitude-Y Spectrum (mm, Hz) \n ')
 %% 
 %the same FFT of position plot in 'dB'
 figure (29)
 plot(f,20*log10(2*abs(dispfft_y(1:NFFT/2))), 'c')%cyan color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 title('Single-Sided Displacement Amplitude-Y Spectrum','fontsize',11,'Fontname','Timesnewroman')
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('20*log10(|U(f)|), (dB)','fontsize',11,'Fontname','Timesnewroman');% dB=20log10(A1/A2) & dB=10log10(Power1/Power2),
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 29. Single-Sided Displacement Amplitude-Y Spectrum (dB, Hz) \n ')
 
 %% PSDs of acceleration FFT and displacement FFTs:
 
 
%PSD of the acceleration FFT:
 xdft_acc_y = fft(accfl_y); %FFT of the acc is selected for PSD evaluation
 xdft_acc_y = xdft_acc_y(1:(L+1)/2);%where L=length(accfl) not length(acc_signal) not length(accfft_y)
 psdx_acc_y = (1/(Fs*L)) * abs(xdft_acc_y).^2; %Power is directly proportional with the square of the absolute value of the position amplitude
 psdx_acc_y(2:end-1) = 2*psdx_acc_y(2:end-1); %dB/Hz
 freq_acc_y = 0:Fs/L:Fs/2;
 % PSD of FFT of Acceleration Plot
 figure (210)
 plot(freq_acc_y,psdx_acc_y,'m')%magenta color
 xlim([0 100])
 grid on
 title_head = sprintf('Periodogram Using FFT of the \n Low-pass Filtered Acceleration-Y Signal');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency ( ((m/s^2)^2)/Hz )','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 fprintf(' \n Figure 210. Periodogram Using FFT of the low-pass filtered Acceleration Signal (((m/s^2)^2)/Hz, Hz) \n ')
 
 
 %the same PSD of Acceleration Plot in 'dB'
 figure (211)
 plot(freq_acc_y,10*log10(psdx_acc_y),'r')%magenta color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 %xlim([0 100])
 grid on
 title_head = sprintf('PSD Using FFT of the \n Low-pass Filtered Acceleration-Y Signal');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 fprintf(' \n Figure 211. PSD Using FFT of the low-pass filtered Acceleration-Y Signal (dB/Hz, Hz) \n ')

 
 % PSD of position signal calculation for position:
 xdft_disp_y = fft(disp_y/1000); %in meters , FFT of the acc is selected for PSD evaluation
 xdft_disp_y = xdft_disp_y(1:(L+1)/2);%where L=length(accfl) not length(acc_signal) not length(accfft)
 psdx_disp_y = (1/(Fs*L)) * abs(xdft_disp_y).^2; %Power is directly proportional with the square of the absolute value of the position amplitude
 psdx_disp_y(2:end-1) = 2*psdx_disp_y(2:end-1); %dB/Hz
 freq_disp_y = 0:Fs/L:Fs/2;
 % PSD of FFT of the position signal Plot
 figure (212)
 plot(freq_disp_y,(psdx_disp_y))%original blue color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 ylim([0 0.3e-7])
 title('Periodogram Using FFT of the Displacement-Y','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency ((m^2)/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 212. Periodogram Using FFT of Displacement-Y ((m^2)/Hz, Hz) \n ')

 
 %the same PSD of Acceleration Plot in 'dB'
 figure (214)
 plot(freq_disp_y,10*log10(psdx_disp_y),'c')%magenta color
 %xlim([0 100])
 title('Periodogram Using FFT of the Displacement-Y','fontsize',11,'Fontname','Timesnewroman')
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 214. PSD Using FFT of Displacement-Y (dB/Hz, Hz) \n ')
%% CONSLUSIVE WELCH PSD OVER ACCELERATION
%VWelch Power Spectral Estimator for Acceleration and displacement signals:
% Create a Welch spectral estimator for acc and disp.
 h_welch= spectrum.welch;
% Welch PSD of Raw and Filtered Acceleration:
 %One-sided PSD of welch of the raw Acceleration-Y signal:
 Hpsd_raw_acc_y=psd(h_welch,acc_intp_y,'Fs',Fs);
 %One-sided PSD of welch of the low-pass filtered Acceleration-Y signal:
 Hpsd_filtered_acc_y=psd(h_welch,accflfft_y,'Fs',Fs);
 % PSD of Welch spectral of the raw Acceleration-Y Signal Plot:
 figure (215)
 plot(Hpsd_raw_acc_y)
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Raw Acceleration-Y Signal \n Welch Power/Freq, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Raw Metrobus Acceleration-Y Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 215. Welch Power Spectrum of Raw Metrobus Acceleration-Y Signal(dB/Hz, Hz) \n ')

  figure (216)
 plot(Hpsd_filtered_acc_y)
 xlim([0 150])
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Metrobus Low-pass Filtered Acc-Y Signal \n Welch Power/Frequency, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Metrobus Low-pass Filtered Acc-Y Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 216. Welch Power Spectrum of Metrobus Filtered Acceleration-Y Signal (dB/Hz, Hz) \n ')

%% FINAL CONCLUSION OVER WELCH PSD OF DISPLACEMENT

%One-sided PSD of welch of Displacement:
 Hpsd_disp_y=psd(h_welch,disp_y,'Fs',Fs);% Calculate and plot the one-sided PSD.
 % PSD of Welch spectral of Displacement Plot:
 figure (217)
 plot(Hpsd_disp_y)
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Metrobus Displacement-Y Welch Power/Frequency, \n (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('One-sided Welch PSD of Metrobus Displacement-Y','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 217. Welch Power Spectrum of Metrobus Displacement-Y(dB/Hz, Hz) \n ')
 
 
 %% FOR Z-AXIS
 
 
 %Z-Direction Filtering & Integrating from Acc to Vel to Displacement
 % acceleration Z-direction signal highpass filtered at 1Hz:
 accf_z=filter(b,a,acc_intp_z);
% acc signal low pass filtered at ~100Hz:
 accfl_z=filter(bb,aa,acc_intp_z);
% acc signal low&high pass filtered between: 4.5 to ~100Hz:
 acc_hp_lp_filtered_z=filter(bb,aa,accf_z); %(4.5Hz-100Hz).
%Integrating acceleration signal to gain displacement:
 veli_z=1000*cumsum(acc_hp_lp_filtered_z)*dt; % filtered acceleration integrated & conv to mm
 velf_z=filter(b,a,veli_z); %integrated velocity highpass filtered
 dispi_z=cumsum(velf_z)*dt; % filtered velocity integrated to displacement
 disp_z=filter(b,a,dispi_z); %integrated displacement is highpass filtered
 
 %Plotting Time varying- Acceleration, Velocity and Displacement Signals:
 %Time varying- Acceleration in m/s^2 and 'g':
 figure (31)
 plot(time_inp, accfl_z,'m')%magenta color
 title('Filtered Acceleration- Z Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Acceleration- Z (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 31. Time-varying Filtered Input Acceleration- Z Signal(m/s^2) \n ')
 
 %Time Varying Filtered & Integrated Velocity- Z:
 figure (32)
 plot(time_inp, accfl_z/(9.81),'m')%where g=9.81 m/s^2, magenta color
 title('Filtered Acceleration- Z Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Acceleration- Z (g)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 32. Time-varying Filtered Input Acceleration- Z Signal(g) \n ')

 %Time varying- Velocity
 figure (33)
 plot(time_inp, velf_z,'g')%green color
 title('Filtered Velocity- Z','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Velocity - Z(mm/sec)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 33. Time-varying Filtered Input Velocity- Z(mm/sec) \n ')

 
  %Time varying- Displacement
 figure (34)
 plot(time_inp, disp_z)%original blue color
 title('Displacement- Z Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Displacement- Z (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 34. Time-varying Input Displacement- Z (mm) \n ')

 
 %% FFTs of acceleration and displacement signals
 
% FFT of Low-pass filtered Acceleration signal with f_cutoff= 100Hz.
 NFFT = 2^nextpow2(L); % Next power of 2 from length of y: it rounds up to the number such that will be as a power of 2
 accflfft_z = fft(accfl_z,NFFT)/L; % FFT of Low-pass filtered Acceleration signal
 f = Fs/2*linspace(0,1,NFFT/2); % FFT freq axis, sampling frequency is Fs=99
 
 %fft of Acceleration plot
 figure (35)
 plot(f,2*abs(accflfft_z(1:NFFT/2)),'m')%magenta color. For Single-Sided fft: 2*(abs value).
 xlim([0 100])
 title_head = sprintf('Single-Sided Low-pass Filtered \n Acceleration- Z Amplitude Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|Y(f)|, (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 35. Single-Sided low-pass Filtered Acceleration- Z Amplitude Spectrum (m/s^2, Hz)\n ')
 
 %the same fft of acc plot in 'dB'
 figure (36)
 plot(f,20*log10(2*abs(accflfft_z(1:NFFT/2))),'r')%red color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 title_head = sprintf('Single-Sided Low-pass Filtered \n Acceleration- Z Amplitude Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('20*log10(|Y(f)|), (dB)','fontsize',11,'Fontname','Timesnewroman');% dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 36. Single-Sided low-pass Filtered Acceleration- Z Amplitude Spectrum (dB, Hz) \n ')
 
 
 
 %FFT of Position signal:
 dispfft_z = fft(disp_z,NFFT)/L; % in (mm)s. FFT of displacement signal
 %fft of position plot:
 figure (37)
 plot(f,2*abs(dispfft_z(1:NFFT/2)))%For Single-Sided fft: 2*(abs value)
 xlim([0 100])
 title_head = sprintf('Single-Sided \n Displacement Amplitude- Z Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|U(f)|, (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 37. Single-Sided Displacement Amplitude- Z Spectrum (mm, Hz) \n ')

%ZOOMED fft of position plot:
 figure (38)
 plot(f,2*abs(dispfft_z(1:NFFT/2)))%For Single-Sided fft: 2*(abs value)
 axis([3.28 11 0 0.005])
 title_head = sprintf('ZOOMED Single-Sided \n Displacement Amplitude- Z Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('ZOOMED Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|U(f)|, (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 38. Zoomed Single-Sided Displacement Amplitude- Z Spectrum (mm, Hz) \n ')
 %% 
 %the same FFT of position plot in 'dB'
 figure (39)
 plot(f,20*log10(2*abs(dispfft_z(1:NFFT/2))), 'c')%cyan color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 title('Single-Sided Displacement Amplitude- Z Spectrum','fontsize',11,'Fontname','Timesnewroman')
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('20*log10(|U(f)|), (dB)','fontsize',11,'Fontname','Timesnewroman');% dB=20log10(A1/A2) & dB=10log10(Power1/Power2),
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 39. Single-Sided Displacement Amplitude- Z Spectrum (dB, Hz) \n ')
 
 %% PSDs of acceleration FFT and displacement FFTs:
 
 
%PSD of the acceleration FFT:
 xdft_acc_z = fft(accfl_z); %FFT of the acc is selected for PSD evaluation
 xdft_acc_z = xdft_acc_z(1:(L+1)/2);%where L=length(accfl) not length(acc_signal) not length(accfft_z)
 psdx_acc_z = (1/(Fs*L)) * abs(xdft_acc_z).^2; %Power is directly proportional with the square of the absolute value of the position amplitude
 psdx_acc_z(2:end-1) = 2*psdx_acc_z(2:end-1); %dB/Hz
 freq_acc_z = 0:Fs/L:Fs/2;
 % PSD of FFT of Acceleration Plot
 figure (310)
 plot(freq_acc_z,psdx_acc_z,'m')%magenta color
 xlim([0 100])
 grid on
 title_head = sprintf('Periodogram Using FFT of the \n Low-pass Filtered Acceleration- Z Signal');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency ( ((m/s^2)^2)/Hz )','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 fprintf(' \n Figure 310. Periodogram Using FFT of the low-pass filtered Acceleration Signal (((m/s^2)^2)/Hz, Hz) \n ')
 
 
 %the same PSD of Acceleration Plot in 'dB'
 figure (311)
 plot(freq_acc_z,10*log10(psdx_acc_z),'r')%magenta color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 %xlim([0 100])
 grid on
 title_head = sprintf('PSD Using FFT of the \n Low-pass Filtered Acceleration- Z Signal');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 fprintf(' \n Figure 311. PSD Using FFT of the low-pass filtered Acceleration- Z Signal (dB/Hz, Hz) \n ')

 
 % PSD of position signal calculation for position:
 xdft_disp_z = fft(disp_z/1000); %in meters , FFT of the acc is selected for PSD evaluation
 xdft_disp_z = xdft_disp_z(1:(L+1)/2);%where L=length(accfl) not length(acc_signal) not length(accfft)
 psdx_disp_z = (1/(Fs*L)) * abs(xdft_disp_z).^2; %Power is directly proportional with the square of the absolute value of the position amplitude
 psdx_disp_z(2:end-1) = 2*psdx_disp_z(2:end-1); %dB/Hz
 freq_disp_z = 0:Fs/L:Fs/2;
 % PSD of FFT of the position signal Plot
 figure (312)
 plot(freq_disp_z,(psdx_disp_z))%original blue color, dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 xlim([0 100])
 ylim([0 1e-7])
 title('Periodogram Using FFT of the Displacement- Z','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency ((m^2)/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 312. Periodogram Using FFT of Displacement- Z ((m^2)/Hz, Hz) \n ')

 
 %the same PSD of Acceleration Plot in 'dB'
 figure (314)
 plot(freq_disp_z,10*log10(psdx_disp_z),'c')%magenta color
 %xlim([0 100])
 title('Periodogram Using FFT of the Displacement- Z','fontsize',11,'Fontname','Timesnewroman')
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 314. PSD Using FFT of Displacement- Z (dB/Hz, Hz) \n ')
%% CONSLUSIVE WELCH PSD OVER ACCELERATION
%VWelch Power Spectral Estimator for Acceleration and displacement signals:
% Create a Welch spectral estimator for acc and disp.
 h_welch= spectrum.welch;
% Welch PSD of Raw and Filtered Acceleration:
 %One-sided PSD of welch of the raw Acceleration- Z signal:
 Hpsd_raw_acc_z=psd(h_welch,acc_intp_z,'Fs',Fs);
 %One-sided PSD of welch of the low-pass filtered Acceleration- Z signal:
 Hpsd_filtered_acc_z=psd(h_welch,accflfft_z,'Fs',Fs);
 % PSD of Welch spectral of the raw Acceleration- Z Signal Plot:
 figure (315)
 plot(Hpsd_raw_acc_z)
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Raw Acceleration- Z Signal \n Welch Power/Freq, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Raw Metrobus Acceleration- Z Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 315. Welch Power Spectrum of Raw Metrobus Acceleration- Z Signal(dB/Hz, Hz) \n ')

  figure (316)
 plot(Hpsd_filtered_acc_z)
 xlim([0 150])
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Metrobus Low-pass Filtered Acc- Z Signal \n Welch Power/Frequency, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Metrobus Low-pass Filtered Acc- Z Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 316. Welch Power Spectrum of Metrobus Filtered Acceleration- Z Signal (dB/Hz, Hz) \n ')

%% FINAL CONCLUSION OVER WELCH PSD OF DISPLACEMENT

%One-sided PSD of welch of Displacement:
 Hpsd_disp_z=psd(h_welch,disp_z,'Fs',Fs);% Calculate and plot the one-sided PSD.
 % PSD of Welch spectral of Displacement Plot:
 figure (317)
 plot(Hpsd_disp_z)
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Metrobus Displacement- Z Welch Power/Frequency, \n (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('One-sided Welch PSD of Metrobus Displacement- Z','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 317. Welch Power Spectrum of Metrobus Displacement- Z(dB/Hz, Hz) \n ')
 
 %% Plotting Time varying- 3D-Acceleration, 3D-Velocity and 3D-Displacement Signals:
 %Time varying- Acceleration XYZ in m/s^2 
figure (33311)
plot(time_inp, accfl_z,'-b' ,'linewidth',1);
hold on
plot(time_inp, accfl_x,'-r' ,'linewidth',1);
plot(time_inp, accfl_y,'-g' ,'linewidth',0.5);
  title('Filtered Acceleration-X-Y-Z Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Acceleration-X-Y-Z (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 legend('Acceleration-Z (m/s^2)','Acceleration-X (m/s^2)','Acceleration-Y (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
   
 grid on
 fprintf(' \n Figure 333.11. Time-varying Filtered Input Acceleration-X Signal(m/s^2) \n ')
 
 
 %Time varying-INTEGRATED AND FILTERED VELOCITY
figure (33313)
plot(time_inp, velf_z,'-b' ,'linewidth',1);
hold on
plot(time_inp, velf_x,'-r' ,'linewidth',1);
plot(time_inp, velf_y,'-g' ,'linewidth',0.5);
  title('INTEGRATED AND FILTERED VELOCITY-X-Y-Z Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Velocity-X-Y-Z (m/s)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 legend('Velocity-Z (m/s)','Velocity-X (m/s)','Velocity-Y (m/s)','fontsize',11,'Fontname','Timesnewroman');
   
 grid on
 fprintf(' \n Figure 333.13. Time-varying INTEGRATED AND FILTERED VELOCITY-X-Y-Z Signal(m/s) \n ')
 
 %Time varying- ıntegrated and filtered Displacement
 figure (33314)
plot(time_inp, disp_z,'-b' ,'linewidth',1);
hold on
plot(time_inp, disp_x,'-r' ,'linewidth',1);
plot(time_inp, disp_y,'-g' ,'linewidth',0.5);
  title('INTEGRATED AND FILTERED Displacement-X-Y-Z Signal','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Time (Sec)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Displacement-X-Y-Z (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 legend('Displacement-Z (mm)','Displacement-X (mm)','Displacement-Y (mm)','fontsize',11,'Fontname','Timesnewroman');
   
 grid on
 fprintf(' \n Figure 333.14. Time-varying INTEGRATED AND FILTERED Displacement-X-Y-Z Signal(mm) \n ')
  %% FFTs of XYZ acceleration and XYZ displacement signals
 
  %fft of Acceleration plot
 figure (33315)
 plot(f,2*abs(accflfft_z(1:NFFT/2)),'-b' ,'linewidth',1);% For Single-Sided fft: 2*(abs value).
 hold on
 plot(f,2*abs(accflfft_x(1:NFFT/2)),'-r' ,'linewidth',1);
 plot(f,2*abs(accflfft_y(1:NFFT/2)),'-g' ,'linewidth',0.5);
 xlim([0 100])% since low pass filtered around 100 Hz
 
 title_head = sprintf('Single-Sided Low-pass Filtered \n Acceleration-XYZ Amplitude Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|Y(f)|, (m/s^2)','fontsize',11,'Fontname','Timesnewroman');
 %legend('Displacement-Z (mm)','Displacement-X (mm)','Displacement-Y (mm)','fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 333.15. Single-Sided low-pass Filtered Acceleration-XYZ Amplitude Spectrum (m/s^2, Hz)\n ')
 
 %the same fft of acc plot in 'dB'
 figure (33316)

 plot(f,20*log10(2*abs(accflfft_z(1:NFFT/2))),'-b' ,'linewidth',1);% For Single-Sided fft: 2*(abs value).
 hold on
 plot(f,20*log10(2*abs(accflfft_x(1:NFFT/2))),'-r' ,'linewidth',1);
 plot(f,20*log10(2*abs(accflfft_y(1:NFFT/2))),'-g' ,'linewidth',0.5);
 xlim([0 100])% since low pass filtered around 100 Hz
 
 title_head = sprintf('Single-Sided Low-pass Filtered \n Acceleration-XYZ Amplitude Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('20*log10(|Y(f)|), (dB)','fontsize',11,'Fontname','Timesnewroman');% dB=20log10(A1/A2) & dB=10log10(Power1/Power2)
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 %legend('Displacement-Z (mm)','Displacement-X (mm)','Displacement-Y (mm)','fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 333.16. Single-Sided low-pass Filtered Acceleration-XYZ Amplitude Spectrum (dB, Hz) \n ')
 
 
 
 %FFT of Position signal:
 
 figure (33317)
 plot(f,2*abs(dispfft_z(1:NFFT/2)),'-b' ,'linewidth',1);% For Single-Sided fft: 2*(abs value).
 hold on
 plot(f,2*abs(dispfft_x(1:NFFT/2)),'-r' ,'linewidth',1);
 plot(f,2*abs(dispfft_y(1:NFFT/2)),'-g' ,'linewidth',0.5);
 xlim([0 100])% since low pass filtered around 100 Hz
 
 title_head = sprintf('Single-Sided \n Displacement Amplitude-XYZ Spectrum');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('|U(f)|, (mm)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 %legend('Displacement-Z (mm)','Displacement-X (mm)','Displacement-Y (mm)','fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 333.17. Single-Sided Displacement Amplitude-XYZ Spectrum (mm, Hz) \n ')
 
%the same FFT of position plot in 'dB'
 figure (33319)
 
 plot(f,20*log10(2*abs(dispfft_z(1:NFFT/2))),'-b' ,'linewidth',1);% For Single-Sided fft: 2*(abs value).
 hold on
 plot(f,20*log10(2*abs(dispfft_x(1:NFFT/2))),'-r' ,'linewidth',1);
 plot(f,20*log10(2*abs(dispfft_y(1:NFFT/2))),'-g' ,'linewidth',0.5);
 xlim([0 100])% since low pass filtered around 100 Hz

 title('Single-Sided Displacement Amplitude-XYZ Spectrum','fontsize',11,'Fontname','Timesnewroman')
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('20*log10(|U(f)|), (dB)','fontsize',11,'Fontname','Timesnewroman');% dB=20log10(A1/A2) & dB=10log10(Power1/Power2),
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 %legend('Displacement-Z (mm)','Displacement-X (mm)','Displacement-Y (mm)','fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 333.19. Single-Sided Displacement Amplitude-XYZ Spectrum (dB, Hz) \n ')

 %% PSDs of acceleration FFT and displacement FFTs:
 
 % PSD of FFT of Acceleration Plot
 figure (333110)
plot(freq_acc_z,psdx_acc_z,'-b' ,'linewidth',1);
 hold on
 plot(freq_acc_x,psdx_acc_x,'-r' ,'linewidth',1);
 plot(freq_acc_y,psdx_acc_y,'-g' ,'linewidth',0.5);

 xlim([0 100])
 ylim([0 6])
 %legend('Displacement-Z (mm)','Displacement-X (mm)','Displacement-Y (mm)','fontsize',11,'Fontname','Timesnewroman');
 grid on
 title_head = sprintf('Periodogram Using FFT of the \n Low-pass Filtered Acceleration-XYZ Signal');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency ( ((m/s^2)^2)/Hz )','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 fprintf(' \n Figure 333.110. Periodogram Using FFT of the low-pass filtered Acceleration-XYZ Signal (((m/s^2)^2)/Hz, Hz) \n ')
 
 %%
 %the same PSD of Acceleration Plot in 'dB'
 figure (333111)
plot(freq_acc_z,10*log10(psdx_acc_z),'-b' ,'linewidth',1);
 hold on
 plot(freq_acc_x,10*log10(psdx_acc_x),'-r' ,'linewidth',1);
 plot(freq_acc_y,10*log10(psdx_acc_y),'-g' ,'linewidth',0.5);
legend('PSD ACC (dB)-Z','PSD ACC (dB)-X','PSD ACC Z (dB)-Y','fontsize',11,'Fontname','Timesnewroman');
 grid on
 title_head = sprintf('PSD Using FFT of the \n Low-pass Filtered Acceleration-XYZ Signal');
 title(title_head,'fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency (Hz)','fontsize',11,'Fontname','Timesnewroman');
 ylabel('Power/Frequency (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 fprintf(' \n Figure 333.111. PSD Using FFT of the low-pass filtered Acceleration-XYZ Signal (dB/Hz, Hz) \n ')

 %% CONSLUSIVE WELCH PSD OVER ACCELERATION

% Welch PSD of Raw and Filtered Acceleration:
 
 % PSD of Welch spectral of the raw Acceleration-X Signal Plot:
 figure (333115)
 %subplot(1,3,1)
 plot(Hpsd_raw_acc_z);
 hold on
 plot(Hpsd_raw_acc_x);
 plot(Hpsd_raw_acc_y);
legend('\color{blue} WELCH PSD Raw Acceleration-Z (dB/Hz)','\color{red} WELCH PSD Raw Acceleration-X (dB/Hz)','\color{green} WELCH PSD Raw Acceleration-Y (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Raw Acceleration-XYZ Signal \n Welch Power/Freq, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Raw Metrobus Acceleration-XYZ Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 333.115. Welch Power Spectrum of Raw Metrobus Acceleration-XYZ Signal(dB/Hz, Hz) \n ')


%% RMS Values of raw Acceleration signal, filtered Acceleration and Disp:
AccRMSvalue_z=rms(acc_signal_z);
fprintf(' \n RMS value of Metrobus Raw Acceleration-Z (m/s^2)= %.3f \n ',AccRMSvalue_z)
fprintf(' \n RMS value of Metrobus Raw Acceleration-Z (dB)= %.3f \n ',20*log10(AccRMSvalue_z))
AccRMSvalue_x=rms(acc_signal_x);
fprintf(' \n RMS value of Metrobus Raw Acceleration-X (m/s^2)= %.3f \n ',AccRMSvalue_x)
fprintf(' \n RMS value of Metrobus Raw Acceleration-X (dB)= %.3f \n ',20*log10(AccRMSvalue_x))
AccRMSvalue_y=rms(acc_signal_y);
fprintf(' \n RMS value of Metrobus Raw Acceleration-Y (m/s^2)= %.3f \n ',AccRMSvalue_y)
fprintf(' \n RMS value of Metrobus Raw Acceleration-Y (dB)= %.3f \n ',20*log10(AccRMSvalue_y))

FilteredAccRMSvalue_z=rms(accfl_z);
fprintf(' \n RMS value of Metrobus Low Pass Filtered Acceleration Z (m/s^2)= %.3f \n ',FilteredAccRMSvalue_z)
fprintf(' \n RMS value of Metrobus Low Pass Filtered Acceleration Z (dB)= %.3f \n ',20*log10(FilteredAccRMSvalue_z))
FilteredAccRMSvalue_x=rms(accfl_x);
fprintf(' \n RMS value of Metrobus Low Pass Filtered Acceleration X (m/s^2)= %.3f \n ',FilteredAccRMSvalue_x)
fprintf(' \n RMS value of Metrobus Low Pass Filtered Acceleration X (dB)= %.3f \n ',20*log10(FilteredAccRMSvalue_x))
FilteredAccRMSvalue_y=rms(accfl_y);
fprintf(' \n RMS value of Metrobus Low Pass Filtered Acceleration Y (m/s^2)= %.3f \n ',FilteredAccRMSvalue_y)
fprintf(' \n RMS value of Metrobus Low Pass Filtered Acceleration Y (dB)= %.3f \n ',20*log10(FilteredAccRMSvalue_y))



DispRMSvalue_z=rms(disp_z);
fprintf(' \n RMS value of Metrobus Displacement Z(mm)= %.3f \n ',DispRMSvalue_z)
fprintf(' \n RMS value of Metrobus Displacement Z(dB)= %.3f \n ',20*log10(DispRMSvalue_z))
DispRMSvalue_x=rms(disp_x);
fprintf(' \n RMS value of Metrobus Displacement X (mm)= %.3f \n ',DispRMSvalue_x)
fprintf(' \n RMS value of Metrobus Displacement X(dB)= %.3f \n ',20*log10(DispRMSvalue_x))
DispRMSvalue_y=rms(disp_y);
fprintf(' \n RMS value of Metrobus Displacement Y(mm)= %.3f \n ',DispRMSvalue_y)
fprintf(' \n RMS value of Metrobus Displacement Y(dB)= %.3f \n ',20*log10(DispRMSvalue_y))

 
 toc %timing ends