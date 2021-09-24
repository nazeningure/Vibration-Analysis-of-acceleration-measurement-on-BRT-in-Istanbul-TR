
 %% NEW METHOD ON WELCH PSD ESTIMATION ACC XYZ
 %pwelch command gives (w/Fs) vs power/(w). Thus, "[2*pi*f/Fs]*(Fs/(2pi)" vs "[power/2*pi*f]*2pi" gives previously estimated standard PSD units. Therefore results can be compared.
 figure (3331)

[PX,FX]=pwelch(acc_intp_x,Fs,'power');
[PY,FY]=pwelch(acc_intp_y,Fs,'power');
[PZ,FZ]=pwelch(acc_intp_z,Fs,'power');

 plot(FZ/(2*pi)*Fs,PZ*2*pi,'-b' ,'linewidth',1)%pwelch command gives (w/Fs) vs power/(w). Thus, "[2*pi*f/Fs]*(Fs/(2pi)" vs "[power/2*pi*f]*2pi" gives previously estimated standard PSD units. Therefore results can be compared.
 hold on
 plot(FX/(2*pi)*Fs,PX*2*pi,'-r' ,'linewidth',1)
 plot(FY/(2*pi)*Fs,PY*2*pi,'-g' ,'linewidth',1.5)
legend('\color{blue} WELCH PSD Raw Acceleration-Z ((m/s^2)^2/Hz)','\color{red} WELCH PSD Raw Acceleration-X ((m/s^2)^2/Hz)','\color{green} WELCH PSD Raw Acceleration-Y ((m/s^2)^2/Hz)','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Raw Acceleration-XYZ Signal \n Welch Power/Freq, ((m/s^2)^2/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Raw Metrobus Acceleration-XYZ Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 333.1. Welch Power Spectrum of Raw Metrobus Acceleration-XYZ Signal((m/s^2)^2/Hz, Hz) \n ')
 
 
  %NEW METHOD ON WELCH PSD ESTIMATION ACC XYZ
 %pwelch command gives (w/Fs) vs power/(w). Thus, "[2*pi*f/Fs]*(Fs/(2pi)" vs "[power/2*pi*f]*2pi" gives previously estimated standard PSD units. Therefore results can be compared.
 [PX,FX]=pwelch(acc_intp_x,Fs,'power');
[PY,FY]=pwelch(acc_intp_y,Fs,'power');
[PZ,FZ]=pwelch(acc_intp_z,Fs,'power');
 
 figure (33311)

 plot(FZ/(2*pi)*Fs,10*log10(PZ*2*pi),'-b' ,'linewidth',1)%pwelch command gives (w/Fs) vs power/(w). Thus, "[2*pi*f/Fs]*(Fs/(2pi)" vs "[power/2*pi*f]*2pi" gives previously estimated standard PSD units. Therefore results can be compared.
 hold on
 plot(FX/(2*pi)*Fs,10*log10(PX*2*pi),'-r' ,'linewidth',1)
 plot(FY/(2*pi)*Fs,10*log10(PY*2*pi),'-g' ,'linewidth',1.5)
legend('\color{blue} WELCH PSD Raw Acceleration-Z (dB/Hz)','\color{red} WELCH PSD Raw Acceleration-X (dB/Hz)','\color{green} WELCH PSD Raw Acceleration-Y (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Raw Acceleration-XYZ Signal \n Welch Power/Freq, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Raw Metrobus Acceleration-XYZ Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 333.11. Welch Power Spectrum of Raw Metrobus Acceleration-XYZ Signal(dB/Hz, Hz) \n ')
 

 
 %% %NEW METHOD ON WELCH PSD ESTIMATION DISP XYZ
 %pwelch command gives (w/Fs) vs power/(w). Thus, "[2*pi*f/Fs]*(Fs/(2pi)" vs "[power/2*pi*f]*2pi" gives previously estimated standard PSD units. Therefore results can be compared.
 figure (3332)

[PX,FX]=pwelch(disp_x,Fs,'power');
[PY,FY]=pwelch(disp_y,Fs,'power');
[PZ,FZ]=pwelch(disp_z,Fs,'power');

 plot(FZ/(2*pi)*Fs,PZ*2*pi,'-b' ,'linewidth',1)%pwelch command gives (w/Fs) vs power/(w). Thus, "[2*pi*f/Fs]*(Fs/(2pi)" vs "[power/2*pi*f]*2pi" gives previously estimated standard PSD units. Therefore results can be compared.
 hold on
 plot(FX/(2*pi)*Fs,PX*2*pi,'-r' ,'linewidth',1)
 plot(FY/(2*pi)*Fs,PY*2*pi,'-g' ,'linewidth',1.5)
legend('\color{blue} WELCH PSD Raw Displacement-Z ((mm^2)/Hz)','\color{red} WELCH PSD Raw Displacement-X ((mm^2)/Hz)','\color{green} WELCH PSD Raw Displacement-Y ((mm^2)/Hz)','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Raw Displacement-XYZ Signal \n Welch Power/Freq, ((mm^2)^2/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Raw Metrobus Displacement-XYZ Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 333.2. Welch Power Spectrum of Raw Metrobus Displacement-XYZ Signal((mm^2)/Hz, Hz) \n ')

  %NEW METHOD ON WELCH PSD ESTIMATION DISP XYZ
 %pwelch command gives (w/Fs) vs power/(w). Thus, "[2*pi*f/Fs]*(Fs/(2pi)" vs "[power/2*pi*f]*2pi" gives previously estimated standard PSD units. Therefore results can be compared.
 figure (33322)
 plot(FZ/(2*pi)*Fs,10*log10(PZ*2*pi),'-b' ,'linewidth',1)%pwelch command gives (w/Fs) vs power/(w). Thus, "[2*pi*f/Fs]*(Fs/(2pi)" vs "[power/2*pi*f]*2pi" gives previously estimated standard PSD units. Therefore results can be compared.
 hold on
 plot(FX/(2*pi)*Fs,10*log10(PX*2*pi),'-r' ,'linewidth',1)
 plot(FY/(2*pi)*Fs,10*log10(PY*2*pi),'-g' ,'linewidth',1.5)
legend('\color{blue} WELCH PSD Raw Displacement-Z (dB/Hz)','\color{red} WELCH PSD Raw Displacement-X (dB/Hz)','\color{green} WELCH PSD Raw Displacement-Y (dB/Hz)','fontsize',11,'Fontname','Timesnewroman');
 xlabel('Frequency, (Hz)','fontsize',11,'Fontname','Timesnewroman');
 yaxis_head = sprintf('Raw Displacement-XYZ Signal \n Welch Power/Freq, (dB/Hz)');
 ylabel(yaxis_head,'fontsize',11,'Fontname','Timesnewroman');
 title('Welch PSD of Raw Metrobus Displacement-XYZ Signal','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on;
 fprintf(' \n Figure 333.22. Welch Power Spectrum of Raw Metrobus Displacement-XYZ Signal(dB/Hz, Hz) \n ')