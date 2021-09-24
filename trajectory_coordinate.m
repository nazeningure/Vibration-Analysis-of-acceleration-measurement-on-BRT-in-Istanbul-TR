clear all; close all; clc;

filename = 'C:\Users\Aleyna\Documents\MATLAB\coordinate_sb.xlsx';
[filepath,name,ext] = fileparts(filename)

%Excel'den datayı varible olarak matlab'e aktarma
%filename = 'BS_YRS1.xlsx'; %path görebilmek için filename üzerinden aynı dosyayı uzantısı ile çektik
data =xlsread(filename);
time =xlsread(filename,'A:A');
alt =xlsread(filename,'B:B');
lat =xlsread(filename,'C:C');
long =xlsread(filename,'D:D');


figure (1)
plot3(lat,long,alt,'b' ,'linewidth',2)
title('YELLOW CONNECTO MERCEDES METROBUS: SOGUTLUCESME-BEYLIKDUZU TRAJECTORY','fontsize',11,'Fontname','Timesnewroman');
xlabel('Lattitude (degree)','fontsize',11,'Fontname','Timesnewroman');
ylabel('Longitude (degree)','fontsize',11,'Fontname','Timesnewroman');
zlabel('Normalized Altitude at Sea Level(m)','fontsize',11,'Fontname','Timesnewroman');
 set(gca,'fontsize',11,'Fontname','Timesnewroman');
 grid on
 fprintf(' \n Figure 1. Metrobus route from Sogutlucesme to Beylikduzu \n ')
 
