% This Script requires the signal processing toolbox 
format compact;                     
close all;                           
clear all;

%% Read Vector
load('.\Elba2019_SeaGrass_Eddy2.mat')
plot(VBvel3)
[x,y] = ginput(2);
x(1) = 1;
x(2) = length(VBvel3);

%% Read the full Vector data
freq = 16;                          % frequency ADV
freqOpt = 4;                        % frequency Optode
freq15 = 1/(15*60);                % frequency 15min intervall
first_sample  = round(x(1));              % first samples to be analyzed
last_sample = round(x(2));             % choose manually, last samples to be analyzed (22.5 hours later)
samples = last_sample-first_sample+1;
ai = 120;                           % averaging interval in seconds

Burst       = VTburst(1,first_sample:last_sample);
Adv_x2      = VBvel1(1,first_sample:last_sample);
Adv_y2      = VBvel2(1,first_sample:last_sample);
Adv_z2     = -VBvel3(1,first_sample:last_sample);
Amp1        = VBamp1(1,first_sample:last_sample);
Amp2        = VBamp2(1,first_sample:last_sample);
Amp3        = VBamp3(1,first_sample:last_sample);
Corr1       = VBcorr1(1,first_sample:last_sample);
Corr2       = VBcorr2(1,first_sample:last_sample);
Corr3       = VBcorr3(1,first_sample:last_sample);
Press       = VTpressure(1,first_sample:last_sample);
Analog1     = VTanain(1,first_sample:last_sample);
Analog2     = VTanain(2,first_sample:last_sample);
Seconds     = (first_sample:last_sample)/freq; % creates an index used as time vector
hours       = Seconds/3600;
plot(Analog1)

%% Read Aanderaa Oxygen and Temperature

filename = '.\aanderaSensor260619190622.log';
[timestampAandT,aandOxy,aandTemp] = lance_aand(filename);
timestampAand = zeros(1,length(timestampAandT));
for n = 1:length(timestampAand)-1
    temp_timestamp = cell2mat(timestampAandT(n));
    tempd = strsplit(temp_timestamp(1:8),'/');
    timestampAand(n) = datenum(strcat('20',tempd(3),'/',tempd(2),'/',tempd(1),temp_timestamp(9:end)));
end
aandTemp(1:10) = NaN;
aandOxy(1:10) = NaN;

plot(aandTemp)
%% Read Par Sensor
filename = 'D:\Seafile\Projects\Elba\Elba2019\Seegrass_Elba2019-2\Logs\parSensor260619190619.log ';
[timestampParT,par] = lance_par(filename);
timestampPar = zeros(1,length(timestampParT));
for n = 1:length(timestampPar)-1
    temp_timestampPar = cell2mat(timestampParT(n));
    tempd = strsplit(temp_timestampPar(1:8),'/');
    timestampPar(n) = datenum(strcat('20',tempd(3),'/',tempd(2),'/',tempd(1),temp_timestampPar(9:end)));
end

plot(par);

%% Calibrate O2
% Pyroscience matlab toolbox required

oxyEddy = zeros(1,length(Analog1));
zero_dphi = 15;     % at 0 mV
max_dphi  = 30;     % at 2500 mV (max_mV)
zero_count = 32323; % count of Vector at zero input
dphi = (Analog1-zero_count)*5000/65536*2/2500*(max_dphi-zero_dphi)+zero_dphi;

timestampEddy = datenum(clockdeploy)+(1:length(dphi))/16/3600/24;
timestampEddyT = datenum(clockdeploy)+(1:length(dphi))/3600/24;
ind = diff(timestampAand)==0;
temperatureEddy = interp1(timestampAand(~ind),aandTemp(~ind),timestampEddy)./1000;
setCalibrationData(27.1,27.1,1027,100,40,20.6)
for n = 1:length(oxyEddy)
    oxTemp = calculateOxygenZYW(dphi(n),temperatureEddy(n),1027,38);
    oxyEddy(n) = oxTemp(1);
end

plot(oxyEddy)

%% function despiking ADV Data
Adv_x1 = despikeADV(Adv_x2');
Adv_y1 = despikeADV(Adv_y2');
Adv_z1 = despikeADV(Adv_z2');

%% planar fit

% function Eddy_Planar_Fit (rotate coordinate system)
[Adv_x,Adv_y,Adv_z,Tilt_x_z_plane,Tilt_y_z_plane,Z_shift] = planar_fit(Adv_x1,Adv_y1,Adv_z1); 

%% downsampling
filterFreq = 0.01;
Optode = oxyEddy;
Electrode = Optode;

% function downsampling to Optode frequency
[adv_x,adv_y,adv_z,electr,optode,seconds] = downsampling( Adv_x,Adv_y,Adv_z,freq,freqOpt,Electrode,Optode',Seconds',samples);

% function running average (according to averaging interval 'ai')
[adv_x_filt,adv_y_filt,adv_z_filt,electr_filt,optode_filt,seconds_filt,adv_x_fluc,adv_y_fluc,adv_z_fluc,electr_fluc,optode_fluc] = running_average(adv_x,adv_y,adv_z,electr,optode,seconds,freqOpt,ai);
%[adv_x_filt,adv_y_filt,adv_z_filt,electr_filt,optode_filt,seconds_filt,adv_x_fluc,adv_y_fluc,adv_z_fluc,electr_fluc,optode_fluc] = lowpass_filter(adv_x,adv_y,adv_z,electr,optode,seconds,freqOpt,filterFreq);

hours_filt   = seconds_filt/3600;
[Dir,Mag]    = cart2pol(adv_x_filt,adv_y_filt); %calculate direction and magnitude
Dir          = Dir/pi/2*360;                    %direction to degree
   
%% BURST ANALYSIS
burst_min       = 30;                                       % define length of individual burst in minutes
burst_length    = freqOpt*60*burst_min;                     % no of samples per burst                 
burst_no        = floor(length(adv_z_fluc)/burst_length);   % no of bursts 
shift           = 47;                                       % uneven number of data points for time shift analysis!! 
window          = 2*floor(burst_length/16);                 % defining the window for the spectral analysis
noverlap        = window/2;
nfft            = window;
fs              = freqOpt;

%Define Variables 
u           = zeros(burst_length,burst_no);
v           = zeros(burst_length,burst_no);
w           = zeros(burst_length,burst_no);
o2opt       = zeros(burst_length,burst_no);
o2opt_mean       = zeros(burst_length,burst_no);
o2elec      = zeros(burst_length,burst_no);
sec         = zeros(burst_length,burst_no);
hour        = zeros(burst_length,burst_no);
u_mean      = zeros(burst_length,burst_no);
v_mean      = zeros(burst_length,burst_no);
w_mean      = zeros(burst_length,burst_no);
u_raw      = zeros(burst_length,burst_no);
v_raw      = zeros(burst_length,burst_no);
w_raw      = zeros(burst_length,burst_no);
o2_raw      = zeros(burst_length,burst_no);
orbit       = zeros(1,burst_no);
dir         = zeros(burst_length,burst_no);
mag         = zeros(burst_length,burst_no);
P_Flux_opt  = zeros(shift,burst_no);  
R_Flux_opt  = zeros(shift,burst_no); 
Mean_flux_max_opt   = zeros(burst_no,11);
flux_max_opt        = zeros(burst_length-shift,burst_no);  
cumflux_max_opt     = zeros(burst_length-shift,burst_no);  
flux_vw             = zeros(burst_length,burst_no);
cumflux_vw          = zeros(burst_length,burst_no);
ampAvg         = zeros(burst_length,burst_no);

% cut out burst time series
for n=1:1:burst_no-2;                       
u(:,n)       = adv_x_fluc((n-1)*burst_length+1:n*burst_length);
v(:,n)       = adv_y_fluc((n-1)*burst_length+1:n*burst_length);
w(:,n)       = adv_z_fluc((n-1)*burst_length+1:n*burst_length);
o2opt(:,n)   = optode_fluc((n-1)*burst_length+1:n*burst_length);
o2opt_mean(:,n)   = optode_filt((n-1)*burst_length+1:n*burst_length);
o2elec(:,n)  = electr_fluc((n-1)*burst_length+1:n*burst_length);
sec(:,n)     = seconds_filt((n-1)*burst_length+1:n*burst_length);
hour(:,n)    = hours_filt((n-1)*burst_length+1:n*burst_length);
dir(:,n)     = Dir((n-1)*burst_length+1:n*burst_length);
mag(:,n)     = Mag((n-1)*burst_length+1:n*burst_length);
u_mean(:,n)  = adv_x_filt((n-1)*burst_length+1:n*burst_length);
v_mean(:,n)  = adv_y_filt((n-1)*burst_length+1:n*burst_length);
w_mean(:,n)  = adv_z_filt((n-1)*burst_length+1:n*burst_length);
u_raw(:,n)   = adv_x((n-1)*burst_length+1+ai*freqOpt/2:n*burst_length+ai*freqOpt/2);
v_raw(:,n)   = adv_y((n-1)*burst_length+1+ai*freqOpt/2:n*burst_length+ai*freqOpt/2);
w_raw(:,n)   = adv_z((n-1)*burst_length+1+ai*freqOpt/2:n*burst_length+ai*freqOpt/2);
ampAvg(:,n)   = Amp1((n-1)*burst_length+1+ai*freqOpt/2:n*burst_length+ai*freqOpt/2);
orbit(n)   = sqrt(2*(var(adv_x_fluc((n-1)*burst_length+1:n*burst_length))+var(adv_z_fluc((n-1)*burst_length+1:n*burst_length))));
o2_raw(:,n)  = optode((n-1)*burst_length+1+ai*freqOpt/2:n*burst_length+ai*freqOpt/2);

for t=1:1:shift
    [R,PP]   = corrcoef(w(((shift+1)/2):1:end-((shift+1)/2),n),o2opt(t:1:end-(shift+1)+t,n));
    R_Flux_opt(t,n)           = R(1,2);     % correlation
    P_Flux_opt(t,n)           = PP(1,2);    % p-value 
end

% finding the time shift with the highest correlation, then calculating flux and cumulative flux
[P_min_opt,I_p]         = min(P_Flux_opt(:,n));
flux_max_opt(:,n)       = w(((shift+1)/2):1:end-((shift+1)/2),n).*o2opt(I_p:1:end-(shift+1)+I_p,n)*86400;
cumflux_max_opt(:,n)    = cumsum(flux_max_opt(:,n));

% calculating the reynolds stress
flux_vw(:,n)            = w(:,n).*v(:,n);
cumflux_vw(:,n)         = cumsum(flux_vw(:,n));

% calculating the cumulative co-spectra
[P,FF]                  = cpsd(o2opt(I_p:1:end-(shift+1)+I_p,n),w(((shift+1)/2):1:end-((shift+1)/2),n),window,noverlap,nfft,freqOpt);
P                       = flipud(real(P));
P                       = cumtrapz(FF,real(P));

% calculating the power spectra for velocities and sensors
[PPu,ff]                = pwelch(u(:,n),window,noverlap,nfft,freqOpt);
Pu(:,n)                 = real(PPu);
[PPv,ff]                = pwelch(v(:,n),window,noverlap,nfft,freqOpt);
Pv(:,n)                 = real(PPv);
[PPw,ff]                = pwelch(w(:,n),window,noverlap,nfft,freqOpt);
Pw(:,n)                 = real(PPw);
[PPo2,ff]               = pwelch(o2opt(:,n),window,noverlap,nfft,freqOpt);
Po2opt(:,n)             = real(PPo2);

corr   = R_Flux_opt(I_p,n);               
mean_flux_shift   = mean(w(((shift+1)/2):1:end-((shift+1)/2),n).*o2opt(I_p:1:end-(shift+1)+I_p,n))*86400; 
mean_flux   = mean(w(:,n).*o2opt(:,n))*86400; 
end
K = ff.^(-5/3).*0.0001;     


%% analysis of entire time series

shift_final = 0; % put number from above
size_fluc = length(adv_z_fluc);
optode_fluc(1:1:shift_final,:)                      =  []; 
adv_z_fluc(size_fluc-shift_final+1:1:size_fluc,:)   =  []; 

[R,PP]   = corrcoef(optode_fluc,adv_z_fluc);                %correlation and probability
   
flux_opt_x = optode_fluc.*adv_x_fluc;
flux_opt_y = optode_fluc.*adv_y_fluc;
flux_opt_z = optode_fluc.*adv_z_fluc;                      %calculate instantaneous flux
mean_flux_opt_x = nanmean(flux_opt_x)*86400;                    %calculate mean flux
mean_flux_opt_y = nanmean(flux_opt_y)*86400;                    %calculate mean flux
mean_flux_opt_z = nanmean(flux_opt_z)*86400;                    %calculate mean flux
cumflux_opt_x = nancumsum(flux_opt_x)/freqOpt;                  %calculate cumulative flux
cumflux_opt_y = nancumsum(flux_opt_y)/freqOpt;                  %calculate cumulative flux
cumflux_opt_z = nancumsum(flux_opt_z)/freqOpt;                  %calculate cumulative flux

%% Figures here shifted and non-shifted time series
figure; bar(mean_flux);
figure; bar(mean_flux_shift);

