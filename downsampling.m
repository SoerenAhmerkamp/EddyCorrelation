function [ adv_x,adv_y,adv_z,electr,optode,seconds] = downsampling( Adv_x,Adv_y,Adv_z,freq,freqOpt,Electrode,Optode,Seconds,samples)
%downsampling to optode frequency

WS          =  round(freq/freqOpt/2)*2;                         
adv_x       = (filter(ones(1,WS)/WS,1,Adv_x));
adv_y       = (filter(ones(1,WS)/WS,1,Adv_y));
adv_z       = (filter(ones(1,WS)/WS,1,Adv_z));
Electr      = (filter(ones(1,WS)/WS,1,Electrode));
Optod       = (filter(ones(1,WS)/WS,1,Optode));

% cutting
adv_x(1:1:WS/2)     = [];  
adv_y(1:1:WS/2)     = []; 
adv_z(1:1:WS/2)     = []; 
Optod (1:1:WS/2)    = []; 
Electr(1:1:WS/2)    = []; 
Seconds(1:1:WS/2)   = [];

seconds     = [Seconds(1):1/freqOpt:Seconds(end)]';
ts_adv_x    = timeseries(adv_x,Seconds);
ts_adv_y    = timeseries(adv_y,Seconds);
ts_adv_z    = timeseries(adv_z,Seconds);
ts_Optod    = timeseries(Optod,Seconds);
ts_Electr   = timeseries(Electr,Seconds);

adv_x = resample(ts_adv_x,seconds);
adv_y = resample(ts_adv_y,seconds);
adv_z = resample(ts_adv_z,seconds);
optode = resample(ts_Optod,seconds);
electr = resample(ts_Electr,seconds);

adv_x = adv_x.data;
adv_y = adv_y.data;
adv_z = adv_z.data;
optode = optode.data;
electr = electr.data;
