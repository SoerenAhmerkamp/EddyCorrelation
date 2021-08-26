function [adv_x_filt,adv_y_filt,adv_z_filt,electr_filt,optode_filt,seconds_filt,adv_x_fluc,adv_y_fluc,adv_z_fluc,electr_fluc,optode_fluc] = running_average(adv_x,adv_y,adv_z,electr,optode,seconds,freqOpt,ai)

% running averaging
   WS2      =    (freqOpt*ai);                   
   adv_x_filt    =    (filter(ones(1,WS2)/WS2,1,adv_x));
   adv_y_filt    =    (filter(ones(1,WS2)/WS2,1,adv_y));
   adv_z_filt    =    (filter(ones(1,WS2)/WS2,1,adv_z));
   optode_filt   =    (filter(ones(1,WS2)/WS2,1,optode));
   electr_filt   =    (filter(ones(1,WS2)/WS2,1,electr));
   seconds_filt   =    (filter(ones(1,WS2)/WS2,1,seconds));
   adv_x_filt(1:1:WS2)    =  [];
   adv_y_filt(1:1:WS2)    =  [];
   adv_z_filt(1:1:WS2)    =  [];
   optode_filt(1:1:WS2)   =  [];
   electr_filt(1:1:WS2)   =  [];
   seconds_filt(1:1:WS2)   =  [];
   
   adv_x2       = adv_x;
   adv_y2       = adv_y;
   adv_z2       = adv_z;
   optode2      = optode;
   electr2      = electr;
   
   adv_x2(1:1:(WS2/2))     = []; 
   adv_y2(1:1:(WS2/2))     = []; 
   adv_z2(1:1:(WS2/2))     = []; 
   optode2(1:1:(WS2/2))    = [];
   electr2(1:1:(WS2/2))    = [];
      
   size=length(adv_z);
   adv_x2((size-WS2+1):1:size-WS2/2)           = []; 
   adv_y2((size-WS2+1):1:size-WS2/2)           = []; 
   adv_z2((size-WS2+1):1:size-WS2/2)           = []; 
   optode2((size-WS2+1):1:size-WS2/2)          = []; 
   electr2((size-WS2+1):1:size-WS2/2)          = []; 
      
   adv_x_fluc = adv_x2-adv_x_filt;
   adv_y_fluc = adv_y2-adv_y_filt;
   adv_z_fluc = adv_z2-adv_z_filt;
   optode_fluc = optode2-optode_filt;
   electr_fluc = electr2-electr_filt;
      
end

