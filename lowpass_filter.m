function [adv_x_filt,adv_y_filt,adv_z_filt,electr_filt,optode_filt,seconds_filt,adv_x_fluc,adv_y_fluc,adv_z_fluc,electr_fluc,optode_fluc] = lowpass_filter(adv_x,adv_y,adv_z,electr,optode,seconds,freqOpt,ai)

% running averaging
   WS2 = (ai/freqOpt);              
   optode(isnan(optode)) = 0;
   Wn = ai/freqOpt;                        
   [b,a] = butter(2,Wn,'low');
   adv_x_filt = filtfilt(b,a,adv_x);
   adv_y_filt = filtfilt(b,a,adv_y);
   adv_z_filt = filtfilt(b,a,adv_z);  
   optode_filt = filtfilt(b,a,optode);  
   %electr_filt = filtfilt(b,a,electr);  
  
   adv_x_fluc = adv_x-adv_x_filt;
   adv_y_fluc = adv_y-adv_y_filt;
   adv_z_fluc = adv_z-adv_z_filt;
   optode_fluc = optode-optode_filt;
   %electr_fluc = electr2-electr_filt;
   
   electr_filt = optode_filt;
   electr_fluc = optode_fluc;
   
   seconds_filt = seconds;
   
   
   %Wn = ai*freqOpt;                        
   %[b,a] = butter(2,Wn,'low');
   %adv_x_filt = filtfilt(b,a,Adv_x1);
   %adv_y_filt = filtfilt(b,a,Adv_x1);
   %adv_y_filt = filtfilt(b,a,Adv_x1);
   %hold on
   %plot(Adv_x1)
   %plot(adv_x_filt) 
      
end

