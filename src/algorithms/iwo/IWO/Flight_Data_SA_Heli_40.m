function [in_H, out_H, data, time, time_RC] = Flight_Data_SA_Heli_40(i)

%-----------------------Loading Flight Data--------------------------------

if i == 1
    load('C:\Users\Pranjal Biswas\Documents\MATLAB\HELI\data\09_01_17.bin-342424.mat') %Test Data
elseif i == 2
    load('C:\Users\Pranjal Biswas\Documents\MATLAB\HELI\data\12_10_16.bin-497185.mat') %Experimental Data
else
    load('C:\Users\Pranjal Biswas\Documents\MATLAB\HELI\data\1.mat') %Test Data
end

%-----------------------Creating Time Space--------------------------------

sz_TIME=size(TIME);
t_start=round(TIME(1,2)/10^6);
t_end=round(TIME(sz_TIME(1,1),2)/10^6);
t_flight=t_end-t_start;

%----------------Calculating Frquency of Required Data---------------------

freq_TIME = sz_TIME(1,1)/t_flight;
sampling_time = 1/freq_TIME;
sz_ATT = size(ATT);
freq_ATT = sz_ATT(1,1)/t_flight;
sampling_time_ATT = 1/freq_ATT;
sz_LPOS = size(LPOS);
freq_LPOS = sz_LPOS(1,1)/t_flight;
sampling_time_LPOS = 1/freq_LPOS;
sz_RC = size(RC);
freq_RC = sz_RC(1,1)/t_flight;
sampling_time_RC = 1/freq_RC;
t=0:sampling_time:t_flight-sampling_time;
t_RC=0:sampling_time_RC:t_flight-sampling_time_RC;

%--------------------Calculating Sampling Factor---------------------------

sampling_factor_ATT = sz_TIME(1,1)/sz_ATT(1,1);
sampling_factor_LPOS = sz_TIME(1,1)/sz_LPOS(1,1);
sampling_factor_RC = sz_TIME(1,1)/sz_RC(1,1);

%----------------------Resampling Flight Data------------------------------

[num,den] = numden(sym(sampling_factor_ATT));
num_ATT = double(num);
den_ATT = double(den);
[num,den] = numden(sym(sampling_factor_LPOS));
num_LPOS = double(num);
den_LPOS = double(den);
[num,den] = numden(sym(sampling_factor_RC));
num_RC = double(num);
den_RC = double(den);
ATT_rs = resample(ATT,num_ATT,den_ATT);
LPOS_rs = resample(LPOS,num_LPOS,den_LPOS);
RC_rs = resample(RC,num_RC,den_RC);

%------------Calculating Frquency of Resampled Data Sets-------------------

sz_ATT_rs=size(ATT_rs);
sz_LPOS_rs=size(LPOS_rs);
sz_RC_rs=size(RC_rs);
freq_ATT_rs=sz_ATT_rs(1,1)/t_flight;
freq_LPOS_rs=sz_LPOS_rs(1,1)/t_flight;
freq_RC_rs=sz_RC_rs(1,1)/t_flight;

%---------------------Signal Processing(Filtering Noise)-------------------

Fs = freq_TIME;                    % Sampling frequency

Fc = 1;                            % Cut Off Frequency
[B,A] = butter(3,Fc/(Fs/2));

f_d_lat = filter(B,A,RC_rs(:,2));
f_d_lon = filter(B,A,RC_rs(:,3));
f_d_ped = filter(B,A,RC_rs(:,5));
f_d_col = filter(B,A,RC_rs(:,7));
input=[f_d_lat f_d_lon f_d_ped f_d_col];

f_u = filter(B,A,LPOS_rs(:,8));
f_v = filter(B,A,LPOS_rs(:,7));
f_p = filter(B,A,ATT_rs(:,9));
f_q = filter(B,A,ATT_rs(:,10));
f_phi = filter(B,A,ATT_rs(:,6));
f_theta = filter(B,A,ATT_rs(:,7));
f_w = filter(B,A,LPOS_rs(:,9));
f_r = filter(B,A,ATT_rs(:,11));
output = [f_u f_v f_p f_q f_phi f_theta f_w f_r];

% % f_d_lat = RC_rs(:,2);
% % f_d_lon = RC_rs(:,3);
% % f_d_ped = RC_rs(:,5);
% % f_d_col = RC_rs(:,7);
% % input=[f_d_lat f_d_lon f_d_ped f_d_col];
% % 
% % f_u = LPOS_rs(:,8);
% % f_v = LPOS_rs(:,7);
% % half_power_bw_p = powerbw(ATT_rs(:,9),Fs);
% % d1 = designfilt('lowpassiir','FilterOrder',12,'HalfPowerFrequency',half_power_bw_p,'DesignMethod','butter');
% % f_p = filtfilt(d1,ATT_rs(:,9));
% % half_power_bw_q = powerbw(ATT_rs(:,10),Fs);
% % d2 = designfilt('lowpassiir','FilterOrder',12,'HalfPowerFrequency',half_power_bw_q,'DesignMethod','butter');
% % f_q = filtfilt(d2,ATT_rs(:,10));
% % f_phi = ATT_rs(:,6);
% % f_theta = ATT_rs(:,7);
% % f_w = LPOS_rs(:,9);
% % f_r = ATT_rs(:,11);
% % output = [f_u f_v f_p f_q f_phi f_theta f_w f_r];

%-------------------------------Output-------------------------------------

in_H = input;
out_H = output;
time = t';
time_RC = t_RC;

out(:,1) = f_u;
out(:,2) = f_v;
out(:,3) = f_p;
out(:,4) = f_q;
out(:,5) = f_phi;
out(:,6) = f_theta;
out(:,7) = 0;
out(:,8) = 0;
out(:,9) = f_w;
out(:,10) = f_r;
out(:,11) = 0;
out(:,12) = 0;
out(:,13) = 0;

data = iddata(out,in_H,0.0081);


end