function [par2]=getTmin2(par)
% Function to create histogram of species cold tolerances (Tmin)
load('path_to_obis_Tmin/OBIS_Thab.mat','obis');
Thab10 = obis.Thab10(obis.nobs>par.nomin);

% PDF of minimum temperature from OBIS
par2.Tmin2 = [-4.5 -1.5 1.5 4.5 7.5 10.5 13.5 16.5 19.5 22.5 25.5 28.5 31.5];% bin edges 
bc = histc(Thab10(:)',par2.Tmin2);% bincounts per bin
par2.Tmw = bc./nanmax(bc(:)); % normalized

% bin mid points: minimum temperature range (-3oC to 30oC)
par2.Tmin =-3:3:30; 

% extract bincounts
par2.Tmw = par2.Tmw(1:end-1); 

clear obis

