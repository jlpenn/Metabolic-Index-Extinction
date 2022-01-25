% Code to unpack model output
clear all; close all; 
 
% Load model output
load('6me_output.mat')
 
%% Time  (year CE)
% Models run to 2100
time = output.time; 
 
% Models run ti 2300
time_extended = output.time_extended;
 
%% Latitude grid 
% 1 degree resolution
lat1deg = output.lat1deg;
 
% 5 degree resolution (for Richness)
lat5deg = output.lat5deg;
 
%% High emissions scenario
 
% CMIP5/6 models
high.models = output.high.models;
 
% Extinction (%) [row = model; column = time] 
high.ext_loss = output.high.ext_loss; % Habitat loss only
high.ext_lossgain = output.high.ext_lossgain; % Habitat loss+gain
 
% Extirpation (%) [row = model; column = time]
high.extr = output.high.extr;
 
% Global air temperature change (deg. C)
high.dT = output.high.dT; 
 
% Map of extirpation (%) in 2100
high.extr_2100_map = output.high.extr_2100_map;
 
% Extirpation (%) versus latitude
high.extr_lat_2100 = output.high.extr_lat_2100; % mean in 2100
high.extr_lat_2100_sd = output.high.extr_lat_2100_sd; % S.D. in 2100
high.extr_lat_2300 = output.high.extr_lat_2300; % mean in 2300
high.extr_lat_2300_sd = output.high.extr_lat_2300_sd; % S.D. in 2300

% Extinction (%) versus latitude 
high.ext_lat_2300 = output.high.ext_lat_2300; % mean in 2300
high.ext_max_lat_2300 = output.high.ext_max_lat_2300; % maximum in 2300
high.ext_min_lat_2300 = output.high.ext_min_lat_2300; % minimum in 2300

% Future richness changes versus time
high.richness = output.high.rich_ens_time; %mean 
high.richness_high = output.high.rich_high_time;% mean+S.D.
high.richness_low = output.high.rich_low_time; % mean - S.D.


 %% Low emissions scenario
 
 % CMIP5/6 models
low.models = output.low.models;
 
% Extinction (%) [row = model; column = time] 
low.ext_loss = output.low.ext_loss; % Habitat loss only
low.ext_lossgain = output.low.ext_lossgain; % Habitat loss+gain


% Extirpation (%) [row = model; column = time]
low.extr = output.low.extr;
 
% Global air temperature change (deg. C)
low.dT = output.low.dT; 
 

% Extirpation (%) versus latitude
low.extr_lat_2100 = output.low.extr_lat; % mean in 2100
low.extr_lat_2100_sd = output.low.extr_lat_sd; % S.D.  in 2100

% Future richness changes versus time
low.richness = output.low.rich_ens_time; % mean extinction 
low.richness_high = output.low.rich_high_time; % mean extinction + S.D.
low.richness_low = output.low.rich_low_time;% mean extinction - S.D.

%% Historical emissions
% Extirpation (%) versus latitude
hist.extr_lat = output.hist.extr_lat; % mean in 2020
hist.extr_lat_sd = output.hist.extr_lat_sd; % S.D.  in 2020


% Marine biological richness vs. latitude (model) 
hist.rich_lat_0m = output.hist.rich_lat_0m; % surface ocean only
hist.rich_lat_500m = output.hist.rich_lat_500m; % 0-500 m
hist.rich_lat_5000m = output.hist.rich_lat_5000m; % 0-5000 m

%% Observations (fossil record, richness, fisheries data)
% Sepkoski fossil record data from Rohde and Muller, Cycles in Fossil Diversity,Nature (2005)
% PBDB fossil record data from Alroy et al., Phanerozoic trends in the global diversity of marine invertebrates, Science (2008)
% Observed richness from Chadhary, Saeedi, Costello, Marine Species Richness is bimodal with latitude: A reply to Fendandez and Marques,Trends in Ecology and Evolution (2017)
% SeaAroundUs.org locations of productive fisheries (1950-2014)
obs = output.obs;
