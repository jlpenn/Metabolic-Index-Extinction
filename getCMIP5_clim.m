%% Make mulit-decadal climatologies of CMIP5 temperature and oxygen anomalies
clear all;close all

tic

% RCP
do.rcp = 2; % rcp 8.5 (8) or rcp 2.6 (2)

% CMIP models
models = {'CESM';'GFDL';'HADL';'IPSL';'MPI'};

for mm = 1:length(models)
    
    % model
    model_name = models{mm}
    
    % sampling years
    if do.rcp == 8
        if strcmp(model_name,'CESM') || strcmp(model_name,'GFDL')
            year0 = 1900; % start year
            years = [2005:2035; 2035:2065; 2070:2100] - year0;
        else
            year0 = 1900; % start year
            years = [2005:2035; 2035:2065; 2085:2115; 2135:2165; 2185:2215; 2235:2265; 2270:2300] - year0;
        end
    else
        if strcmp(model_name,'MPI') || strcmp(model_name,'IPSL')
            year0 = 1911; % start year
            years = [2005:2035; 2035:2065; 2070:2100] - year0;
        else
            year0 = 1900; % start year
            years = [2005:2035; 2035:2065; 2070:2100] - year0;
        end
    end
    %% Set up environmental variables + measurements
    
    % CMIP5 files
    if do.rcp == 8
        
        % load annual climate anomalies from CMIP (anomaly relative to time-mean 1850-1900; dimensions: time,depth,latitude,longitude)
        o2_path = sprintf('path_to_cmip5_rcp85/cmip5_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip5_rcp85/cmip5_dT_%s.nc',model_name);
        scenario = sprintf('Business as usual: RCP 8.5')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        DO2(DO2<-1e3)=nan;
        DT(DT<-1e3)=nan;
        
        % make climatology
        for tt = 1:(size(years,1)+1)
            if tt == (size(years,1)+1)
                tclim = size(DT,1)-10:size(DT,1);
                dO2(tt,:,:,:) = nanmean(DO2(tclim,:,:,:));
                dT(tt,:,:,:) = nanmean(DT(tclim,:,:,:));
                par.time(tt) = nanmean(tclim);
            else
                tclim = years(tt,:);
                dO2(tt,:,:,:) = nanmean(DO2(tclim,:,:,:));
                dT(tt,:,:,:) = nanmean(DT(tclim,:,:,:));
                par.time(tt) = nanmean(tclim);
            end
        end
        
        clear DT DO2
        
    elseif do.rcp == 2
        % load annual climate anomalies from CMIP (anomaly relative to time-mean 1850-1900; dimensions: time,depth,latitude,longitude)
        o2_path = sprintf('path_to_cmip5_rcp26/cmip5_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip5_rcp26/cmip5_dT_%s.nc',model_name);
        scenario = sprintf('Best case: RCP 2.6')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        DO2(DO2<-1e3)=nan;
        DT(DT<-1e3)=nan;
        
        % make climatology
        for tt = 1:(size(years,1)+1)
            if tt == (size(years,1)+1)
                tclim = size(DT,1)-10:size(DT,1);
                dO2(tt,:,:,:) = nanmean(DO2(tclim,:,:,:));
                dT(tt,:,:,:) = nanmean(DT(tclim,:,:,:));
                par.time(tt) = nanmean(tclim);
            else
                tclim = years(tt,:);
                dO2(tt,:,:,:) = nanmean(DO2(tclim,:,:,:));
                dT(tt,:,:,:) = nanmean(DT(tclim,:,:,:));
                par.time(tt) = nanmean(tclim);
            end
        end
        
        clear DT DO2
    else
        sprintf('Error: RCP not recognized')
        keyboard
    end

% save out
out_name = sprintf('%s_dT_dO2.mat',model_name);
if do.rcp == 8
    out_path = sprintf('path_to_cmip5_rcp85/%s',out_name);
elseif do.rcp == 2
    out_path = sprintf('path_to_cmip5_rcp26/%s',out_name);
else
    sprintf('Error: RCP not recognized')
end
save(out_path,'dT','dO2','par','-v7.3');
clear dT dO2
end
