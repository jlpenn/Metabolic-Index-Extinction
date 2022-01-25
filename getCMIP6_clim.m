%% Make multi-decadal climatologies of CMIP6 temperature and oxygen anomalies
clear all;close all

tic

% RCP
do.rcp = 2; % rcp 8.5 (8), rcp 4.5 (4) rcp 2.6 (2)

% CMIP models
if do.rcp == 8
    models ={'CanESM';'CSIRO';'GFDL';'IPSL';'MPI';'CNRM';'UKESM';'NorESM';'MIROC';'MRI'};
else
    models ={'CanESM';'CSIRO';'GFDL';'IPSL';'MPI';'CNRM';'UKESM';'NorESM';'MIROC'};
end

for mm = 1:length(models)
    
    % model
    model_name = models{mm}
    
    % sampling years
    if do.rcp == 8 || do.rcp == 2
        if strcmp(model_name,'CSIRO') || strcmp(model_name,'GFDL') || strcmp(model_name,'MPI') || strcmp(model_name,'CNRM') || strcmp(model_name,'UKESM') || strcmp(model_name,'NorESM') || strcmp(model_name,'MIROC') || strcmp(model_name,'MRI')
            year0 = 1900; % start year
            years = [2005:2035; 2035:2065; 2070:2100] - year0;
        else
            year0 = 1900; % start year
            years = [2005:2035; 2035:2065; 2085:2115; 2135:2165; 2185:2215; 2235:2265; 2270:2300] - year0;
        end
    else % SSP4.5
         year0 = 1900; % start year
         years = [2005:2035; 2035:2065; 2070:2100] - year0;
    end
    %% Set up environmental variables + measurements
    
    % CMIP5 files
    if do.rcp == 8
        
        % load annual climate anomalies from cmip (anomaly is relative to time-mean 1850-1900; dimensions: time,depth,latitude,longitude)
        o2_path = sprintf('path_to_cmip6_ssp85/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp85/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('High emissions: SSP8.5')
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
        
    elseif do.rcp == 4
        
        % load annual climate anomalies from cmip (anomaly is relative to time-mean 1850-1900; dimensions: time,depth,latitude,longitude)
        o2_path = sprintf('path_to_cmip6_ssp45/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp45/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('Mid emissions: SSP 4.5')
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
        
        % load annual climate anomalies from cmip (anomaly is relative to time-mean 1850-1900)
        o2_path = sprintf('path_to_cmip6_ssp26/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp26/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('Low emissions: SSP 2.6')
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
    out_path = sprintf('path_to_cmip6_ssp85/%s',out_name);
elseif do.rcp == 4
    out_path = sprintf('path_to_cmip6_ssp45/%s',out_name);
elseif do.rcp == 2
    out_path = sprintf('path_to_cmip6_ssp26/%s',out_name);
else
    sprintf('Error: RCP not recognized')
end
save(out_path,'dT','dO2','par','-v7.3');
toc 
clear dT dO2 par
end
