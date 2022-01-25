%% Compute richness of animal species versus latitude from WOA data
clear all;close all

% CMIP5 vs. CMIP6
do.cmip = 6;

% Remove simulated historical climate anomalies from WOA data? (1 = yes; 0 = no)
do.ref = 0;

% Bianchi correct?
do.bianchi = 1;

% Temperature vs O2 components 
do.const = 0; % Fix Temp,O2 to constant value?

% include black/caspian seas? (1 = include; 0 = exclude)
do.sea = 0;

% extinction by basin?
do.basin = 0;

% Latitude and depth intervals for G-div 
upar.yres = 5;% latitude bin resolution (degrees) 
upar.zres = 2:14;% depth interval (5-500 m)

% RCP 
rcp = 8;

% CMIP models
if do.cmip == 5
    models = {'HADL';'MPI';'IPSL';'CESM';'GFDL'};
elseif do.cmip == 6 && rcp == 8
    models = {'GFDL'};%{'CanESM';'CSIRO';'GFDL';'IPSL';'MPI';'CNRM';'UKESM';'NorESM';'MIROC';'MRI'};
elseif do.cmip == 6 && rcp == 2
    models = {'CanESM';'CSIRO';'GFDL';'IPSL';'MPI';'CNRM';'UKESM';'NorESM';'MIROC'};
end

for rr = 1:length(rcp)
    for mm = 1:length(models)
        do.rcp = rcp(rr)
        model_name = models{mm}
        
        % compute extinction
        [pre,par,do] = getRichness_func(do,model_name,upar);
        
        % save out
        out_name = sprintf(Richness_WOA_500m.mat',model_name);
        
        if do.cmip == 5
            if do.rcp == 8
                out_path = sprintf('output/Richness/%s',out_name);
            elseif do.rcp == 2
                 out_path = sprintf('output/Richness/%s',out_name);
            else
                sprintf('Error: RCP not recognized')
            end
        else
            if do.rcp == 8
                out_path = sprintf('output/Richness/%s',out_name);
            elseif do.rcp == 4
                out_path = sprintf('output/Richness/%s',out_name);
            elseif do.rcp == 2
                out_path = sprintf('output/Richness/%s',out_name);
            else
                sprintf('Error: RCP not recognized')
            end
        end 
        save(out_path,'pre','par','do','-v7.3');
        sprintf('% Richness vs. latitude computed...',model_name)
        
        clear pre par
        
       
    end
end

sprintf('It would not be much of a universe if the ones you loved were not in it. -Stephen Hawking')
