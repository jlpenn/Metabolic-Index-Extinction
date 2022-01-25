%% Compute 3d metabolic index, global aerobic habitat, and spatial patterns of extinction/extirpation risks of marine animals from CMIP5/6 Projections
clear all;close all

% CMIP5 vs. CMIP6
do.cmip = 6;

% anomaly fields from 30-yr climatologies (1) vs single years (0)?
do.clim = 1;

% Remove simulated historical climate trends from WOA data? (yes = 1; no = 0)
do.ref = 1;

% Bianchi correct?
do.bianchi = 1;

% Temperature vs O2 components 
do.const = 0; % Fix temperature or O2 at a constant value (1) in the future or vary both (0)?

% include black/caspian seas? (1 = include; 0 = exclude)
do.sea = 0;

% extinction by basin? (compute extinction globally (0) or by ocean basin (1)
do.basin = 0;

% Emissions scenario (2 or 8) 
rcp = 8;

% global extinction parameters
upar.Vcrit = 50; % critical habitat loss threshold (Vcrit; percent)
upar.zhab = 14; % Maximum habitat depth (500 m) (z level 14);

% CMIP models
if do.cmip == 5
    models = {'HADL';'MPI';'IPSL';'CESM';'GFDL'};
elseif do.cmip == 6 && rcp == 8
    models = {'CanESM';'CSIRO';'GFDL';'IPSL';'MPI';'CNRM';'UKESM';'NorESM';'MIROC';'MRI'};
elseif do.cmip == 6 && rcp == 2
    models = {'CanESM';'CSIRO';'GFDL';'IPSL';'MPI';'CNRM';'UKESM';'NorESM';'MIROC'};
end

for rr = 1:length(rcp)
    for mm = 1:length(models)
        do.rcp = rcp(rr)
        model_name = models{mm}
        
        % compute extinction
        [fut,pre,par,do] = get6extinct_func(do,model_name,upar);
        
        % save out
        out_name = sprintf('Extinct_%s_WOAref_Vcrit50.mat',model_name);
        
        if do.cmip == 5
            if do.rcp == 8
                out_path = sprintf('output/CMIP5/rcp85/%s',out_name);
            elseif do.rcp == 2
                out_path = sprintf('output/CMIP5/rcp26/%s',out_name);
            else
                sprintf('Error: RCP not recognized')
            end
        else
            if do.rcp == 8
                out_path = sprintf('output/CMIP6/ssp85/%s',out_name);
            elseif do.rcp == 4
                out_path = sprintf('output/CMIP6/ssp45/%s',out_name);
            elseif do.rcp == 2
                out_path = sprintf('output/CMIP6/ssp26/%s',out_name);
            else
                sprintf('Error: RCP not recognized')
            end
        end 
        save(out_path,'pre','fut','par','do','-v7.3');
        sprintf('%s The sixth extinction saved...',model_name)
        
        clear fut pre par 
        
    end
end

sprintf('It would not be much of a universe if the ones you loved were not in it. -Stephen Hawking')
