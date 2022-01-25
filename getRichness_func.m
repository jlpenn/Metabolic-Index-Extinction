function [pre,par,do] = getRichness_func(do,model_name,upar)
%% Function to compute richness of species versus latitude from WOA data 
 
tic

% Fixed Metabolic Index parameters
Spec.doB = 0; % biomass normalized index?
Spec.doTref = 1; % use reference Temperature?
Spec.gscheme = 'po2'; % partial pressure
Spec.mfit = 0; % biomass exponent
Spec.Tref = 15; % reference Temperature oC

%% Set up environmental variables + measurements
% WOA grid, temperature, O2 climatologies
load('path_to_WOA_data.mat','reWOA'); % regridded to CMIP dimensions

if do.cmip == 5
    % CMIP5 files
    if do.rcp == 8
        
        % load climate anomalies
        o2_path = sprintf('path_to_cmip5_rcp85/cmip5_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip5_rcp85/cmip5_dT_%s.nc',model_name);
        scenario = sprintf('Business as usual: RCP 8.5')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        dO2 = DO2;dO2(dO2<-1e3)=nan;
        dT = DT;dT(dT<-1e3)=nan;
        clear DT DO2
        
    elseif do.rcp == 2
        
        % load climate anomalies
        o2_path = sprintf('path_to_cmip5_rcp26/cmip5_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip5_rcp26/cmip5_dT_%s.nc',model_name);
        scenario = sprintf('Best case: RCP 2.6')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        dO2 = DO2;dO2(dO2<-1e3)=nan;
        dT = DT;dT(dT<-1e3)=nan;
        clear DT DO2
    else
        sprintf('Error: RCP not recognized')
        keyboard
    end
elseif do.cmip == 6
    % CMIP6 files
    if do.rcp == 8
        
        % load climate anomalies
        o2_path = sprintf('path_to_cmip6_ssp85/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp85/ssp585/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('High emissions: RCP 8.5')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        dO2 = DO2;dO2(dO2<-1e3)=nan;
        dT = DT;dT(dT<-1e3)=nan;
        clear DT DO2
            
    
    elseif do.rcp == 4
        
        % load climate anomalies
        o2_path = sprintf('path_to_cmip6_ssp45/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp45/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('Mid emissions: RCP 4.5')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        dO2 = DO2;dO2(dO2<-1e3)=nan;
        dT = DT;dT(dT<-1e3)=nan;
        clear DT DO2
        
    elseif do.rcp == 2
        
        % load climate anomalies
        o2_path = sprintf('path_to_cmip6_ssp26/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp26/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('Low emissions: RCP 2.6')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        dO2 = DO2;dO2(dO2<-1e3)=nan;
        dT = DT;dT(dT<-1e3)=nan;
        clear DT DO2
    else
        sprintf('Error: RCP not recognized')
        keyboard
    end
    
end

% Compute reference anomalies for WOA years (1955 to 2014)
if do.ref == 1
    O2ref = squeeze(nanmean(dO2(55:114,:,:,:),1));
    Tref = squeeze(nanmean(dT(55:114,:,:,:),1));
    
    % make same land/ocn mask for T and O2
    if strcmp(model_name,'CSIRO') || strcmp(model_name,'MRI')
        mask1 = O2ref.*nan;mask1(~isnan(O2ref))=1;
        mask2 = Tref.*nan;mask2(~isnan(Tref))=1;
        mask3 = mask1.*mask2;
        O2ref = O2ref.*mask3;
        Tref = Tref.*mask3;
        clear mask1 mask2 mask3
    end
else
    O2ref = 0;
    Tref = 0;
end
    
%% Metabolic index PDFs 
getPDF;

% Tmin PDF from OBIS minimum temperature data
par.nomin = 100; % minimum obis observations for Tmin
[par2] = getTmin2(par);
par.Tmin = par2.Tmin;
par.Tmw = par2.Tmw;


% Species weights
Eow = par.Eow;
Aow = par.Aow;
Pcw = par.Pcw;
Tmw = par.Tmw;


%% Compute Index
% Annual MI
% WOA annual average
if do.bianchi == 1
    reWOA.o2 = 1.009*reWOA.o2 - 2.523; 
    reWOA.o2(reWOA.o2<0)=0;
end
woa_o2 = squeeze(nanmean(reWOA.o2,1));
woa_temp = squeeze(nanmean(reWOA.temp,1));

% remove unconnected seas
if do.sea < 1
    Bmask = zeros(33,180,360)+1;
    
    % mask out Black Sea
    Bmask(reWOA.Bmask==6)=NaN;
    
    % mask out Caspian Sea
    Bmask(reWOA.Bmask==7) = NaN;
    
    % remove seas
    woa_o2 = woa_o2.*Bmask;
    woa_temp = woa_temp.*Bmask;
    
end
clear reWOA.o2 reWOA.temp

% pre-allocate space for speed
pre.rich = zeros(33,180,360);
pre.rich_lat = zeros(1,36);

%% Pre-industrial O2 and Temp fields -> Gamma diversity

% Pre-industrial Richness and Species Masks
sprintf('Computing Pre-Industrial Gamma-diversity (richness)...')

% Pre-industrial O2 and Temp (WOA - WOA anomaly since pre-industrial)
pre.o2 = woa_o2-O2ref; pre.o2(pre.o2<0)=0; % prevent negative O2
pre.temp = woa_temp-Tref;pre.temp(pre.temp<-2)=-2; % freezing point of seawater @ surface (where watermasses form)

% Ecophyiotypes
for pp = 1:length(par.Tmin)
    for ii = 1:length(par.phicrit)
        for ss =1:length(par.Eo)
            for jj = 1:length(par.Ao)
                
                % Traits
                Tmin = par.Tmin(pp); 
                phi_crit = par.phicrit(ii);
                Spec.Eo = par.Eo(ss);
                Spec.Ao = par.Ao(jj);
                
                % Metabolic index
                phi=metab_index(Spec,pre.temp,pre.o2,reWOA.z3d);
                
                % Presence/Absence Mask
                mask = phi*0;mask(phi>=phi_crit & pre.temp>=Tmin)=1;
                
                % Pre-industrial Richness
                pre.rich = pre.rich+(mask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj));
                
                % latitudes
                lat = -90:upar.yres:90;% binned latitude
                
                % sum species over latitude
                for yy = 1:length(lat)-1
                    
                    if yy == 1
                        ymin = 1;
                        ymax = 5;
                    else
                        ymin = ymax + 1;
                        ymax = ymax + 5;
                    end
                    
                    % species presence in latitude band?
                    mask_l = nansum(nansum(nansum(mask(upar.zres,ymin:ymax,:))));
                    
                    % make binary
                    mask_l = mask_l>0;
                    
                    % sum gamma diversity
                    pre.rich_lat(yy) = pre.rich_lat(yy)+(mask_l.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj));
                    
                end
                
            end
        end
    end
end

sprintf('done!')
toc

% save latitude grid and time
par.lat = lat;
par.latm = -87.5:5:87.5;
par.runtime = toc;   
