function [fut,pre,par,do] = get6extinct_func(do,model_name,upar)
%% Function to compute 3d metabolic index, global aerobic habitat, and spatial patterns of extinction/extirpation risks of marine animals from CMIP5/6 Projections
tic

% Fixed Metabolic Index parameters
Spec.doB = 0; % biomass normalized metabolic index?
Spec.doTref = 1; % use reference Temperature?
Spec.gscheme = 'po2'; % use partial pressure for O2 supply 
Spec.mfit = 0; % biomass exponent (not used)
Spec.Tref = 15; % reference Temperature deg. C

%% Set up environmental variables + measurements
% WOA grid, temperature, O2 climatologies
load('path_to_WOA_data.mat','reWOA'); % regridded to CMIP dimensions

% Load CMIP data to create reference anomalies for baseline climate (WOA)
if do.cmip == 5
    % CMIP5 files
    if do.rcp == 8
        
        % load annual climate anomalies from cmip (anomaly relative to time-mean 1850-1900)
        o2_path = sprintf('path_to_cmip5_rcp85/cmip5_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip5_rcp85/cmip5_dT_%s.nc',model_name);
        scenario = sprintf('Business as usual: RCP 8.5')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        dO2 = DO2;dO2(dO2<-1e3)=nan;
        dT = DT;dT(dT<-1e3)=nan;
        clear DT DO2
        
    elseif do.rcp == 2
        
        % load annual climate anomalies from cmip (anomaly relative to time-mean 1850-1900)
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
        
        % load annual climate anomalies from cmip (anomaly relative to time-mean 1850-1900)
        o2_path = sprintf('path_to_cmip6_ssp85/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp85/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('High emissions: SSP 8.5')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        dO2 = DO2;dO2(dO2<-1e3)=nan;
        dT = DT;dT(dT<-1e3)=nan;
        clear DT DO2
            
    
    elseif do.rcp == 4
        
        % load annual climate anomalies from cmip (anomaly relative to time-mean 1850-1900)
        o2_path = sprintf('path_to_cmip6_ssp45/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp45/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('Mid emissions: SSP 4.5')
        ncload(t_path,'DT')
        ncload(o2_path,'DO2')
        dO2 = DO2;dO2(dO2<-1e3)=nan;
        dT = DT;dT(dT<-1e3)=nan;
        clear DT DO2
        
    elseif do.rcp == 2
        
        % load annual climate anomalies from cmip (anomaly relative to time-mean 1850-1900)
        o2_path = sprintf('path_to_cmip6_ssp26/cmip6_dO2_%s.nc',model_name);
        t_path = sprintf('path_to_cmip6_ssp26/cmip6_dT_%s.nc',model_name);
        scenario = sprintf('Low emissions: SSP 2.6')
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

% Load 30-year climate anomalies
if do.clim 
    clear dT dO2
    if do.cmip == 5
        if do.rcp == 8
            % multi-decadal anomalies
            clim_path = sprintf('path_to_cmip5_rcp85/%s_dT_dO2.mat',model_name);
            load(clim_path,'dT','dO2');
        else
            % multi-decadal anomalies
            clim_path = sprintf('path_to_cmip5_rcp26/%s_dT_dO2.mat',model_name);
            load(clim_path,'dT','dO2');
        end
    elseif do.cmip == 6
        if do.rcp == 8
            % multi-decadal anomalies
            clim_path = sprintf('path_to_cmip6_ssp85/%s_dT_dO2.mat',model_name);
            load(clim_path,'dT','dO2');
        elseif do.rcp == 4
            % multi-decadal anomalies
            clim_path = sprintf('path_to_cmip6_ssp45/%s_dT_dO2.mat',model_name);
            load(clim_path,'dT','dO2');
        elseif do.rcp == 2
            % multi-decadal anomalies
            clim_path = sprintf('path_to_cmip6_ssp26/%s_dT_dO2.mat',model_name);
            load(clim_path,'dT','dO2');
        end
    end
end
    
%% Metabolic index PDFs 
getPDF;


% Tmin PDF 
par.nomin = 100; % minimum obis observations for Tmin
[par2] = getTmin2(par);
par.Tmin = par2.Tmin; 
par.Tmw = par2.Tmw;


% Species weights
Eow = par.Eow;
Aow = par.Aow;
Pcw = par.Pcw;
Tmw = par.Tmw;

%% Extinction parameters
par.Vcrit = upar.Vcrit; % critical habitat volume threshold (Vcrit,%)
par.zhab = upar.zhab ; % maximum habitat depth; 500 m (z-level 14) 
par.zmax = par.zhab; % maximum depth for spatial averaging; 

%% Compute Metabolic Index
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
pre.V = zeros(length(par.Tmin),length(par.phicrit),length(par.Eo),length(par.Ao))*nan;

% Basin masks
if do.basin
    
    % Atlantic
    ATLmask = reWOA.Bmask*nan;ATLmask(reWOA.Bmask == 2) = 1;
    ATLmask(reWOA.Bmask == 5) = 1; % includes Med. Sea
    
    % Pacific
    PACmask = reWOA.Bmask*nan;PACmask(reWOA.Bmask == 1) = 1;
    
    % Indian
    INDmask = reWOA.Bmask*nan;INDmask(reWOA.Bmask == 3) = 1;
    INDmask(reWOA.Bmask == 8) = 1; % includes Med. Sea
    
    % Arctic
    ARCmask = reWOA.Bmask*nan;ARCmask(reWOA.Bmask == 4) = 1;
   
    % initialize ocean volumes by basin
    pre.Vpac = zeros(length(par.Tmin),length(par.phicrit),length(par.Eo),length(par.Ao))*nan;
    pre.Vatl = zeros(length(par.Tmin),length(par.phicrit),length(par.Eo),length(par.Ao))*nan;
    pre.Vind = zeros(length(par.Tmin),length(par.phicrit),length(par.Eo),length(par.Ao))*nan;
    pre.Varc = zeros(length(par.Tmin),length(par.phicrit),length(par.Eo),length(par.Ao))*nan;
end
%% Pre-industrial O2 and Temp fields

% Pre-industrial Richness and Species Masks
sprintf('Computing Pre-Industrial Richness...')

% Pre-industrial O2 and Temp (WOA - WOA anomaly since pre-industrial)
pre.o2 = woa_o2-O2ref; pre.o2(pre.o2<0)=0; % prevent negative O2
pre.temp = woa_temp-Tref;pre.temp(pre.temp<-2)=-2; % freezing point of seawater @ surface (where watermasses form)

% Ecophysiotypes
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
                
                % Habitat volume 
                pre.V(pp,ii,ss,jj) = nansum(nansum(nansum(mask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
            
                % Basin-scale habitat volumes
                if do.basin
                    pre.Vatl(pp,ii,ss,jj) = nansum(nansum(nansum(ATLmask(2:par.zhab,:,:).*mask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                    pre.Vpac(pp,ii,ss,jj) = nansum(nansum(nansum(PACmask(2:par.zhab,:,:).*mask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                    pre.Vind(pp,ii,ss,jj) = nansum(nansum(nansum(INDmask(2:par.zhab,:,:).*mask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                    pre.Varc(pp,ii,ss,jj) = nansum(nansum(nansum(ARCmask(2:par.zhab,:,:).*mask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                end
            end
        end
    end
end
sprintf('done!')
toc

%% Climate change and Extinction
sprintf('The 6th Extinction...')

% sampling years
if do.clim == 0
    if do.rcp == 2 && strcmp(model_name,'IPSL')
        time = [1 39 89 139 189];% starts in year 1912
    elseif do.rcp == 2 && strcmp(model_name,'MPI')
        time = [1 39 89 139 189];% starts in year 1912
    else
        time = 0:50:size(dT,1);time(1)=1;% starts in year 1901
    end
else
    time = 1:size(dT,1);
end

for yy = 1:length(time) % CMIP year
    
    % year
    tt = time(yy);
    tt
    
    % reset extinction 3d
    extr = zeros(33,180,360);
    rextm = zeros(33,180,360);
    extm = zeros(33,180,360);
    rexts = zeros(33,180,360);
    exts = zeros(33,180,360);
    
    % Sum CMIP anomalies to reference climatology
    o2 = squeeze(dO2(tt,:,:,:)) + pre.o2;o2(o2<0)=0;
    temp = squeeze(dT(tt,:,:,:)) + pre.temp;temp(temp<-2)=-2;
    
    % Ecophysiotypes
    for pp = 1:length(par.Tmin)
        for ii = 1:length(par.phicrit)
            for ss = 1:length(par.Eo)
                for jj = 1:length(par.Ao)
                    
                    % traits
                    Tmin = par.Tmin(pp);
                    phi_crit = par.phicrit(ii);
                    Spec.Eo = par.Eo(ss);
                    Spec.Ao = par.Ao(jj);
                    
                    % Pre-industrial Metabolic Index
                    phi=metab_index(Spec,pre.temp,pre.o2,reWOA.z3d);
                    
                    % Pre-industrial Presence/Absence Mask
                    Pmask = phi*0;Pmask(phi>=phi_crit & pre.temp>=Tmin)=1;
                    
                    % Future Metabolic index
                    phi=metab_index(Spec,temp,o2,reWOA.z3d);
                    
                    % Future Presence/Absence Mask
                    Fmask = phi*0;Fmask(phi>=phi_crit & temp>=Tmin)=1;
                    
                    % Extirpation
                    extr = extr + (Pmask-Fmask).*((Pmask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj))./pre.rich);
               
                    % Habitat volume [mobile, sessile]
                    fut.Vm(pp,ii,ss,jj) = nansum(nansum(nansum(Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                    fut.Vs(pp,ii,ss,jj) = nansum(nansum(nansum(Pmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                    
                    % Basin-scale volume changes
                    if do.basin
                        % Future habitat volume [mobile, sessile]
                        % Atlantic - habitat volumes
                        fut.Vm_atl(pp,ii,ss,jj) = nansum(nansum(nansum(ATLmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                        fut.Vs_atl(pp,ii,ss,jj) = nansum(nansum(nansum(ATLmask(2:par.zhab,:,:).*Pmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                        
                        % Atlantic - Mobile habitat Loss (% change)
                        atl_delVm = -100*(fut.Vm_atl(pp,ii,ss,jj)-pre.Vatl(pp,ii,ss,jj))./pre.Vatl(pp,ii,ss,jj);
                        atl_delVm(atl_delVm<-100)=-100;
                        atl_Habloss_m(pp,ii,ss,jj) = atl_delVm;
                        
                        % Atlantic - Sessile habitat Loss (% change)
                        atl_delVs = -100*(fut.Vs_atl(pp,ii,ss,jj)-pre.Vatl(pp,ii,ss,jj))./pre.Vatl(pp,ii,ss,jj);
                        atl_delVs(atl_delVs<-100)=-100;
                        atl_Habloss_s(pp,ii,ss,jj) = atl_delVs;
                        
                        % Pacific - habitat volumes
                        fut.Vm_pac(pp,ii,ss,jj) = nansum(nansum(nansum(PACmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                        fut.Vs_pac(pp,ii,ss,jj) = nansum(nansum(nansum(PACmask(2:par.zhab,:,:).*Pmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                        
                        % Pacific - Mobile habitat Loss (% change)
                        pac_delVm = -100*(fut.Vm_pac(pp,ii,ss,jj)-pre.Vpac(pp,ii,ss,jj))./pre.Vpac(pp,ii,ss,jj);
                        pac_delVm(pac_delVm<-100)=-100;
                        pac_Habloss_m(pp,ii,ss,jj) = pac_delVm;
                        
                        % Pacific - Sessile habitat Loss (% change)
                        pac_delVs = -100*(fut.Vs_pac(pp,ii,ss,jj)-pre.Vpac(pp,ii,ss,jj))./pre.Vpac(pp,ii,ss,jj);
                        pac_delVs(pac_delVs<-100)=-100;
                        pac_Habloss_s(pp,ii,ss,jj) = pac_delVs;
                        
                        % Indian
                        fut.Vm_ind(pp,ii,ss,jj) = nansum(nansum(nansum(INDmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                        fut.Vs_ind(pp,ii,ss,jj) = nansum(nansum(nansum(INDmask(2:par.zhab,:,:).*Pmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                        
                        % Indian - Mobile habitat Loss (% change)
                        ind_delVm = -100*(fut.Vm_ind(pp,ii,ss,jj)-pre.Vind(pp,ii,ss,jj))./pre.Vind(pp,ii,ss,jj);
                        ind_delVm(ind_delVm<-100)=-100;
                        ind_Habloss_m(pp,ii,ss,jj) = ind_delVm;
                        
                        % Indian - Sessile habitat Loss (% change)
                        ind_delVs = -100*(fut.Vs_ind(pp,ii,ss,jj)-pre.Vind(pp,ii,ss,jj))./pre.Vind(pp,ii,ss,jj);
                        ind_delVs(ind_delVs<-100)=-100;
                        ind_Habloss_s(pp,ii,ss,jj) = ind_delVs;
                        
                        % Arctic
                        fut.Vm_arc(pp,ii,ss,jj) = nansum(nansum(nansum(ARCmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                        fut.Vs_arc(pp,ii,ss,jj) = nansum(nansum(nansum(ARCmask(2:par.zhab,:,:).*Pmask(2:par.zhab,:,:).*Fmask(2:par.zhab,:,:).*reWOA.v3d(2:par.zhab,:,:))));
                        
                        % Arctic - Mobile habitat Loss (% change)
                        arc_delVm = -100*(fut.Vm_arc(pp,ii,ss,jj)-pre.Varc(pp,ii,ss,jj))./pre.Varc(pp,ii,ss,jj);
                        arc_delVm(arc_delVm<-100)=-100;
                        arc_Habloss_m(pp,ii,ss,jj) = arc_delVm;
                        
                        % Arctic - Sessile habitat Loss (% change)
                        arc_delVs = -100*(fut.Vs_arc(pp,ii,ss,jj)-pre.Varc(pp,ii,ss,jj))./pre.Varc(pp,ii,ss,jj);
                        arc_delVs(arc_delVs<-100)=-100;
                        arc_Habloss_s(pp,ii,ss,jj) = arc_delVs;
                    end
                    
                    %% Mobile habitat Loss (% change)
                    delVm = -100*(fut.Vm(pp,ii,ss,jj)-pre.V(pp,ii,ss,jj))./pre.V(pp,ii,ss,jj);
                    delVm(delVm<-100)=-100; % set maximum habitat gain to 100% (this has no effect on results)
                    Habloss_m(pp,ii,ss,jj) = delVm;
                    
                    % Global extinction (mobile species)
                    if delVm < par.Vcrit % Not extinct
                    
                        Rmaskm = phi*0;Rmaskm(phi>=phi_crit & temp>=Tmin)=1;% local extirpation
                        Emaskm = Pmask; % survivors
                    
                    elseif delVm >= par.Vcrit % Extinct
                        
                        Rmaskm = Pmask*0;% total extirpation
                        Emaskm = Pmask*0;% extinction
                    
                    else % zero habitat species (no extinction)
                        Rmaskm = phi*0;Rmaskm(phi>=phi_crit & temp>=Tmin)=1;% local extirpation
                        Emaskm = Pmask; % not extinct
                    end
                    
                    % Regional Extinction (mobile)
                    rextm = rextm + (Pmask-Rmaskm).*((Pmask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj))./pre.rich);
                    
                    % Global Extinction (mobile)
                    extm = extm + (Pmask-Emaskm).*((Pmask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj))./pre.rich);
                    
                    
                    %%  Sessile habitat Loss (% change)
                    delVs = -100*(fut.Vs(pp,ii,ss,jj)-pre.V(pp,ii,ss,jj))./pre.V(pp,ii,ss,jj);
                    delVs(delVs<-100)=-100;% set maximum habitat gain to 100% (this has no effect on results)
                    Habloss_s(pp,ii,ss,jj) = delVs;
                    
                    % Global extinction (sessile species)
                    if delVs < par.Vcrit % Not extinct
                        
                        Rmasks = phi*0;Rmasks(phi>=phi_crit & temp>=Tmin)=1;% local extirpation
                        Emasks = Pmask; % survival
                    
                    elseif delVs >= par.Vcrit % Extinct
                        
                        Rmasks = Pmask*0;% total extirpation
                        Emasks = Pmask*0;% extinction
                    
                    else % zero habitat species (no extinction)
                        Rmasks = phi*0;Rmasks(phi>=phi_crit & temp>=Tmin)=1;% local extirpation
                        Emasks = Pmask; % not extinct
                    end
                    
                    % Regional Extinction (sessile)
                    rexts = rexts + (Pmask-Rmasks).*((Pmask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj))./pre.rich);
                    
                    % Global Extinction (sessile)
                    exts = exts + (Pmask-Emasks).*((Pmask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj))./pre.rich);
                    
                end
            end
        end
    end
    % Average over Z,Y,X
    idx = 2:par.zmax;
    tocn = temp*nan;tocn(isnan(temp)==0)=1;
    eocn = extr*nan; eocn(isnan(extr)==0)=1;
    
    % global means % 
    % future temperature and o2
    fut.Tu(yy) = squeeze(nansum(nansum(nansum(temp(idx,:,:).*reWOA.v3d(idx,:,:).*tocn(idx,:,:))))./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*tocn(idx,:,:)))));
    fut.Ou(yy) = squeeze(nansum(nansum(nansum(o2(idx,:,:).*reWOA.v3d(idx,:,:).*tocn(idx,:,:))))./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*tocn(idx,:,:)))));
    
    % pre-industrial temperature and o2
    pre.Tu(yy) = squeeze(nansum(nansum(nansum(pre.temp(idx,:,:).*reWOA.v3d(idx,:,:).*tocn(idx,:,:))))./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*tocn(idx,:,:)))));
    pre.Ou(yy) = squeeze(nansum(nansum(nansum(pre.o2(idx,:,:).*reWOA.v3d(idx,:,:).*tocn(idx,:,:))))./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*tocn(idx,:,:)))));
    
    % extirpation and richness
    fut.extru(yy) = 100*squeeze(nansum(nansum(nansum(extr(idx,:,:).*reWOA.v3d(idx,:,:).*eocn(idx,:,:),1),2),3)./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*eocn(idx,:,:)))));
    pre.richu(yy) = squeeze(nansum(nansum(nansum(pre.rich(idx,:,:).*reWOA.v3d(idx,:,:).*tocn(idx,:,:),1),2),3)./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*tocn(idx,:,:)))));
    
    % regional,global extinction (mobile)
    fut.rextmu(yy) = 100*squeeze(nansum(nansum(nansum(rextm(idx,:,:).*reWOA.v3d(idx,:,:).*eocn(idx,:,:),1),2),3)./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*eocn(idx,:,:)))));
    fut.extmu(yy) = 100*squeeze(nansum(nansum(nansum(extm(idx,:,:).*reWOA.v3d(idx,:,:).*eocn(idx,:,:),1),2),3)./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*eocn(idx,:,:)))));
    
    % regional,global extinction (sessile)
    fut.rextsu(yy) = 100*squeeze(nansum(nansum(nansum(rexts(idx,:,:).*reWOA.v3d(idx,:,:).*eocn(idx,:,:),1),2),3)./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*eocn(idx,:,:)))));
    fut.extsu(yy) = 100*squeeze(nansum(nansum(nansum(exts(idx,:,:).*reWOA.v3d(idx,:,:).*eocn(idx,:,:),1),2),3)./nansum(nansum(nansum(reWOA.v3d(idx,:,:).*eocn(idx,:,:)))));
    
    % 3d fields 
    fut.extr(yy,:,:,:) = extr; % extirpation
    fut.rextm(yy,:,:,:) = rextm; % regional extinction (mobile) 
    fut.extm(yy,:,:,:) = extm; % global extinction (mobile)
    fut.rexts(yy,:,:,:) = rexts; % regional extinction (sessile)
    fut.exts(yy,:,:,:) = exts; % global extinction (sessile)
    
    % Habitat loss (mobile)
    fut.delVm(yy,:,:,:,:) = Habloss_m; 
    
    % Habitat loss (sessile)
    fut.delVs(yy,:,:,:,:) = Habloss_s; 
    
    % Basin-scale habitat loss volumes
    if do.basin
        % Atlantic habitat loss 
        fut.delVm_atl(yy,:,:,:,:) = atl_Habloss_m;
        fut.delVs_atl(yy,:,:,:,:) = atl_Habloss_s;
        
        % Pacific habitat loss 
        fut.delVm_pac(yy,:,:,:,:) = pac_Habloss_m;
        fut.delVs_pac(yy,:,:,:,:) = pac_Habloss_s;
        
        % Indian habitat loss 
        fut.delVm_ind(yy,:,:,:,:) = ind_Habloss_m;
        fut.delVs_ind(yy,:,:,:,:) = ind_Habloss_s;
        
        % Arctic habitat loss 
        fut.delVm_arc(yy,:,:,:,:) = arc_Habloss_m;
        fut.delVs_arc(yy,:,:,:,:) = arc_Habloss_s;
    end
    
    
    % print out extirpation
    fprintf('Extirpation = %g\n',fut.extru(yy));
    
end

% time
% store temporary parameters
temp_par = par;
load(clim_path,'par');
if do.cmip == 5
    if do.rcp == 2 && strcmp(model_name,'IPSL')
        temp_par.time=par.time+1911;
    elseif do.rcp == 2 && strcmp(model_name,'MPI')
        temp_par.time=par.time+1911;
    else
        temp_par.time=par.time+1900;
    end
else
    temp_par.time=par.time+1900;
end
par = temp_par;

par.runtime = toc;   
