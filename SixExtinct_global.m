%% Global extinction as a fraction of the entire global community from CMIP5/6 Projections 
clear all;close all

models = {'CanESM';'CSIRO';'GFDL';'IPSL';'MPI';'CNRM';'UKESM';'NorESM';'MIROC';'MRI'};
for mm = 1:length(models)
    
model_name = models{mm};
out_name = sprintf('Extinct_%s_WOAref_Vcrit70.mat',model_name); % load global habitat volume changes for extinction
out_path = sprintf('output/CMIP6/ssp85/%s',out_name);
load(out_path)

% Extinction threshold
par.Vcrit = 70;

% simplify notation
Tmw = par.Tmw;
Pcw = par.Pcw;
Eow = par.Eow;
Aow = par.Aow;

%% Pre-industrial richness
% Ecophyiotypes
pre.richi = 0;
for pp = 1:length(par.Tmin)
    for ii = 1:length(par.phicrit)
        for ss =1:length(par.Eo)
            for jj = 1:length(par.Ao)
                
                % Pre-industrial volume? 
                pmask = pre.V(pp,ii,ss,jj)>0; %(1 = presence; 0 = absence)
                
                % Pre-industrial Richness
                pre.richi = pre.richi+(pmask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj));% sum number of species with pre-industrial habitat
               
            end
        end
    end
end
sprintf('done!')

%% Aggregated Extinction 
fut.extsi = zeros(1,length(fut.Tu));
fut.extmi = zeros(1,length(fut.Tu));

for tt = 1:length(fut.Tu)
for pp = 1:length(par.Tmin)
    for ii = 1:length(par.phicrit)
        for ss =1:length(par.Eo)
            for jj = 1:length(par.Ao)
                
                % Pre-industrial volume?
                pmask = pre.V(pp,ii,ss,jj)>0; %(1 = presence; 0 = absence)
                
                % Extinction: Habitat loss only 
                smask = fut.delVs(tt,pp,ii,ss,jj)>=par.Vcrit; %(1 = extinct; 0 = survivor)
                fut.extsi(tt) = fut.extsi(tt)+((smask.*pmask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj))./pre.richi); 
                
                % Extinction: Habitat loss + gain
                mmask = fut.delVm(tt,pp,ii,ss,jj)>=par.Vcrit; %(1 = extinct; 0 = survivor)
                fut.extmi(tt) = fut.extmi(tt)+((mmask.*pmask.*Tmw(pp).*Pcw(ii).*Eow(ss).*Aow(jj))./pre.richi);
                
            end
        end
    end
end
end
sprintf('done!')

% save results 
output.richi = pre.richi;
output.extmi = 100*fut.extmi;
output.extsi = 100*fut.extsi;
output.warm = fut.Tu-pre.Tu;

out_name = sprintf('Extinct_%s_WOAref_Vcrit70_global.mat',model_name);
out_path = sprintf('output/CMIP6/ssp85/%s',out_name);
save(out_path,'output','par','-v7.3');
clear output par 
end 
