% Code to read the IUCN red list endangered species regional assessments and compute the fraction threatened with extinction
% env: including marine, kingdom: animalia 
clear;close all;

% Exclusion switches
do.ectotherm = 1; % restrict to water-breathers?
do.marter = 1; % include marine/terr. species?
do.genus_ext = 1; % compute extinction of genera?
do.three_phyla = 0; % restrict to well-sampled Met. Index phyla (Arthropoda, Chordata, Mollusca) 
do.scope = 0; % remove global assessments 

regions = {'Arabian_Sea';'Caribean';'Gulf_of_Mexico';'Mediterranean';'Persian_Gulf';'C_Africa';'E_Africa';'Europe';'N_Africa';'NE_Africa';'S_Africa';'W_Africa'} % 'Pan_Africa';
for rr = 1:length(regions)
this_region = regions{rr}

if rr < 6
    % status data
    status_path = sprintf('path_to_regional_assessments_%s.mat',this_region);
    load(status_path);
    
    % taxonomy data
    taxonomy_path = sprintf('path_to_regional_taxonomy_%s.mat',this_region);
    load(taxonomy_path);
else
    % status data
    status_path = sprintf('path_to_regional_assessments_%s.mat.mat',this_region);
    load(status_path);

    % taxonomy data
    taxonomy_path = sprintf('path_to_regional_taxonomy_%s.mat',this_region);
    load(taxonomy_path);
end

% parse variables
% taxonomy variables 
species2 = string(table2array(taxonomy(:,2)));
kingdom = string(table2array(taxonomy(:,3)));
phylum = string(table2array(taxonomy(:,4)));
order = string(table2array(taxonomy(:,5)));
class = string(table2array(taxonomy(:,6)));
family = string(table2array(taxonomy(:,7)));
genus = string(table2array(taxonomy(:,8)));

% assessment variables
species = string(table2array(assessments(:,3)));
status = string(table2array(assessments(:,4)));
env = string(table2array(assessments(:,17)));
scope = string(table2array(assessments(:,end)));

% remove table
clear assessments taxonomy

%% clean data

% convert statuses from prior to 2001 to latest version
status(strcmp(status,'Lower Risk/near threatened')) = 'Near Threatened';% lower risk/near-threatened to NT
status(strcmp(status,'Lower Risk/conservation dependent')) = 'Near Threatened';% lower risk/conservation-depenendent to NT
status(strcmp(status,'Lower Risk/least concern')) = 'Least Concern';% lower risk/least-concern to LC
status(strcmp(status,'Not Applicable')) = 'Data Deficient';

% Environment selection
if do.marter == 0
    % marine environments
    mar1 = strcmp(env,'Marine');
    mar2 = strcmp(env,'Freshwater (=Inland waters)|Marine');
    % combine marine environments
    mar = mar1+mar2;
else % include terrestrial+marine
    mar1 = strcmp(env,'Marine');
    mar2 = strcmp(env,'Freshwater (=Inland waters)|Marine');
    mar3 = strcmp(env,'Terrestrial|Freshwater (=Inland waters)|Marine');
    mar4 = strcmp(env,'Terrestrial|Marine');
    
    % combine marine environments
    mar = mar1+mar2+mar3+mar4;
end
mar=mar>0;

% restrict to water-breathers?
if do.ectotherm 
    % restrict to animals
    animal = strcmp(kingdom,'ANIMALIA') & ~strcmp(class,'REPTILIA') & ~strcmp(class,'AMPHIBIA') & ~strcmp(class,'AVES') &  ~strcmp(class,'INSECTA') &  ~strcmp(class,'MAMMALIA') &  ~strcmp(class,'DIPLOPODA');
    
    % restrict to marine animals
    mar_animal = mar.*animal;
    mar_animal = mar_animal>0;

else
    % restrict to animals
    animal = strcmp(kingdom,'ANIMALIA');
    
    % restrict to marine animals
    mar_animal = mar.*animal;
    mar_animal = mar_animal>0;
    
end



% restrict to well-sampled Met. Index phyla?
if do.three_phyla
    % restrict to arthropods, Chordates, molluscs
    three_phyla = strcmp(phylum,'CHORDATA') | strcmp(phylum,'ARTHROPODA') | strcmp(phylum,'MOLLUSCA');
    
    % restrict to marine animals
    mar_animal = mar_animal.*three_phyla;
    mar_animal = mar_animal>0;
end

% restrict to regional assessments
if do.scope
    for idx = 1:length(scope)
        reg_scope(idx) = sum(strcmp(split(scope(idx)),'Global')) + sum(strcmp(split(scope(idx)),'Global,'));
    end
    reg_mask = reg_scope'*0;reg_mask(reg_scope'<1)=1;
    mar_animal = mar_animal.*(reg_mask);
    mar_animal = mar_animal>0;
    clear reg_scope reg_mask
end

% save marine species 
mar_species = species(mar_animal);
mar_class = class(mar_animal);
mar_status = status(mar_animal);
reg_idx = find(~strcmp(mar_status,'Not Evaluated') & ~strcmp(mar_status,'Data Deficient'));
output.mar_spec(rr,1:length(reg_idx)) = mar_species(reg_idx);
output.mar_class(rr,1:length(reg_idx)) = mar_class(reg_idx);

%% Percent species endangerment

% total marine animal species per category
output.ne(rr) = nansum(strcmp(status(mar_animal),'Not Evaluated'));
output.dd (rr)= nansum(strcmp(status(mar_animal),'Data Deficient'));
output.lc(rr) = nansum(strcmp(status(mar_animal),'Least Concern'));
output.nt(rr) = nansum(strcmp(status(mar_animal),'Near Threatened'));
output.vu(rr) = nansum(strcmp(status(mar_animal),'Vulnerable'));
output.en(rr) = nansum(strcmp(status(mar_animal),'Endangered'));
output.cr(rr) = nansum(strcmp(status(mar_animal),'Critically Endangered'));
output.ew(rr) = nansum(strcmp(status(mar_animal),'Extinct in the Wild'));
output.ex(rr) = nansum(strcmp(status(mar_animal),'Extinct'));
output.rex(rr) = nansum(strcmp(status(mar_animal),'Regionally Extinct'));

% total species across categories
output.total_dd(rr) = output.dd(rr)+output.lc(rr)+output.nt(rr)+output.vu(rr)+output.en(rr)+output.cr(rr)+output.ew(rr)+output.ex(rr)+output.rex(rr);

% total species - DD 
output.total(rr) = output.lc(rr)+output.nt(rr)+output.vu(rr)+output.en(rr)+output.cr(rr)+output.ew(rr)+output.ex(rr)+output.rex(rr);

% threatened species (VU+)
output.threat_vu(rr) = output.vu(rr)+output.en(rr)+output.cr(rr)+output.ew(rr)+output.ex(rr)+output.rex(rr);

% threatened species (NT+) 
output.threat_nt(rr) = output.nt(rr)+output.vu(rr)+output.en(rr)+output.cr(rr)+output.ew(rr)+output.ex(rr)+output.rex(rr);

% percent endangered (low-estimate)
output.endanger_low(rr) = 100*output.threat_vu(rr)./output.total(rr);

% percent endangered (high-estimate)
output.endanger_high(rr) = 100*output.threat_nt(rr)./output.total(rr);


%% Percent genus endangerment 
if do.genus_ext
    
    % initialize genus counts
    output_g.genera(rr) = 0;
    output_g.vu_ext_gen_high(rr)  = 0;
    output_g.vu_ext_gen_low(rr) = 0;
    output_g.nt_ext_gen_high(rr)  = 0;
    output_g.nt_ext_gen_low(rr) = 0;
    
    % unique marine animal generua
    genus_unq = unique(genus(mar_animal));
    
    % loop across genera
    for ii = 1:length(genus_unq)
        
        % species in this genus
        this_genus = strcmp(genus,genus_unq(ii));
        
        % marine species in this genus
        this_genus_mar = this_genus.*mar_animal;
        
        % convert to logical
        this_genus_mar = this_genus_mar>0;
        
        % total # species per category per genus
        output_g.ne(rr,ii) = nansum(strcmp(status(this_genus_mar),'Not Evaluated'));
        output_g.dd(rr,ii) = nansum(strcmp(status(this_genus_mar),'Data Deficient'));
        output_g.lc(rr,ii) = nansum(strcmp(status(this_genus_mar),'Least Concern'));
        output_g.nt(rr,ii) = nansum(strcmp(status(this_genus_mar),'Near Threatened'));
        output_g.vu(rr,ii) = nansum(strcmp(status(this_genus_mar),'Vulnerable'));
        output_g.en(rr,ii) = nansum(strcmp(status(this_genus_mar),'Endangered'));
        output_g.cr(rr,ii) = nansum(strcmp(status(this_genus_mar),'Critically Endangered'));
        output_g.ew(rr,ii) = nansum(strcmp(status(this_genus_mar),'Extinct in the Wild'));
        output_g.ex(rr,ii) = nansum(strcmp(status(this_genus_mar),'Extinct'));
        output_g.rex(rr,ii) = nansum(strcmp(status(this_genus_mar),'Regionally Extinct'));
        
       
        % add to genus total count, if all species are NOT 'DD'
        if output_g.dd(rr,ii) ~= nansum(this_genus_mar) && output_g.ne(rr,ii) ~= nansum(this_genus_mar) 
            output_g.genera(rr) = output_g.genera(rr) + 1;
        end
        
        % Extinction using VU+ threshold
        output_g.vu_thresh_high(rr,ii) = output_g.vu(rr,ii)+output_g.en(rr,ii)+output_g.cr(rr,ii)+output_g.ew(rr,ii)+output_g.ex(rr,ii)+output_g.rex(rr,ii);
        output_g.vu_thresh_low(rr,ii) = output_g.lc(rr,ii)+output_g.nt(rr,ii);
        
        % extinction - higher estimate (any species with >=VU => genus is extinct)
        if output_g.vu_thresh_high(rr,ii)>0
            output_g.vu_ext_gen_high(rr) =  output_g.vu_ext_gen_high(rr) + 1;
        end
        
        % extinction - lower estimate (any species with <=NT => is not extinct)
        if  (output_g.vu_thresh_low(rr,ii)<1) && (output_g.dd(rr,ii) ~= sum(this_genus_mar)) && (output_g.ne(rr,ii) ~= sum(this_genus_mar))
            output_g.vu_ext_gen_low(rr) =  output_g.vu_ext_gen_low(rr) + 1;
        end
        
        % Extinction using NT+ threshold
        output_g.nt_thresh_high(rr,ii) = output_g.nt(rr,ii)+output_g.vu(rr,ii)+output_g.en(rr,ii)+output_g.cr(rr,ii)+output_g.ew(rr,ii)+output_g.ex(rr,ii)+output_g.rex(rr,ii);
        output_g.nt_thresh_low(rr,ii) = output_g.lc(rr,ii);
        
        % extinction - higher estimate (any speices with >=NT => genus is extinct)
        if output_g.nt_thresh_high(rr,ii)>0
            output_g.nt_ext_gen_high(rr) =  output_g.nt_ext_gen_high(rr) + 1;
        end
        
        % extinction - lower estimate (any speices with ==LC => is not extinct)
        if  (output_g.nt_thresh_low(rr,ii)<1) && (output_g.dd(rr,ii) ~= sum(this_genus_mar)) && (output_g.ne(rr,ii) ~= sum(this_genus_mar))
            output_g.nt_ext_gen_low(rr) =  output_g.nt_ext_gen_low(rr) + 1;
        end
    end  
    
    % Extinction using VU+ threshold
    % percent endangered (low-estimate)
    output_g.vu_endanger_low(rr) = 100*output_g.vu_ext_gen_low(rr)./output_g.genera(rr);
    % percent endangered (high-estimate)
    output_g.vu_endanger_high(rr) = 100*output_g.vu_ext_gen_high(rr)./output_g.genera(rr);
    
    
    % Extinction using NT+ threshold
    % percent endangered (low-estimate)
    output_g.nt_endanger_low(rr) = 100*output_g.nt_ext_gen_low(rr)./output_g.genera(rr);
    % percent endangered (high-estimate)
    output_g.nt_endanger_high(rr) = 100*output_g.nt_ext_gen_high(rr)./output_g.genera(rr);
end


end 


% unique species across regions 
output.mar_spec_unq = unique(output.mar_spec(~ismissing(output.mar_spec)));
output.mar_class_unq = unique(output.mar_class(~ismissing(output.mar_class)));


%% Average across regions

% minimum # of species per regions
thresh = [0 10 100 500];

% average regions
for tt = 1:length(thresh)
    low(tt) = nanmean(output.endanger_low(output.total>thresh(tt)));
    low_sd(tt) = nanstd(output.endanger_low(output.total>thresh(tt)));
    high(tt) = nanmean(output.endanger_high(output.total>thresh(tt)));
    high_sd(tt) = nanstd(output.endanger_high(output.total>thresh(tt)));
end

plot(thresh,low,'r')
hold on
scatter(thresh,low,'r','filled')
plot(thresh,low+low_sd,'--r')
plot(thresh,low-low_sd,'--r')
plot(thresh,high,'b')
scatter(thresh,high,'b','filled')
plot(thresh,high+high_sd,'--b')
plot(thresh,high-high_sd,'--b')
ylim([0 100])







    
    
    
    
    
    
    
    
    
    
    
    
    
    












