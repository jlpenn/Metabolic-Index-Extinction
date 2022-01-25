%% Code to tally threat factors for IUCN threatened species 
clear all;close all;

% get endangered species list
get_iucn_marine

% threat matrix (size: species, threats)
output.threat_fac_nt = zeros(length(output.endanger_list_nt),12);
output.threat_fac_vu = zeros(length(output.endanger_list_vu),12);

% skip missing files 
if do.three_phyla == 1
    skip = [423 426 429];
elseif do.complete == 1
    skip = [805 807 810];
elseif do.complete == 0
    skip = [829 832 835];
elseif do.complete == 3;
    skip = 674;
end

% Count threats by species (NT+)
for idx =1:length(output.endanger_list_nt)
    if sum(idx == skip)
    elseif sum(idx == skip)
    elseif sum(idx == skip)
    else
    % file name
    split_name = split(output.endanger_list_nt(idx));
    species_name = sprintf('iucn_marine_species_threats_%s_%s.csv',split_name(1),split_name(2));
    
    % file path
    species_path = sprintf('path_to_species_data/%s',species_name);
    
    % code to read csv files
    [code, title, timing, scope, severity, score, invasive] =textread(species_path,'%s %s %s %s %s %s %s','delimiter',',','headerlines',1,'emptyvalue',NaN);
    
    % Threat factors
    %residential & commerical development
    output.threat_fac_nt(idx,1)=sum(contains(code,'1.'))>0;
    
    % agriculture & aquaculture
    output.threat_fac_nt(idx,2)=sum(contains(code,'2.'))>0;
    
    % energy production & mining
    output.threat_fac_nt(idx,3)=sum(contains(code,'3.'))>0;
    
    % transportation & service corridoors
    output.threat_fac_nt(idx,4)=sum(contains(code,'4.'))>0;
    
    % biological resource use
    output.threat_fac_nt(idx,5)=sum(contains(code,'5.'))>0;
    
    % human intrusions & disturbance
    output.threat_fac_nt(idx,6)=sum(contains(code,'6.'))>0;
    
    % natural system modification
    output.threat_fac_nt(idx,7)=sum(contains(code,'7.'))>0;
    
    % invasive species, disease, genes
    output.threat_fac_nt(idx,8)=sum(contains(code,'8.'))>0;
    
    % pollution
    output.threat_fac_nt(idx,9)=sum(contains(code,'9.'))>0;
    
    % geological events
    output.threat_fac_nt(idx,10)=sum(contains(code,'10.'))>0;
    
    % climate change
    output.threat_fac_nt(idx,11)=sum(contains(code,'11.'))>0;
    
    % other
    output.threat_fac_nt(idx,12)=sum(contains(code,'12.'))>0;
    end 
    
end    

% summary statistics (NT+)
output.total_threat_nt = sum(sum(output.threat_fac_nt,2)>0);
output.threat_pct_nt = 100*sum(output.threat_fac_nt)./output.total_threat_nt; 

% skip missing files 
if do.three_phyla == 1
    skip = [302 304 307];
elseif do.complete == 1
    skip = [515 517 520];
elseif do.complete == 0
    skip = [534 536 539];
elseif do.complete == 3;
    skip = 404;
end

% Count threats by species (VU+)
for idx =1:length(output.endanger_list_vu)
    if sum(idx == skip)
    elseif sum(idx == skip)
    elseif sum(idx == skip)
    else
    % file name
    split_name = split(output.endanger_list_vu(idx));
    species_name = sprintf('iucn_marine_species_threats_%s_%s.csv',split_name(1),split_name(2));
    
    % file path
    species_path = sprintf('path_to_species_data/%s',species_name);
    
    % code to read csv files
    [code, title, timing, scope, severity, score, invasive] =textread(species_path,'%s %s %s %s %s %s %s','delimiter',',','headerlines',1,'emptyvalue',NaN);
    
    % Threat factors
    %residential & commerical development
    output.threat_fac_vu(idx,1)=sum(contains(code,'1.'))>0;
    
    % agriculture & aquaculture
    output.threat_fac_vu(idx,2)=sum(contains(code,'2.'))>0;
    
    % energy production & mining
    output.threat_fac_vu(idx,3)=sum(contains(code,'3.'))>0;
    
    % transportation & service corridoors
    output.threat_fac_vu(idx,4)=sum(contains(code,'4.'))>0;
    
    % biological resource use
    output.threat_fac_vu(idx,5)=sum(contains(code,'5.'))>0;
    
    % human intrusions & disturbance
    output.threat_fac_vu(idx,6)=sum(contains(code,'6.'))>0;
    
    % natural system modification
    output.threat_fac_vu(idx,7)=sum(contains(code,'7.'))>0;
    
    % invasive species, disease, genes
    output.threat_fac_vu(idx,8)=sum(contains(code,'8.'))>0;
   
    % pollution
    output.threat_fac_vu(idx,9)=sum(contains(code,'9.'))>0;
    
    % geological events
    output.threat_fac_vu(idx,10)=sum(contains(code,'10.'))>0;
    
    % climate change
    output.threat_fac_vu(idx,11)=sum(contains(code,'11.'))>0;
    
    % other
    output.threat_fac_vu(idx,12)=sum(contains(code,'12.'))>0;
    end 
    
end    

% summary statistics (VU+)
output.total_threat_vu = sum(sum(output.threat_fac_vu,2)>0);
output.threat_pct_vu = 100*sum(output.threat_fac_vu)./output.total_threat_vu; 

%% Climate threat (% of all species)

output.climate_vu = (output.threat_pct_vu(11)/100).*output.endanger_low;
output.climate_nt = (output.threat_pct_nt(11)/100).*output.endanger_high;
