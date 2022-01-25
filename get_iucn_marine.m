% Code to read the IUCN red list endangered species and compute the fraction threatened to extinction
% env: including marine, kingdom: animalia 
clear;close all;

% Exclusion switches
do.ectotherm = 1; % restrict to water-breathers?
do.marter = 1; % include marine/terr. species?
do.match = 1; % match precision (1 = only exact matches)
do.pop = 0; % exclude populations (same species, multiple assessments)
do.phyla = 1; % do by phylum?
do.genus_ext = 1; % compute extinction of genera?
do.species_list = 1; % extract list of endandgered species names
do.three_phyla = 0; % restrict to well-sampled Met. Index phyla (Arthropoda, Chordata, Mollusca) 
do.complete = 0; % restrict to comprehensively assessed classes (0 = no; 1= yes; 2 = restrict to non-comprehensive)

%% load the Red list
load('path_to_iucn_data.mat')

% parse variables
tax = table2array(iucnspecieslistwormsmatched(:,1));
kingdom = string(table2array(iucnspecieslistwormsmatched(:,2)));
phylum = string(table2array(iucnspecieslistwormsmatched(:,3)));
class = string(table2array(iucnspecieslistwormsmatched(:,4)));
order = string(table2array(iucnspecieslistwormsmatched(:,5)));
family = string(table2array(iucnspecieslistwormsmatched(:,6)));
genus = string(table2array(iucnspecieslistwormsmatched(:,7)));
species = string(table2array(iucnspecieslistwormsmatched(:,8)));
infra_rank = table2array(iucnspecieslistwormsmatched(:,9));
infra_name = table2array(iucnspecieslistwormsmatched(:,10));
population = table2array(iucnspecieslistwormsmatched(:,11));
status = string(table2array(iucnspecieslistwormsmatched(:,12)));
env = string(table2array(iucnspecieslistwormsmatched(:,13)));
accepted = string(table2array(iucnspecieslistwormsmatched(:,18)));
match = string(table2array(iucnspecieslistwormsmatched(:,19)));
match_type = string(table2array(iucnspecieslistwormsmatched(:,20)));

% remove table
clear iucnspecieslistwormsmatched

%% clean data
% convert statuses from prior to 2001 to latest version
status(strcmp(status,'LR/nt')) = 'NT';% lower risk/near-threatened to NT
status(strcmp(status,'LR/cd')) = 'NT';% lower risk/conservation-depenendent to NT
status(strcmp(status,'LR/lc')) = 'LC';% lower risk/least-concern to LC


% Environment selection
if do.marter == 0
    % marine environments
    mar1 = strcmp(env,'Marine');
    mar2 = strcmp(env,'Marine/Brackish');
    mar3 = strcmp(env,'Marine/Brackish/Freshwater');
    mar4 = strcmp(env,'Marine/Freshwater');
    
    % combine marine environments
    mar = mar1+mar2+mar3+mar4;
else
    mar1 = strcmp(env,'Marine');
    mar2 = strcmp(env,'Marine/Brackish');
    mar3 = strcmp(env,'Marine/Brackish/Freshwater');
    mar4 = strcmp(env,'Marine/Freshwater');
    mar5 = strcmp(env,'Marine/Terrestrial');
    mar6 = strcmp(env,'Marine/Terrestrial/Brackish');
    mar7 = strcmp(env,'Marine/Terrestrial/Brackish/Freshwater');
    mar8 = strcmp(env,'Marine/Terrestrial/Freshwater');
    % combine marine environments
    mar = mar1+mar2+mar3+mar4+mar5+mar6+mar7+mar8;
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

% precise matching?
if do.match 
    match_exact = strcmp(match_type,'matchtype: exact');
    mar_animal = mar_animal.*match_exact;
    mar_animal = mar_animal>0;
end

% exclude non-global assessments?
if do.pop 
    pop_unq = strcmp(population,'');
    mar_animal = mar_animal.*pop_unq;
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

% restrict to comprehensively assessed taxa?
if do.complete == 1 
    % restrict to classes with complete sampling of groups 
    complete_taxa = strcmp(class,'ACTINOPTERYGII') | strcmp(class,'ANTHOZOA') | strcmp(class,'CHONDRICHTHYES') | strcmp(class,'GASTROPODA') | strcmp(class,'HYDROZOA') | strcmp(class,'MALACOSTRACA') | strcmp(class,'MYXINI') ;
   
    % restrict to marine animals
    mar_animal = mar_animal.*complete_taxa;
    mar_animal = mar_animal>0;
elseif do.complete == 2
    % restrict to classes with complete sampling of species 
    complete_taxa = ~strcmp(class,'ACTINOPTERYGII') & ~strcmp(class,'ANTHOZOA') & ~strcmp(class,'CHONDRICHTHYES') & ~strcmp(class,'GASTROPODA') & ~strcmp(class,'HYDROZOA') & ~strcmp(class,'MALACOSTRACA') & ~strcmp(class,'MYXINI') ;
   
    % restrict to marine animals
    mar_animal = mar_animal.*complete_taxa;
    mar_animal = mar_animal>0;
elseif do.complete == 3 
   % restrict to completely sampled groups
%     complete_taxa = strcmp(order,'CARCHARHINIFORMES') | strcmp(order,'HETERODONTIFORMES') | strcmp(order,'HEXANCHIFORMES') | strcmp(order,'LAMNIFORMES') | strcmp(order,'ORECTOLOBIFORMES') | strcmp(order,'PRISTIOPHORIFORMES')... 
%         | strcmp(order,'SQUALIFORMES') | strcmp(order,'SQUATINIFORMES') | strcmp(order,'MYLIOBATIFORMES') | strcmp(order,'RAJIFORMES') | strcmp(order,'RHINOPRISTIFORMES') | strcmp(order,'TORPEDINIFORMES') | strcmp(order,'TORPEDINIFORMES')...
%         | strcmp(family,'POMACANTHIDAE') | strcmp(family,'POMACANTHIDAE')  | strcmp(family,'CHAETODONTIDAE') | strcmp(family,'CHAETODONTIDAE') | strcmp(genus,'Megalops') | strcmp(family,'ELOPIDAE') | strcmp(family,'SCARIDAE') | strcmp(family,'ACIPENSERIDAE')...
%         | strcmp(genus,'Alphestes') | strcmp(genus,'Anyperodon') | strcmp(genus,'Aethaloperca') | strcmp(genus,'Cephalopholis') | strcmp(genus,'Cephalopholis') | strcmp(genus,'Chromileptes') | strcmp(genus,'Dermatolepis') | strcmp(genus,'Epinephelus') ...
%         | strcmp(genus,'Gonioplectrus')  | strcmp(genus,'Gracila') | strcmp(genus,'Hyporthodus')  | strcmp(genus,'Mycteroperca') | strcmp(genus,'Paranthias') | strcmp(genus,'Plectropomus') | strcmp(genus,'Saloptia') | strcmp(genus,'Triso') | strcmp(genus,'Variola')...
%         | strcmp(family,'LABRIDAE') | strcmp(genus,'Allothunnus') | strcmp(genus,'Auxis') | strcmp(genus,'Euthynnus') | strcmp(genus,'Katsuwonus') | strcmp(genus,'Thunnus');
%    
    filename = 'path_to_comprehensive_taxa.xlsx';
    % load taxa names 
    [n, t, ad] = xlsread(filename);
    
    % intialize vector
    complete_taxa = zeros(length(phylum),1);
    
    % taxanomic levels
    for cc = 2:5
        tax_unq = unique(t(:,cc));
        for rr = 2:length(tax_unq)-1
            tax = tax_unq(rr);
            if cc == 2
                tax_comp = strcmp(class,tax);
                complete_taxa = complete_taxa + tax_comp;
            elseif cc == 3
                tax_comp = strcmp(order,tax);
                complete_taxa = complete_taxa + tax_comp;
            elseif cc == 4
                tax_comp = strcmp(family,tax);
                complete_taxa = complete_taxa + tax_comp;
            elseif cc == 5
                tax_comp = strcmp(genus,tax);
                complete_taxa = complete_taxa + tax_comp;
            end
        end
    end 
    
    % restrict to marine animals
   
    complete_taxa = complete_taxa > 0;
    mar_animal = mar_animal.*complete_taxa;
    mar_animal = mar_animal>0;

end


%% Percent species endangerment

% total marine animal species per category
output.ne = nansum(strcmp(status(mar_animal),'NE'));
output.dd = nansum(strcmp(status(mar_animal),'DD'));
output.lc = nansum(strcmp(status(mar_animal),'LC'));
output.nt = nansum(strcmp(status(mar_animal),'NT'));
output.vu = nansum(strcmp(status(mar_animal),'VU'));
output.en = nansum(strcmp(status(mar_animal),'EN'));
output.cr = nansum(strcmp(status(mar_animal),'CR'));
output.ew = nansum(strcmp(status(mar_animal),'EW'));
output.ex = nansum(strcmp(status(mar_animal),'EX'));

% total species across categories
output.total_dd = output.dd+output.lc+output.nt+output.vu+output.en+output.cr+output.ew+output.ex;

% total species - DD 
output.total = output.lc+output.nt+output.vu+output.en+output.cr+output.ew+output.ex;

% threatened species (VU+)
output.threat_vu = output.vu+output.en+output.cr+output.ew+output.ex;

% threatened species (NT+) 
output.threat_nt = output.nt+output.vu+output.en+output.cr+output.ew+output.ex;

% percent endangered (low-estimate)
output.endanger_low = 100*output.threat_vu./output.total;

% percent endangered (high-estimate)
output.endanger_high = 100*output.threat_nt./output.total;

%% endangered species names (NT+)
if do.species_list 
    
    % marine species,statuses
    mar_species = species(mar_animal);
    mar_status = status(mar_animal);
   
    % endangered species counter 
    endngr_count = 1;
 
    % Endangered (NT+)
    for idx = 1:length(mar_status)
        if strcmp(mar_status(idx),'NT')
            output.endanger_list_nt(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'VU')
            output.endanger_list_nt(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'EN')
            output.endanger_list_nt(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'CR')
            output.endanger_list_nt(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'EW')
            output.endanger_list_nt(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'EX')
            output.endanger_list_nt(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        end
    end
       
    % endangered species counter 
    endngr_count = 1;
 
    % Endangered (VU+)
    for idx = 1:length(mar_status)
        if strcmp(mar_status(idx),'VU')
            output.endanger_list_vu(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'EN')
            output.endanger_list_vu(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'CR')
            output.endanger_list_vu(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'EW')
            output.endanger_list_vu(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        elseif strcmp(mar_status(idx),'EX')
            output.endanger_list_vu(endngr_count) = mar_species(idx);
            endngr_count = endngr_count+1;
        end
    end
end 
%% Percent species endangerment per phylum
if do.phyla
    % unique marine animal phyla
    phyla_unq = unique(phylum(mar_animal));
    
    % loop across phyla 
    for ii = 1:length(phyla_unq)
        
        % species in this phylum
        this_phyla = strcmp(phylum,phyla_unq(ii));
        
        % restrict to marine species
        this_phyla_mar = this_phyla.*mar_animal;
        
        % convert to logical
        this_phyla_mar = this_phyla_mar>0;
        
        % total species per category per phylum
        output_p.ne(ii) = nansum(strcmp(status(this_phyla_mar),'NE'));
        output_p.dd(ii) = nansum(strcmp(status(this_phyla_mar),'DD'));
        output_p.lc(ii) = nansum(strcmp(status(this_phyla_mar),'LC'));
        output_p.nt(ii) = nansum(strcmp(status(this_phyla_mar),'NT'));
        output_p.vu(ii) = nansum(strcmp(status(this_phyla_mar),'VU'));
        output_p.en(ii) = nansum(strcmp(status(this_phyla_mar),'EN'));
        output_p.cr(ii) = nansum(strcmp(status(this_phyla_mar),'CR'));
        output_p.ew(ii) = nansum(strcmp(status(this_phyla_mar),'EW'));
        output_p.ex(ii) = nansum(strcmp(status(this_phyla_mar),'EX'));
        
        % total species across categories
        output_p.total_dd(ii) = output_p.dd(ii)+output_p.lc(ii)+output_p.nt(ii)+output_p.vu(ii)+output_p.en(ii)+output_p.cr(ii)+output_p.ew(ii)+output_p.ex(ii);
        
        % total species - DD
        output_p.total(ii) = output_p.lc(ii)+output_p.nt(ii)+output_p.vu(ii)+output_p.en(ii)+output_p.cr(ii)+output_p.ew(ii)+output_p.ex(ii);
        
        % threatened species (VU+)
        output_p.threat_vu(ii) = output_p.vu(ii)+output_p.en(ii)+output_p.cr(ii)+output_p.ew(ii)+output_p.ex(ii);
        
        % threatened species (NT+)
        output_p.threat_nt(ii) = output_p.nt(ii)+output_p.vu(ii)+output_p.en(ii)+output_p.cr(ii)+output_p.ew(ii)+output_p.ex(ii);
    
        % percent endangered (low-estimate)
        output_p.endanger_low(ii) = 100*output_p.threat_vu(ii)./output_p.total(ii);

        % percent endangered (high-estimate)
        output_p.endanger_high(ii) = 100*output_p.threat_nt(ii)./output_p.total(ii);
        
        % total genera
        output_p.genera(ii) = length(unique(genus(this_phyla_mar)));

    end
end

%% Percent genus endangerment 
if do.genus_ext
    
    % initialize genus counts
    output_g.genera = 0;
    output_g.vu_ext_gen_high  = 0;
    output_g.vu_ext_gen_low = 0;
    output_g.nt_ext_gen_high  = 0;
    output_g.nt_ext_gen_low = 0;
    
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
        output_g.ne(ii) = nansum(strcmp(status(this_genus_mar),'NE'));
        output_g.dd(ii) = nansum(strcmp(status(this_genus_mar),'DD'));
        output_g.lc(ii) = nansum(strcmp(status(this_genus_mar),'LC'));
        output_g.nt(ii) = nansum(strcmp(status(this_genus_mar),'NT'));
        output_g.vu(ii) = nansum(strcmp(status(this_genus_mar),'VU'));
        output_g.en(ii) = nansum(strcmp(status(this_genus_mar),'EN'));
        output_g.cr(ii) = nansum(strcmp(status(this_genus_mar),'CR'));
        output_g.ew(ii) = nansum(strcmp(status(this_genus_mar),'EW'));
        output_g.ex(ii) = nansum(strcmp(status(this_genus_mar),'EX'));
        
       
        % add to genus total count, if all species are NOT 'DD'
        if output_g.dd(ii) ~= nansum(this_genus_mar) && output_g.ne(ii) ~= nansum(this_genus_mar)
            output_g.genera = output_g.genera + 1;
        end
        
        % Extinction using VU+ threshold
        output_g.vu_thresh_high(ii) = output_g.vu(ii)+output_g.en(ii)+output_g.cr(ii)+output_g.ew(ii)+output_g.ex(ii);
        output_g.vu_thresh_low(ii) = output_g.lc(ii)+output_g.nt(ii);
       
        % extinction - higher estimate (any speices with >=VU => genus is extinct)
        if output_g.vu_thresh_high(ii)>0
            output_g.vu_ext_gen_high =  output_g.vu_ext_gen_high + 1;
        end
        
        % extinction - lower estimate (any speices with <=NT => is not extinct)
        if  (output_g.vu_thresh_low(ii)<1) && (output_g.dd(ii) ~= sum(this_genus_mar)) && (output_g.ne(ii) ~= nansum(this_genus_mar))
           
            output_g.vu_ext_gen_low =  output_g.vu_ext_gen_low + 1;
        end
        
        % Extinction using NT+ threshold
        output_g.nt_thresh_high(ii) = output_g.nt(ii)+output_g.vu(ii)+output_g.en(ii)+output_g.cr(ii)+output_g.ew(ii)+output_g.ex(ii);
        output_g.nt_thresh_low(ii) = output_g.lc(ii);
        
        % extinction - higher estimate (any speices with >=NT => genus is extinct)
        if output_g.nt_thresh_high(ii)>0
            output_g.nt_ext_gen_high =  output_g.nt_ext_gen_high + 1;
        end
        
        % extinction - lower estimate (any speices with ==LC => is not extinct)
        if  (output_g.nt_thresh_low(ii)<1) && (output_g.dd(ii) ~= sum(this_genus_mar)) && (output_g.ne(ii) ~= nansum(this_genus_mar))
            output_g.nt_ext_gen_low =  output_g.nt_ext_gen_low + 1;
        end
    end  
    
    % Extinction using VU+ threshold
    % percent endangered (low-estimate)
    output_g.vu_endanger_low = 100*output_g.vu_ext_gen_low./output_g.genera;
    % percent endangered (high-estimate)
    output_g.vu_endanger_high = 100*output_g.vu_ext_gen_high./output_g.genera;
    
    
    % Extinction using NT+ threshold
    % percent endangered (low-estimate)
    output_g.nt_endanger_low = 100*output_g.nt_ext_gen_low./output_g.genera;
    % percent endangered (high-estimate)
    output_g.nt_endanger_high = 100*output_g.nt_ext_gen_high./output_g.genera;
end

