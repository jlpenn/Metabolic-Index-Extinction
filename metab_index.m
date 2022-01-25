function phi = metab_index(Par,temp,o2,z3d,spind)
% function to compute the metabolic index 
if nargin<5
    spind=1;
end

Inan1=isnan(temp);
Inan2=isnan(o2);
pres=z3d;% depth
o2=max(o2,1); % set minimum o2 to 1 micromol/liter
salt=temp*0+35;  % set to constant salinity, which has negligible effect
Kh = O2solubility_edit(temp,salt,pres);  % O2 solubility

switch Par.gscheme
    case 'ko2'
        o2=o2.*kappa;
    case 'po2'
        o2=o2./Kh;
    case 'kpo2'
        o2=o2./Kh.*kappa;
end

Ao=Par.Ao(spind);
Eo=Par.Eo(spind);
kb=8.6173324e-5;
tk=273.15;

if Par.doTref ==0;
    expt=exp(-Eo./(kb*(temp+tk))); % No reference Temperature
else
    Tdif = (1./(temp+tk))-(1./(Par.Tref+tk));% Reference Temperature
    expt=exp((-Eo./kb).*Tdif); 
end

if Par.doB == 1;
    B =Par.B(spind);
    mexpo=Par.mfit(spind);
    phi=Ao*o2./(B^mexpo*expt); % biomass included
else
    phi=Ao*o2./(expt); % biomass normalized
end
phi(Inan1)=nan;
phi(Inan2)=nan;
