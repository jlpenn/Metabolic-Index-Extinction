function [Kh] = O2solubility_edit(T,S,z)

%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes solubility of O2 in seawater  
% The solubility (Kh) is based on Henry's Law, but includes a
% correction for the effect of hydrostatic pressure of the water column.
% This allows the partial pressure of O2 to be computed at depth from
%   pO2=O2./Kh;
%
% The relevant equation is (see Enns et al., J. Phys. Chem. 1964)
%  d(ln p)/dP = V/RT
% where p = partial pressure of O2, P = hydrostatic pressure
% V = partial molar volume of O2, R = gas constant, T = temperature
%
%%%%%%%%%%%%%%%%%%%%%%

%% constants

XiO2=0.209; % Mean atmospheric O2 mixing ratio
Patm=1; % Atm pressure

V=32e-6; % partial molar volume of O2 (m3/mol)
R=8.31; % Gas constant [J/mol/K]

db2Pa=1e4; % convert pressure: decibar to Pascal


%% Pressure correction

P=sw_pres(z,z*0); % seawater pressure [db] !! Warning - z*0 neglects gravity differences w/ latitude
dP=P*db2Pa;
pCor=exp(-V*dP./(R*(T+273.15)));

rho=sw_dens(S,T,P); % seawater density [kg/m3]

%% Solubility with pressure correction
Kh0=O2sol(S,T)/(Patm*XiO2); % solubility at surface (atm) pressure [mmol/m3/atm]
Kh=Kh0.*pCor.*rho.*1e-3;
