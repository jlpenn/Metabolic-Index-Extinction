%% Summary statistics of Metabolic Index traits 
% (data: https://github.com/cadeutsch/Metabolic-Index-Traits)

% Ao (atm^-^1)
Amn = 2.99398;% mean of log(Ao)
Asd = 0.479501;% standard dev. of log(Ao)

% Eo (eV)
Emn = 0.398792;% mean of Eo
Esd = 0.269235;% standard dev. of Eo

% Phi-crit
Pmn = 1.12854;% mean of log(phi-crit)
Psd = 0.423362;% standard dev. of log(phi-crit)

%% Ao PDF
binranges = 5:5:55; % 5 to 55 1/atm (range of significant Ao values) 
pd = makedist('LogNormal','mu',Amn,'sigma',Asd);% make Ao distribution
x = binranges;% sample same as data
y = pdf(pd,x);% compute PDF
mod = y./max(y);% normaize so max is 1

% save weights
par.Aow = mod;
par.Ao = binranges;

%%  Eo PDF
binranges = -0.1:0.1:1.1;% -0.1 to 1.1 eV (range of significant Eo values)
pd = makedist('Normal','mu',Emn,'sigma',Esd);
x = binranges;
y = pdf(pd,x);
mod = y./max(y);

% save weights
par.Eow = mod;
par.Eo = binranges;

%% Phi-crit PDF 
binranges =1.5:0.5:7; % 1.5-7 active:resting metabolic rate (range of Phi-crit)
pd = makedist('LogNormal','mu',Pmn,'sigma',Psd);% fit distribution
x = binranges;% sample same as data
y = pdf(pd,x);% compute PDF
mod = (y./max(y));% normalize so max is 1

% save weights
par.Pcw = mod;
par.phicrit = binranges;

