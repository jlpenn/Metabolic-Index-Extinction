%% Code to extract model output and reproduce the main figures of 'Avoiding Ocean Mass Extinction from Climate Warming"
clear; close all;
get_6me_output;

%% Figure 1 - panel A
% Fossil record from Rohde and Muller, Cycles in Fossil Diversity,Nature (2005)
load('/rraid1/jpenn/Data/6me/Sepkoski_fossil_data.mat')
ext = obs.ext; % Extinction proportion (%) for each 1 myr interval from 542 Ma to present

% PDF of extinction proportion
binranges = 0:2.5:100; % bin edges for PDF
binmid = 1.25:2.5:98.75;% bin midpoints for PDF 
[bincounts,e] = histcounts(ext,binranges);% bin counts

% (1-CDF) of extinction proportion
ext_levels = 0:1:100; % extinction levels for (1-CDF) 
for zz = 1:length(ext_levels)
    thisc= ext_levels(zz);% extinction threshold
    ext_cdf(zz) =nansum(ext>=thisc);% sum time-periods with extinctions above threshold
end

% normalize counts to total
ext_cdf=ext_cdf./nanmax(ext_cdf(:));

% Make a matrix of (1-CDF) vs. temperature (independent relationship)
ext_cdfM =zeros(length(ext_levels),10);% size = ext_levels x 10 (arbitrary for plotting)
for mod = 1:10
     ext_cdfM(:,mod)=ext_cdf;% fill matrix
end

% Plot Fig. 1A
figure(1)
set(gcf,'units','inches','position',[2,2,6,4])
subplot(1,5,1:2);

% Shade CDF
% x/y axis variables for plotting
xaxis = 0:nanmax(bincounts)/9:nanmax(bincounts);% temperature axis (length = size(ext_cdfM,2));
yaxis = ext_levels;% extinction proportion
v = [.05 .25 .5 .95];% 1-percentiles for plotting
hold on
[C,h,Cf]=contourf(xaxis,yaxis,ext_cdfM,v);% contour CDF
for q=1:length(h)
    set(h(q),'LineStyle','none')
end
shading flat
hold on
colormap(flipud(spring(5)))

text(117.1,2.5,'5th','color',[0 0 0],'fontsize',11)
text(117.1,8.5,'50th','color',[0 0 0],'fontsize',11)
text(117.1,14.5,'75th','color',[0 0 0],'fontsize',11)
text(117.1,28.5,'95th','color',[0 0 0],'fontsize',11)

% Plot mass extinctions
% Late Ordivician
plot([0,117],[45.5 45.5],'color',[.7 .7 .7],'linewidth',1)

% Late Devonian
plot([0,117],[18 18],'color',[.7 .7 .7],'linewidth',1)

% Late Permian
plot([0,117],[62 62],'color',[.7 .7 .7],'linewidth',1)

% Triassic
plot([0,117],[28.85 28.85],'color',[.7 .7 .7],'linewidth',1)

% Cretaceous
plot([0,117],[39 39],'color',[.7 .7 .7],'linewidth',1)

% Plot extinction PDF
barh(binmid,bincounts,1,'facecolor',[.7 .7 .7],'facealpha',.5)
box on
set(gca,'fontsize',14)
xlabel('Past Occurrences')
ylabel('Extinction (%)')

clearvars -except output high low hist obs 

%% Figure 1 - panel B
figure(2)
xsplit = 8;

% Permian extinction 
ptr.tas = 18.8;% atmospheric delta T (deg. C)
ptr.em = 40; % Extinction (%) (100% dispersal)
ptr.es = 60; % Extinction (%) (0% dispersal)

% Ensemble average 
% rcp 8.5
w8 = high.dT;
w8(1:9,5)=nan;% remove last data point from average(repeated)
w8(10:14,end)=nan;% remove last data point from average(repeated)
es8 = high.ext_loss; % Extinction (0% dispersal)
es8(1:9,5)=nan;% remove last data point from average(repeated)
es8(10:14,end)=nan;% remove last data point from average(repeated)
em8 = high.ext_lossgain; % Extinction (100% dispersal)
em8(1:9,5)=nan;% remove last data point from average(repeated)
em8(10:14,end)=nan;% remove last data point from average(repeated)

% rcp 2.6
w2 = low.dT;
w2(:,5)=nan;% remove last data point from average(repeated)
es2 = low.ext_loss;% Extinction (0% dispersal)
es2(:,5)=nan;% remove last data point from average(repeated)
em2 = low.ext_lossgain;% Extinction (100% dispersal)
em2(:,5)=nan;% remove last data point from average(repeated)


% rcp 8.5 warming (average models in 5 deg intervals)
t = 0:5:20;

for ii = 1:length(t)
    if ii ==1
        idx = w8<=t(ii);
        es8mn(ii) = nanmean(es8(idx));
        em8mn(ii) = nanmean(em8(idx));
        w8mn(ii) = nanmean(w8(idx));
    else
        idx = w8<=t(ii) & w8>t(ii-1);
        es8mn(ii) = nanmean(es8(idx));
        em8mn(ii) = nanmean(em8(idx));
        w8mn(ii) = nanmean(w8(idx));
    end
end

w8mn_save = w8mn; % save for Fig 3

% rcp 2.6 warming (average models in 1 deg intervals (finer resolution b/c narrower ∆T range))
t = 0:1:20;

for ii = 1:length(t)
    if ii ==1
        idx = w2<=t(ii);
        es2mn(ii) = nanmean(es2(idx));
        em2mn(ii) = nanmean(em2(idx));
        w2mn(ii) = nanmean(w2(idx));
    else
        idx = w2<=t(ii) & w2>t(ii-1);
        es2mn(ii) = nanmean(es2(idx));
        em2mn(ii) = nanmean(em2(idx));
        w2mn(ii) = nanmean(w2(idx));
    end
end

% 2020 average 
% warming (average across scenarios)
w2_2020 = low.dT(:,2);
w8_2020 = high.dT(:,2);
w2020 = [w2_2020;w8_2020];
w2020= nanmean(w2020);
w2020 = [0 w2020];


% no dispersal extinction (average across scenarios)
es2_2020 = low.ext_loss(:,2);
es8_2020 = high.ext_loss(:,2);
es2020 = [es2_2020;es8_2020];
es2020= nanmean(es2020);
es2020 = [0 es2020];

% dispersal extinction (average across scenarios)
em2_2020 = low.ext_lossgain(:,2);
em8_2020 = high.ext_lossgain(:,2);
em2020 = [em2_2020;em8_2020];
em2020= nanmean(em2020);
em2020 = [0 em2020];

% Plot figure 1B
set(gcf,'units','inches','position',[2,2,6,4])
haxis(1)=subplot(1,5,1:3);

% Mass Extinctions 
% Late Ordivician
plot([0,20],[45.5 45.5],'color',[.7 .7 .7],'linewidth',1)
hold on
% Late Devonian
plot([0,20],[18 18],'color',[.7 .7 .7],'linewidth',1)

% Late Permian
plot([0,20],[62 62],'color',[.7 .7 .7],'linewidth',1)

% Triassic
plot([0,20],[28.85 28.85],'color',[.7 .7 .7],'linewidth',1)

% Cretaceous
plot([0,20],[39 39],'color',[.7 .7 .7],'linewidth',1)


% IUCN endangered species
patch_x = [0 20 20 0];
patch_y = [9.4 9.4 13.8 13.8]; % threatened + endangered species 
p=patch(patch_x(:),patch_y(:),'k');set(p,'facecolor',[1 .7 .4],'edgecolor',[0 0 0],'facealpha',1,'edgealpha',0);hold on


scatter(ptr.tas(end), ptr.em(end),70,'k','filled')
scatter(ptr.tas(end), ptr.es(end),70,'k','filled')

% High emissions future
for mod = 1:14
     plot((high.dT(mod,:)),high.ext_loss(mod,:),'color',[.8 0.4 0.4],'linewidth',1);
     plot((high.dT(mod,:)),high.ext_lossgain(mod,:),'color',[.8 0.4 0.4],'linewidth',1);
end

% Highlight historical emissions
for mod = 1:14
     plot((high.dT(mod,1:2)),high.ext_loss(mod,1:2),'color',[.4 0.4 0.4],'linewidth',1);
     plot((high.dT(mod,1:2)),high.ext_lossgain(mod,1:2),'color',[.4 0.4 0.4],'linewidth',1);
end

% ensemble average
plot(w8mn,es8mn,'linewidth',3,'color',[.7 0 0])
plot(w8mn,em8mn,'linewidth',3,'color',[.7 0 0])

% Low emissions future
for mod = 1:10
     plot((low.dT(mod,:)),low.ext_loss(mod,:),'color',[.4 .4 .8],'linewidth',1);
     plot((low.dT(mod,:)),low.ext_lossgain(mod,:),'color',[.4 .4 .8],'linewidth',1);
end


for mod = 1:10
     plot((low.dT(mod,1:2)),low.ext_loss(mod,1:2),'color',[.4 .4 .4],'linewidth',1);
     plot((low.dT(mod,1:2)),low.ext_lossgain(mod,1:2),'color',[.4 .4 .4],'linewidth',1);
end
% ensemble average
plot(w2mn,es2mn,'linewidth',3,'color',[0 0 .7])
plot(w2mn,em2mn,'linewidth',3,'color',[.0 0 .7])

% 2020 average across scenarios
plot(w2020,es2020,'linewidth',3,'color',[0 0 0])
plot(w2020,em2020,'linewidth',3,'color',[0 0 0])

% RCP8.5 2100 warming
x = [nanmin(high.dT(:,5)) nanmax(high.dT(:,4))];
y = [100,100];
plot(x,y,'color',[.7 0 0],'linewidth',3)
text(x(1)-1.5,y(1)+2.5,'2100','color','k','fontsize',10)
scatter(nanmean(high.dT(:,4)),100,70,[.7 0 0],'filled')

%RCP2.6: 2100 warming
x = [nanmin(low.dT(:,end-1)) nanmax(low.dT(:,end-1))];
y = [100,100];
plot(x,y,'color',[0 0 .7],'linewidth',3)
scatter(nanmean(low.dT(:,end-1)),100,70,[0 0 .7],'filled')


% 2020 warming
x = [nanmin(low.dT(:,2)) nanmax(high.dT(:,2))];
y = [100,100];
plot(x,y,'k','linewidth',3)
text(x(1)-1,y(1)+2.5,'2020','color','k','fontsize',10)
scatter(1.28,100,70,'k','filled')

xlim([0 xsplit])
ylim([0 100])
box on
set(gcf,'color',[1 1 1])
set(gca,'fontsize',16,'xtick',[0 1 2 3 5 7])
ylabel('Extinction (%)')
xlabel('∆ Global Temperature ({\circ}C)')

% Beyond 2100 

haxis(2)=subplot(1,5,4:5);
pos = get(haxis(2), 'Position');
pos(1) = .5798;
set(haxis(2),'Position',pos)

% Mass Extinctions 
% Late Ordivician
plot([0,20],[45.5 45.5],'color',[.7 .7 .7],'linewidth',1)
hold on
% Late Devonian
plot([0,20],[18 18],'color',[.7 .7 .7],'linewidth',1)

% Late Permian
plot([0,20],[62 62],'color',[.7 .7 .7],'linewidth',1)

% Triassic
plot([0,20],[28.85 28.85],'color',[.7 .7 .7],'linewidth',1)

% Cretaceous
plot([0,20],[39 39],'color',[.7 .7 .7],'linewidth',1)


% IUCN endangered species
patch_x = [0 20 20 0];
patch_y = [9.4 9.4 13.8 13.8]; % threatened + endangered species 
p=patch(patch_x(:),patch_y(:),'k');set(p,'facecolor',[1 .7 .4],'edgecolor',[0 0 0],'facealpha',1,'edgealpha',0);hold on

% Permian extinction risks
scatter(ptr.tas, ptr.em,70,'k','filled')
scatter(ptr.tas, ptr.es,70,'k','filled')
xp = [ptr.tas ptr.tas];
yp = [ptr.em ptr.es];
plot(xp,yp,'k','linewidth',3)

% High emissions future
for mod = 1:14
     plot((high.dT(mod,:)),high.ext_loss(mod,:),'color',[.8 0.4 0.4],'linewidth',1);%[1 .95 .45],'linewidth',3)
     plot((high.dT(mod,:)),high.ext_lossgain(mod,:),'color',[.8 0.4 0.4],'linewidth',1);%[1 .95 .45],'linewidth',3)
end

% ensemble average
plot(w8mn,es8mn,'linewidth',3,'color',[.7 0 0])
plot(w8mn,em8mn,'linewidth',3,'color',[.7 0 0])

%2300 warming
x = [nanmin(high.dT(:,end)) nanmax(high.dT(:,end))];
y = [100,100];
plot(x,y,'color',[.7 0 0],'linewidth',3)
text(x(1)+1,y(1)+2.5,'2300','color','k','fontsize',10)
scatter(nanmean(high.dT(:,end)),100,70,[.7 0 0],'filled')


set(gca,'fontsize',16,'xtick',[10 15 20],'Yticklabels',[])
ylabel('')
xlabel('')
xlim([xsplit 20])
ylim([0 100])
box on
set(gcf,'color',[1 1 1])

clearvars -except output high low hist obs es8mn em8mn w8mn_save

%% Figure 1- panel C
figure(3)
xsplit = 8; % ∆T for xaxis split

% Permian extinction 
ptr.tas = 18.8;% atmospheric delta T (deg. C)
ptr.e = 62; % Extirpation (%)

% rcp 8.5
w8 = high.dT;
w8(1:9,5)=nan; % remove last data point from average(repeated)
w8(10:14,end)=nan; % remove last data point from average(repeated)
e8 = high.extr;
e8(1:9,5)=nan; % remove last data point from average(repeated)
e8(10:14,end)=nan; % remove last data point from average(repeated)

% rcp 2.6
w2 = low.dT;
w2(:,5)=nan; % remove last data point from average (repeated)
e2 = low.extr;
e2(:,5)=nan; % remove last data point from average (repeated)


% rcp 8.5 warming (average models in 5 deg intervals)
t = 0:5:20;

for ii = 1:length(t)
    if ii ==1
        idx = w8<=t(ii);
        e8mn(ii) = nanmean(e8(idx));
        w8mn(ii) = nanmean(w8(idx));
    else
        idx = w8<=t(ii) & w8>t(ii-1);
        e8mn(ii) = nanmean(e8(idx));
        w8mn(ii) = nanmean(w8(idx));
    end
end

% rcp 2.6 warming (average models in 1 deg intervals (finer resolution b/c narrower ∆T range))
t = 0:1:20;

for ii = 1:length(t)
    if ii ==1
        idx = w2<=t(ii);
        e2mn(ii) = nanmean(e2(idx));
        w2mn(ii) = nanmean(w2(idx));
    else
        idx = w2<=t(ii) & w2>t(ii-1);
        e2mn(ii) = nanmean(e2(idx));
        w2mn(ii) = nanmean(w2(idx));
    end
end

% 2020 average 

% warming (average across scenarios)
w2_2020 = low.dT(:,2);
w8_2020 = high.dT(:,2);
w2020 = [w2_2020;w8_2020];
w2020= nanmean(w2020);
w2020 = [0 w2020];


%  extirpation (average across scenarios)
e2_2020 = low.extr(:,2);
e8_2020 = high.extr(:,2);
e2020 = [e2_2020;e8_2020];
e2020= nanmean(e2020);
e2020 = [0 e2020];


% Plot Figure 1C
set(gcf,'units','inches','position',[2,2,6,4])
haxis(1)=subplot(1,5,1:3);

% dummy for legend
plot(ptr.e(:),-1000*ptr.e(:),'color',[0 0 0],'linewidth',3)
hold on
plot(ptr.e(:),-1000*ptr.e(:),'color',[0.1 0.1 1],'linewidth',3)
plot(ptr.e(:),-1000*ptr.e(:),'color',[.7 0 0],'linewidth',3)

% Endangered species 
patch_x = [0 20 20 0];
patch_y = [11 11 15.9 15.9];
p=patch(patch_x(:),patch_y(:),'k');set(p,'facecolor',[1 0.5 0],'edgecolor',[0 0 0],'facealpha',.5,'edgealpha',0);hold on

% High emissions future
for mod = 1:14
     plot(high.dT(mod,:),high.extr(mod,:),'color',[.8 0.4 0.4],'linewidth',1);
end
for mod = 1:9
     plot(high.dT(mod,:),high.extr(mod,:),'color',[.8 0.4 0.4],'linewidth',1);
end

% High emissions historical
for mod = 1:14
     plot(high.dT(mod,1:2),high.extr(mod,1:2),'color',[.4 0.4 0.4],'linewidth',1);
end
for mod = 1:9
     plot(high.dT(mod,1:2),high.extr(mod,1:2),'color',[.4 0.4 0.4],'linewidth',1);
end

% ensemble average
plot(w8mn,e8mn,'color',[.7 0 0],'linewidth',4)

% Low emissions future
for mod = 1:10
     plot(low.dT(mod,:),low.extr(mod,:),'color',[0.4 .4 .8],'linewidth',1);%[1 .95 .45],'linewidth',3)
end

% Low emissions historical
for mod = 1:10
     plot(low.dT(mod,1:2),low.extr(mod,1:2),'color',[0.4 0.4 0.4],'linewidth',1);%[1 .95 .45],'linewidth',3)
end

% 2.6 ensemble average
plot(w2mn,e2mn,'color',[0.1 0.1 1],'linewidth',4)

% 2020 average across scenarios
plot(w2020,e2020,'linewidth',4,'color',[0 0 0])

% RCP8.5 2100 warming
x = [nanmin(high.dT(:,5)) nanmax(high.dT(:,4))];
y = [100,100];
plot(x,y,'color',[.7 0 0],'linewidth',3)
text(x(1)-1.5,y(1)+2.5,'2100','color','k','fontsize',10)
scatter(nanmean(high.dT(:,4)),100,70,[.7 0 0],'filled')

% RCP2.6: 2100 warming
x = [nanmin(low.dT(:,end-1)) nanmax(low.dT(:,end-1))];
y = [100,100];
plot(x,y,'color' ,[0 0 .7],'linewidth',3)
scatter(nanmean(low.dT(:,end-1)),100,70,[0 0 .7],'filled')

% 2020 warming
x = [nanmin(low.dT(:,2)) nanmax(high.dT(:,2))];
y = [100,100];
plot(x,y,'k','linewidth',3)
text(x(1)-1,y(1)+2.5,'2020','color','k','fontsize',10)
scatter(1.28,100,70,'k','filled')

xlim([0 xsplit])
ylim([0 100])
set(gca,'fontsize',16,'xtick',[0 1 2 3 5 7])
ylabel('Extirpation (%)')
xlabel('∆ Global Temperature ({\circ}C)')
text(14,70,'Permian (model)','fontsize',12)
text(20.1,14,'Modern','color','k','fontsize',11)
box on
set(gcf,'color',[1 1 1])
legend('Historical','2.6 W/m^2','8.5 W/m^2')

% Beyond 2100
haxis(2)=subplot(1,5,4:5);
pos = get(haxis(2), 'Position');
pos(1) = .5798;
set(haxis(2),'Position',pos)

% Endangered species 
patch_x = [0 20 20 0];
patch_y = [11 11 15.9 15.9];
p=patch(patch_x(:),patch_y(:),'k');set(p,'facecolor',[1 0.5 0],'edgecolor',[0 0 0],'facealpha',.5,'edgealpha',0);hold on

% High emissions future
for mod = 1:14
     plot(high.dT(mod,:),high.extr(mod,:),'color',[.8 0.4 0.4],'linewidth',1);
end
for mod = 1:9
     plot(high.dT(mod,:),high.extr(mod,:),'color',[.8 0.4 0.4],'linewidth',1);
end

% ensemble average
plot(w8mn,e8mn,'color',[.7 0 0],'linewidth',4)

% "The Great Dying" + extirpation errorbar
scatter(ptr.tas(end),ptr.e(end),70,'k','filled')
xx = [ptr.tas(end) ptr.tas(end)];
yy = [ptr.e(end)-21.6 ptr.e(end)+21.6];
plot(xx,yy,'color','k','linewidth',2)

xx = [ptr.tas(end)-.2 ptr.tas(end)+.2];
yy = [ptr.e(end)-21.6 ptr.e(end)-21.6];
plot(xx,yy,'color','k','linewidth',2)

xx = [ptr.tas(end)-.2 ptr.tas(end)+.2];
yy = [ptr.e(end)+21.6 ptr.e(end)+21.6];
plot(xx,yy,'color','k','linewidth',2)

scatter(ptr.tas(end),ptr.e(end),70,'k')


% 2300 warming
x = [nanmin(high.dT(:,end)) nanmax(high.dT(:,end))];
y = [100,100];
plot(x,y,'color',[.7 0 0],'linewidth',3)
text(x(1)+1,y(1)+2.5,'2300','color','k','fontsize',10)
scatter(nanmean(high.dT(:,end)),100,70,[.7 0 0],'filled')

xlim([xsplit 20])
ylim([0 100])
set(gca,'fontsize',16)
box on
set(gcf,'color',[1 1 1])
set(gca,'fontsize',16,'xtick',[10 15 20],'Yticklabels',[])

clearvars -except output high low hist obs es8mn em8mn w8mn_save
%% Figure 2
figure(4)
% Panel 2A

% Ocean/land mask
coasts = obs.coasts;

% SeaAroundUs.org locations (lat,lon) of productive fisheries (1950-2014)
Fishlat = obs.Fishlat;
Fishlon = obs.Fishlon;

% Average climate models for 2100
Extr_2100 = output.high.extr_2100_map;

% LONG,LAT
% Model coordinates
x = 0.5:1:359.5;
y = -89.5:1:89.5;

% Coastal coordinates
par.x = 0.25:0.5:359.75;
par.y2 = 0.25:.5:180;
 
pcolor(x,y,squeeze(Extr_2100));hold on ;shading flat;
caxis([0 100])
colorbar
set(colorbar,'fontsize',14,'Ticks',[0 25 50 75 100],'location','south')
cmap = colormap(flipud(hot));
colormap(cmap)
hold on
v = [1,1];
contour(par.x,par.y2-90,coasts,v,'k','linewidth',1)
scatter(Fishlon(:),Fishlat(:),2.5,'b','filled')

set(gca,'color',[0.7 0.7 0.7],'fontsize',14)
set(gca,'fontsize',14,'Ytick',[-89 -60 -30 0 30 60 89],'Yticklabels',{'90^oS','60^oS','30^oS','0','30^oN','60^oN','90^oN'},'fontsize',14)
set(gca,'Xtick',[1 60 120 180 240 300 359])
set(gca,'Xticklabels',{'0','60^oE','120^oE','180^o','120^oW','60^oW','0'},'fontsize',14)
box on
xlabel('Longitude','fontsize',14)
ylabel('Latitude','fontsize',14)
set(gcf,'color',[1 1 1])
text(435,25,'Extirpation (%)','fontsize',14,'rotation',-90)
ylim([-78 89])

clearvars -except output high low hist obs es8mn em8mn w8mn_save

%% Panel 2B
% Make figure 2B
figure(5)
set(gcf,'units','inches','position',[2,2,4,5.5])

lat = -89.5:1:89.5;

% plot extirpation vs. latitude
plot(hist.extr_lat,lat,'color',[0 0 0],'linewidth',3)
hold on
plot(low.extr_lat_2100,lat,'color',[0.1 0.1 1],'linewidth',3)
plot(high.extr_lat_2100,lat,'color',[1 .7 0],'linewidth',3)
plot(high.extr_lat_2300,lat,'color',[1 0 0],'linewidth',3)

% High emissions
% 2300 
mx3 = high.extr_lat_2300+high.extr_lat_2300_sd;
min3 = high.extr_lat_2300-high.extr_lat_2300_sd;
patch_x = [lat(isnan(mx3)==0) fliplr(lat(isnan(mx3)==0))];
patch_y = [min3(isnan(mx3)==0) fliplr(mx3(isnan(mx3)==0))];
p=patch(patch_y,patch_x,'k');set(p,'facecolor',[.9 0.7 0.7],'edgecolor',[1 0 0],'facealpha',.5,'edgealpha',1);hold on
plot(high.extr_lat_2300,lat,'color',[1 0 0],'linewidth',3)

% 2100
mx2 = high.extr_lat_2100+high.extr_lat_2100_sd;
min2 = high.extr_lat_2100-high.extr_lat_2100_sd;
patch_x = [lat(isnan(mx2)==0) fliplr(lat(isnan(mx2)==0))];
patch_y = [min2(isnan(mx2)==0) fliplr(mx2(isnan(mx2)==0))];
p=patch(patch_y,patch_x,'k');set(p,'facecolor',[1 .7 0],'edgecolor',[1 .7 0],'facealpha',.5,'edgealpha',1);hold on
plot(high.extr_lat_2100,lat,'color',[1 .7 0],'linewidth',3)

% Low emissions
mx4 = low.extr_lat_2100+low.extr_lat_2100_sd;
min4 = low.extr_lat_2100-low.extr_lat_2100_sd;
patch_x = [lat(isnan(mx4)==0) fliplr(lat(isnan(mx4)==0))];
patch_y = [min4(isnan(mx4)==0) fliplr(mx4(isnan(mx4)==0))];
p=patch(patch_y,patch_x,'k');set(p,'facecolor',[.7 .7 .9],'edgecolor',[0 0 0.8],'facealpha',.5,'edgealpha',0.5);hold on
plot(low.extr_lat_2100,lat,'color',[0.1 0.1 1],'linewidth',3)

% Historical 
mx1 = hist.extr_lat+hist.extr_lat_sd;
min1 = hist.extr_lat-hist.extr_lat_sd;
patch_x = [lat(isnan(mx1)==0) fliplr(lat(isnan(mx1)==0))];
patch_y = [min1(isnan(mx1)==0) fliplr(mx1(isnan(mx1)==0))];
p=patch(patch_y,patch_x,'k');set(p,'facecolor',[.7 .7 .7],'edgecolor',[0 0 0],'facealpha',.5,'edgealpha',0.5);hold on
plot(hist.extr_lat,lat,'color',[0 0 0],'linewidth',3)

set(gca,'fontsize',14);
set(gcf,'color',[1 1 1]);
xlim([0 100])
ylabel('Latitude','fontsize',14)
set(gca,'Ytick',[-89 -60 -30 0 30 60 88])
set(gca,'Yticklabels',{'90^oS','60^oS','30^oS','0','30^oN','60^oN','90^oN'},'fontsize',14)

xlabel('Extirpation (%)','fontsize',14)
ylim([-75 88])
box on

legend('2020','2100 (2.6)','2100 (8.5)','2300 (8.5)')

clearvars -except output high low hist obs es8mn em8mn w8mn_save lat

%% Panel 2C
figure(6)
set(gcf,'units','inches','position',[2,2,4,5.5])

% Model historical richness pattern vs. latitude (normalized)
rich0 = hist.rich_lat_0m; % maximum summation depth of 0 meters
rich500 = hist.rich_lat_500m;% maximum summation depth of 500 meters
rich5000 = hist.rich_lat_5000m;% maximum summation depth of 5000 meters

% Observed richness from Chadhary, Saeedi, Costello, Marine Species Richness is bimodal with latitude: A reply to Fendandez and Marques,Trends in Ecology and Evolution (2017)
richo = obs.richo; % observed richness (number of species) vs. latitude (central estimate)
richh = obs.richh; % upper errorbar
richl = obs.richl; % lower errorbar

% Maximum richness in observations
norm = nanmax(richo(:));

% Normalize obs to max(obs)
richo = richo./norm;
richh = richh./norm;
richl = richl./norm;

% Plot 
% for legend
plot(rich500,output.lat5deg,'k','linewidth',3)
hold on
scatter(richo,output.lat5deg(2:end-1),100,[0 0 0],'filled')

% Model 
patch_x = [output.lat5deg fliplr(output.lat5deg)];
patch_y = [nanmin(rich0,rich5000) fliplr(nanmax(rich0,rich5000))];
p=patch(patch_y,patch_x,'k');set(p,'facecolor',[.7 .7 .7],'edgecolor',[0 0 0],'facealpha',.5,'edgealpha',0);hold on
hold on
plot(rich500,output.lat5deg,'k','linewidth',3)

% Observations
scatter(richo,output.lat5deg(2:end-1),100,[0 0 0],'filled')
hold on
h=errorbar_x(richo,output.lat5deg(2:end-1),richl,richh,'.k');
set(h,'linewidth',1.5)
legend('Model','Observations')
set(gca,'Ytick',[-89 -60 -30 0 30 60 88])
set(gca,'Yticklabels',{'90^oS','60^oS','30^oS','0','30^oN','60^oN','90^oN'},'fontsize',14)

ylabel('Latitude','fontsize',14)
xlabel('Marine Biological Richness')
set(gca,'fontsize',14)
xlim([0 1.5])
ylim([-75 88])
box on

clearvars -except output high low hist obs es8mn em8mn w8mn_save lat
%% Panel 2D
figure(7)
set(gcf,'units','inches','position',[2,2,4,5.5])

plot(high.ext_lat_2300,lat,'color',[.1 .8 .1],'linewidth',3)
hold on

ext_min = high.ext_min_lat_2300;
ext_max = high.ext_max_lat_2300;
patch_x = [lat(~isnan(ext_min)) fliplr(lat(~isnan(ext_min)))];
patch_y = [ext_min(~isnan(ext_min)) fliplr(ext_max(~isnan(ext_min)))];
p=patch(patch_y,patch_x,'g');set(p,'facecolor',[0.1 .8 0.1],'edgecolor',[0.1 .8 0.1],'facealpha',.3,'edgealpha',1);hold on
plot(high.ext_lat_2300,lat,'color',[.1 .8 .1],'linewidth',3)

box on
ylim([-75 88])
xlim([0 100])
ylabel('Latitude','fontsize',14)
xlabel('Extinction (%)','fontsize',14)

set(gca,'fontsize',14)
set(gca,'Ytick',[-89 -60 -30 0 30 60 88])
set(gca,'Yticklabels',{'90^oS','60^oS','30^oS','0','30^oN','60^oN','90^oN'},'fontsize',14)

clearvars -except output high low hist obs es8mn em8mn w8mn_save
%%  Figure 3 

% future time
ens.time=output.time_extended;

% Ensemble average
% High emissions
r8.ext_loss = high.ext_loss;
r8.ext_lossgain = high.ext_lossgain;
r8.ext_loss(1:9,5) = nan; % remove duplicate time point
r8.ext_lossgain(1:9,5) = nan; % remove duplicate time point
sixme = (r8.ext_loss+r8.ext_lossgain)/2;% average dispersal scenarios
sixme = nanmean(sixme);% average models
six = [r8.ext_lossgain;r8.ext_loss];
sixhigh = sixme+nanstd(six);% mean + standard deviation
sixlow = sixme-nanstd(six); % mean - standard deviaion


% Low emissions
r2.ext_loss = zeros(10,9)*nan;
r2.ext_lossgain = zeros(10,9)*nan;
r2.ext_loss(1:8,1:5) = low.ext_loss(1:8,:);
r2.ext_lossgain(1:8,1:5) = low.ext_lossgain(1:8,:);
r2.ext_loss(1:8,5) = nan; % remove duplicate time point
r2.ext_lossgain(1:8,5) = nan; % remove duplicate time point
r2.ext_loss(9,:) = output.low.ext_loss_2300(1,:);% add 2300 data
r2.ext_loss(10,:) = output.low.ext_loss_2300(2,:);% add 2300 data
r2.ext_lossgain(9,:) = output.low.ext_lossgain_2300(1,:);% add 2300 data
r2.ext_lossgain(10,:) = output.low.ext_lossgain_2300(2,:);% add 2300 data

sixme2 = (r2.ext_loss+r2.ext_lossgain)/2;% average dispersal
sixme2 = nanmean(sixme2);% average models
six2 = [r2.ext_lossgain;r2.ext_loss];
sixhigh2 = sixme2+nanstd(six2);% mean + standard deviation
sixlow2 = sixme2-nanstd(six2); % mean - standard deviaion

 
% create ∆temperature axis from Fig. 1 data
warm_interp = 0:1:20;% warming query points
sixme_taxis = (es8mn+em8mn)./2;% average extinction-∆T relationship across dispersal scenarios (Fig. 1)
sixme_interp2 = interp1(w8mn_save,sixme_taxis,warm_interp);% interpolate relationship across ∆T

% retrieve extinction levels at ∆T = 0,1,3,5,10,15 deg C 
% warming
widx2 = [warm_interp(1) warm_interp(2) warm_interp(4) warm_interp(6) warm_interp(11) warm_interp(16)];
% extinction
sidx2 = [sixme_interp2(1) sixme_interp2(2) sixme_interp2(4) sixme_interp2(6) sixme_interp2(11) sixme_interp2(16)];


% Sepkoski fossil record data from Rohde and Muller, Cycles in Fossil Diversity,Nature (2005)
time = obs.time; % Ma before present
gen1 = obs.gen1; % richness from all marine genera time ranges


% PBDB fossil record data from Alroy et al., Phanerozoic trends in the global diversity of marine invertebrates, Science (2008)
time2 = obs.time2;% Ma before present
pbdb = obs.pbdb; % richness from PBDB occurrences (w/standardized sampling)

% Normalize Sepkoski time-mean to PBDB time-mean
nconst = nanmean(gen1(:))./nanmean(pbdb(:));% normalization constant 
genn= gen1./nconst;% convert mean of sepkoski to mean of pbdb 

% fit 3rd order polynomial fit to Sepkoski constrained to earliest Cambrian levels 
psep = polyfix(time(~isnan(genn)),genn(~isnan(genn)),3,time(end-37),genn(end-37));
% secular trend over Phanerozoic from Sepkoski
modsep = psep(1).*time.^3 + psep(2).*time.^2 + psep(3).*time+psep(4);

% 3rd order polynomial fit to PBDB constrained to earliest Cambrian levels in Sepkoski
ppbdb = polyfix(time2(~isnan(pbdb)),pbdb(~isnan(pbdb)),3,time(end-37),genn(end-37));
% secular trend over Phanerozoic from PBDB
modpbdb = ppbdb(1).*time.^3 + ppbdb(2).*time.^2 + ppbdb(3).*time+ppbdb(4);

% detrend Sepkoski and add PBDB secular trend
gencor = genn-modsep+modpbdb; 

% Normalize to the recent
gencor= gencor./gencor(1);% Sepkoski data (w/PBDB trend) normalize to the Recent
gen1n = gen1./gen1(1);% Raw Sepkoski data normalized to the Recent

% Initial richness for future richness projections (Ro)
Ro = 1;

% Plot Fig. 3
figure
set(gcf,'units','inches','position',[2,2,8,4])

% normalization constant for future timescale
tfac = 10^1; % years

% Future patch
patch_x = [0 0 40 40];
patch_y = [0 1.2 1.2 0];
p=patch(patch_x,patch_y,'g');set(p,'facecolor',[1 1 0],'edgecolor',[.8 .8 0],'facealpha',.15,'edgealpha',0);hold on

% Past patch
patch_x = [-550 -550 0 0];
patch_y = [0 1.2 1.2 0];
p=patch(patch_x,patch_y,'g');set(p,'facecolor',[0 1 1],'edgecolor',[.8 .8 0],'facealpha',.1,'edgealpha',0);hold on

% Fossil data patch
patch_x = [fliplr(-time(~isnan(gen1n))') (-time(~isnan(gen1n))')];
patch_y = [fliplr((gen1n(~isnan(gen1n)))') (gencor(~isnan(gen1n)))'];
p=patch(patch_x,patch_y,'g');set(p,'facecolor',[.7 .7 .7],'edgecolor',[.7 .7 .7],'facealpha',.3,'edgealpha',0);hold on

% Plot Fossil data 
plot(-time,gencor,'color',[.7 .7 .7],'linewidth',2)
hold on
plot(-time,gen1n,'color',[.7 .7 .7],'linewidth',2)

% Future richness from extinction projections: Richness(t) = Ro*(1-Extinction)

% High emissions 
R_sixhigh = Ro.*(1-(sixhigh./100));% mean extinction+SD 
R_sixlow = Ro.*(1-(sixlow./100));% mean extinction-SD
R_six = Ro.*(1-(sixme./100));% mean extinction

patch_x = [(ens.time-1900)./tfac fliplr((ens.time-1900)./tfac)];
patch_y = [(R_sixhigh) fliplr(R_sixlow)];
p=patch(patch_x,patch_y,'g');set(p,'facecolor',[0.8 0 0],'edgecolor',[.8 0 0],'facealpha',.3,'edgealpha',0);hold on
plot((ens.time-1900)./tfac,R_six,'color',[.8 0 0],'linewidth',3) 

% Low emissions 
R_sixhigh2 = Ro.*(1-(sixhigh2./100));% mean extinction+SD
R_sixlow2 = Ro.*(1-(sixlow2./100));% mean extinction-SD
R_six2 = Ro.*(1-(sixme2./100));% mean extinction

patch_x = [(ens.time-1900)./tfac fliplr((ens.time-1900)./tfac)];
patch_y = [(R_sixhigh2) fliplr(R_sixlow2)];
p=patch(patch_x,patch_y,'g');set(p,'facecolor',[0 0 0.8],'edgecolor',[0 0 .8],'facealpha',.3,'edgealpha',0);hold on
plot((ens.time-1900)./tfac,R_six2,'color',[0 0 .8],'linewidth',3) 


% Delineate the Big 5 extinctions 
% Cretaceous
xx = [-65.5 -65.5];
yy = [0 .4];
plot(xx,yy,'-.','color',[.7 .7 .7])

% Triassic
xx = [-199.6 -199.6];
yy = [0 .42];
plot(xx,yy,'-.','color',[.7 .7 .7])

% Permian 
xx = [-251 -251];
yy = [0 .32];
plot(xx,yy,'-.','color',[.7 .7 .7])

% Devonian
xx = [-374 -374];
yy = [0 .4];
plot(xx,yy,'-.','color',[.7 .7 .7])

% Ordovician
xx = [-443.7 -443.7];
yy = [0 .45];
plot(xx,yy,'-.','color',[.7 .7 .7])

% Make Geologic timeline
xx = [-550 ((ens.time(end)-1900)./tfac)+0.1*((ens.time(end)-1900)./tfac)];
yy = [-.1 -.1];
plot(xx,yy,'k')
yy = [0 0];
plot(xx,yy,'k')

% Neogene (23 - 1.8 Ma)
xx = [-23 -23];
yy = [-.1 0];
plot(xx,yy,'k')
xx = [-1.8 -1.8];
yy = [-.1 0];
plot(xx,yy,'k')
text(-20,-.05,'N','fontsize',14)

% Quaternery (1.8 - Present)
xx = [0 0];
yy = [-.1 0];
plot(xx,yy,'k')

% Neogene (23 - 1.8 Ma)
xx = [-23 -23];
yy = [-.1 0];
plot(xx,yy,'k')
xx = [-1.8 -1.8];
yy = [-.1 0];
plot(xx,yy,'k')
text(-20,-.05,'N','fontsize',14)

% Pliocene (65.5 - 23 Ma)
xx = [-65.5 -65.5];
yy = [-.1 0];
plot(xx,yy,'k')
text(-55,-.05,'Pg','fontsize',14)

% Cretaceous (145.5 - 65.5 Ma)
xx = [-145.5 -145.5];
yy = [-.1 0];
plot(xx,yy,'k')
text(-115,-.05,'K','fontsize',14)

% Jurrasic(199.6 - 145.5 Ma)
xx = [-199.6 -199.6];
yy = [-.1 0];
plot(xx,yy,'k')
text(-175,-.05,'J','fontsize',14)

% Triassic (251 - 199.6 Ma)
xx = [-251 -251];
yy = [-.1 0];
plot(xx,yy,'k')
text(-235,-.05,'Tr','fontsize',14)

% Permian (299 - 251 Ma)
xx = [-299 -299];
yy = [-.1 0];
plot(xx,yy,'k')
text(-280,-.05,'P','fontsize',14)

% Carboniferous (359.2 - 299 Ma)
xx = [-359.2 -359.2];
yy = [-.1 0];
plot(xx,yy,'k')
text(-335,-.05,'C','fontsize',14)

% Devonian (416 - 359.2  Ma)
xx = [-416 -416];
yy = [-.1 0];
plot(xx,yy,'k')
text(-395,-.05,'D','fontsize',14)


% Silurian (443.7 - 416  Ma)
xx = [-443.7 -443.7];
yy = [-.1 0];
plot(xx,yy,'k')
text(-437,-.05,'S','fontsize',14)

% Ordivician (488 - 443.7  Ma)
xx = [-488 -488];
yy = [-.1 0];
plot(xx,yy,'k')
text(-472,-.05,'O','fontsize',14)

% Cambrian (542 - 488 Ma)
xx = [-542 -542];
yy = [-.1 0];
plot(xx,yy,'k')
text(-530,-.05,'Cm','fontsize',14)

% Anthropocene 
text(1,-.05,'Anthropocene','fontsize',10.7)

% Right-axis ∆Temperature scale 
% Convert extinction to richness for Temperature scale
waxis2 = Ro.*(1-(sidx2./100));

% Make tick mark for each ∆T
for ii = 1:length(waxis2)
    xx = [40-7,40];
    yy = [waxis2(ii) waxis2(ii)];
    plot(xx,yy,'color',[.7 .7 .7],'linewidth',1)
    
    % deltaT labels
    if ii == 1
        text(xx(2)+7,yy(2),'0')
    elseif ii == 2
        text(xx(2)+7,yy(2),'1')
    elseif ii == 3
        text(xx(2)+7,yy(2),'3')
    elseif ii == 4
        text(xx(2)+7,yy(2),'5')
    elseif ii == 5
        text(xx(2)+7,yy(2),'10')
    else
        text(xx(2)+7,yy(2),'15')
    end
end

xlimmax= ((ens.time(end)-1900)./tfac)+0.02*((ens.time(end)-1900)./tfac);
xlim([-542 xlimmax])
ylim([-.1 1.2])
box on
set(gca,'fontsize',14)
set(gca,'xticklabel',[500,400,300,200,100,0])
set(gca,'ytick',[0,0.25,.5,.75,1])
ylabel('Marine Biological Richness ')
xlabel('Past (millions of years ago')
text(xlimmax,-.17,'2300','fontsize',14)
text(0,-.275,'Future','fontsize',14)
text(75,1,'∆ Global Temperature','fontsize',14,'rotation',-90)

