%%% BIOLUMINESCENCE PROXIES PROCESSING %%%
% 
% Reference: Messié, M., I. Shulman, S. Martini and S.D.H. Haddock (2019). 
% Using fluorescence and bioluminescence sensors to characterize auto- and heterotrophic plankton communities. 
% Progress in Oceanography, 171, 76-92, doi:10.1016/j.pocean.2018.12.010.


load('example_datasets.mat')


%% BIOLUMINESCENCE PROCESSING
% proxies: dinoflagellates, larvaceans, copepods, small jellies
% Requires the function peakfinder.m, to download at https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder
%
% [BP,proxies]=bl_proxies_biolum(BP,window,flash_threshold,envelope_mini,varargin)
%
% demo corresponding to Fig. 3a:

[BP,proxies]=bl_proxies_biolum(BP);

flash_threshold=1E11;			% needs to match input to bl_proxies_biolum (this is the default value).
xBP = (BP.time-BP.time(1))*3600*24;												% elapsed time, in seconds
xP = (proxies.time-proxies.time(1))*3600*24;									% elapsed time, in seconds

% bldemos_fig1
figure
pos=get(gcf,'Position'); pos(3)=pos(3)*2; set(gcf,'Position',pos)				% make the figure larger
axes('Position',[0.07 0.55 0.7 0.4],'XLim',[0 max(xBP)]), hold on
	h1=plot(xBP,BP.biolum,'k');													% original time series
	h2=plot(xBP,BP.med_bgrd,'b','LineWidth',1.5);								% median background
	h3=plot(xBP,BP.min_bgrd,'c'); plot(xBP,BP.max_bgrd,'c')						% envelope
	ilowflash = BP.iflash & BP.biolum-BP.med_bgrd <= flash_threshold;			% time steps with low-intensity flashes
	ihighflash = BP.iflash & BP.biolum-BP.med_bgrd > flash_threshold;			% time steps with high-intensity flashes
	h4=plot(xBP(ilowflash),BP.biolum(ilowflash),'Color',[0.75 0.75 0],'LineStyle','none','marker','p');	% low flashes
	h5=plot(xBP(ihighflash),BP.biolum(ihighflash),'Color','r','LineStyle','none','marker','p');			% high flashes
	ylabel('Bioluminescence (ph s^{-1})')
	title('Example of bioluminescence proxies processing')
	legend([h1,h2,h3,h4,h5],{'bioluminescence','med-background','envelope','low-intensity flashes','high-intensity flashes'},...
		'Position',[0.79 0.7 0.2 0.15],'FontSize',9)
	% set(gca,'XLim',[50.5 110.5])		% to zoom on the part displayed in Fig. 2a
axes('Position',[0.07 0.1 0.7 0.4],'XLim',[0 max(xBP)],'YLim',[0 1]), hold on
	h1=plot(xP,proxies.dinoflagellates./max(proxies.dinoflagellates),'c','LineWidth',1.5);		% dinoflagellates proxy
	h2=plot(xP,proxies.larvaceans./max(proxies.larvaceans),'Color',[0.75 0.75 0],'LineWidth',1.5);	% larvacean proxy
	h3=plot(xP,proxies.copepods./max(proxies.copepods),'r','LineWidth',1.5);					% copepods proxy
	h4=plot(xP,proxies.jellies./max(proxies.jellies),'Color',[0.9 0.5 0.9],'LineWidth',1.5);	% small jellies proxy
	ylabel('Normalized proxies')
	legend([h1,h2,h3,h4],{'dinoflagellates','larvaceans','copepods','jellies'},'Position',[0.79 0.25 0.2 0.15],'FontSize',9)
xlabel('Elapsed time (seconds)')




%% FLUORESCENCE/BIOLUMINESCENCE PROCESSING
% proxies: autotrophic dinoflagellates (adinos), heterotrophic dinoflagellates (hdinos), other phytoplankton (aother)
%
% proxies=bl_proxies_fluobiolum(fluo,bgrd_BL,ratioAdinos,calfactor)
%
% demo corresponding to Fig. 7 (ratioAdinos and calfactor were computed from the entire AUV dataset)
% The colorbar used in the paper is "dense" from the cmocean toolbox (https://matplotlib.org/cmocean/)

ratioAdinos=1.26E13;
calfactor=1.37E-3;
proxies=bl_proxies_fluobiolum(AUV.fluo,AUV.bgrd_BL,ratioAdinos,calfactor);

% bldemos_fig2
figure
subplot(3,1,1)			% h-dinos proxy
	scatter(AUV.time,AUV.depth,5,proxies.hdinos,'fill')
	set(gca,'YDir','reverse','YLim',[0 60],'YTick',0:20:60,'XLim',[min(AUV.time) max(AUV.time)])
	datetick('x','keeplimits'), caxis([0 1])
	ylabel('Depth (m)'), title('Heterotrophic dinoflagellate proxy')
subplot(3,1,2)			% a-dinos proxy
	scatter(AUV.time,AUV.depth,5,proxies.adinos,'fill')
	set(gca,'YDir','reverse','YLim',[0 60],'YTick',0:20:60,'XLim',[min(AUV.time) max(AUV.time)])
	datetick('x','keeplimits'), caxis([0 1])
	ylabel('Depth (m)'), title('Autotrophic dinoflagellate proxy')
subplot(3,1,3)			% a-other proxy
	scatter(AUV.time,AUV.depth,5,proxies.aother,'fill')
	set(gca,'YDir','reverse','YLim',[0 60],'YTick',0:20:60,'XLim',[min(AUV.time) max(AUV.time)])
	datetick('x','keeplimits'), caxis([0 1])
	ylabel('Depth (m)'), title('Other phytoplankton proxy (diatoms)')
colorbar('Position',[0.93 0.15 0.01 0.75])
mat_color=colormap('gray'); colormap(flipud(mat_color(1:end-2,:)))

