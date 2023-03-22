function varargout=bl_proxies_biolum(BP,window,flash_threshold,envelope_mini,window_proxies)

% BL_PROXIES_BIOLUM: calculates the bioluminescence background and flashes, and the corresponding proxies. 
% The method is intended for 60Hz data and requires the function peakfinder.m (https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder).
% For fast computing, the window smoothing method (bl_window_smoothing.m) assumes that the data is equally spaced. This is not always the case 
% (e.g., a few missing data points near the beginning of the examples time series) and for larger data gaps, a different method should be used 
% such as a loop (much slower). Note that the fast method uses a lot of RAM which can be an issue for very long time series.
%
% The program:
% 1- computes med-background (med_bgrd here) by median window-averaging (median), representing the mean dinoflagellates bioluminescence
% 2- computes min-background (here min_bgrd) similarly but with a minimum instead of median function,
% 3- calculates the envelope (range of variation of dinoflagellates bioluminescence) by symmetry across med-background
% 4- finds flashes in the time series and keeps those above the envelope
% 5- generates 1Hz proxies time series
%
% bl_proxies_biolum(BP,window,flash_threshold,envelope_mini,window_proxies)
% 	BP = bl_proxies_biolum(BP); is identical to BP = bl_proxies_biolum(BP,5,1E11,1.5E10,15);
% 	For an example, see bl_demos.
%
% INPUTS: 
% 	BP: structure containing the bioluminescence time series (.time, .biolum, .flow). 
%		Time in matlab units, biolum in ph/s, flow in L/s (assumed to be 0.35 if missing). 
%	window: time window (in s) used for window-sliding (default 5 s)
%	flash_threshold: threshold separating low and high intensity flashes (default 1E11 ph/s)
%	envelope_mini: minimum value for the envelope (max_bgrd - med_bgrd) to avoid very dim flashes when the background is low (default 1.5E10 ph/s)
%	window_proxies: window over which the zooplankton proxies are computed (nb of flashes, max flash intensity) (default 15 s)
%
% OUTPUTS: 
% BP as the input structure with additional fields:
% 	.med_bgrd: median background (med-background), in biolum units (ph/s).
% 	.min_bgrd: minimum background (min-background), in biolum units (ph/s).
%	.max_bgrd: upper limit of the envelope, obtained by symmetry as 2*med_bgrd-min_bgrd, in biolum units (ph/s)
%		The envelope is delimited by min_bgrd and max_bgrd. 
%	.iflash: indices for flashes above the envelope (both low- and high-intensity)
% proxies as a new 1Hz structure
%	.time				1Hz time (matlab format)
% 	.dinoflagellate		dinoflagellate proxy (min-background expressed in ph/L)
%	.larvaceans			larvaceans proxy (low-intensity flashes within a window_proxies window in flashes/L)
%	.copepods			copepods proxy (high-intensity flashes within a window_proxies window in flashes/L)
%	.jellies			jellies proxy (maximum flash intensity within a window_proxies window in ph/s)
% 
% Monique Messié, 2018, MBARI
% Reference: Messié, M., I. Shulman, S. Martini and S.D.H. Haddock (2019). 
% Using fluorescence and bioluminescence sensors to characterize auto- and heterotrophic plankton communities. 
% Progress in Oceanography, 171, 76-92, doi:10.1016/j.pocean.2018.12.010.


% Reading input data
if nargin<5 || isempty(window_proxies), window_proxies=15; end
if nargin<4 || isempty(envelope_mini), envelope_mini=1.5E10; end
if nargin<3 || isempty(flash_threshold), flash_threshold=1E11; end
if nargin<2 || isempty(window), window=5; end
if nargin<1, error('Give BP'), end
if ~isfield(BP,'flow'), BP.flow = BP.time*0+0.35; end			% if flow is not given, assumes 350 mL/s

% Checking that input time steps are regular and 60Hz (within 1%), and issuing a warning if they are not
time_diff=BP.time(2:end)-BP.time(1:end-1); 
if abs(max(time_diff)-min(time_diff))/max(time_diff)>0.01
	disp('!!! WARNING !!! the window smoothing method will assume that time steps are regular!') 
end
time_interval=median(time_diff);
if abs(time_interval-1/60/3600/24)/time_interval>0.01
	disp('!!! WARNING !!! the program assumes time resolution is 60Hz!')
end


% -------------------------- Proxies processing -------------------------- %

% 1- Compute med-background (median window-sliding method)
BP.med_bgrd=bl_window_smoothing(BP.biolum,BP.time,window/3600/24,'median');		% median window smoothing (possible to use movmedian instead)
BP.med_bgrd=bl_window_smoothing(BP.med_bgrd,BP.time,window/3600/24,'mean');		% additional smoothing (possible to use movmean instead)

% 2- Compute min-background (minimum window-sliding method)
BP.min_bgrd=bl_window_smoothing(BP.biolum,BP.time,window/3600/24,'min'); 		% min window smoothing
BP.min_bgrd=bl_window_smoothing(BP.min_bgrd,BP.time,window/3600/24,'mean');		% additional smoothing (needed to avoid a box-like time series)

% 3- Calculate the envelope delimited by min_bgrd and symmetrical across med_bgrd
BP.max_bgrd=2*BP.med_bgrd-BP.min_bgrd;
ilow=BP.max_bgrd-BP.med_bgrd<envelope_mini; 
BP.max_bgrd(ilow)=BP.med_bgrd(ilow)+envelope_mini;

% 4- Identify flashes
iok=~isnan(BP.biolum); 										% only work on data points
indices=1:length(BP.biolum); indices=indices(iok);			% generate a corresponding set of indices marking time steps
iiflash1=peakfinder(BP.biolum(iok),1E8); 					% find all flashes with a minimum peak height of 1E8 ph/s (only on data points)
iiflash=indices(iiflash1); 									% get the corresponding indices in the original time series
BP.iflash = BP.time*0; BP.iflash(iiflash)=1; BP.iflash=logical(BP.iflash);	% keep flashes as a logical vector (0 = no flash, 1 = flash) 
BP.iflash(BP.biolum<BP.max_bgrd)=0; 						% remove flashes below the envelope (max_bgrd), expected to be generated by dinoflagellates

% 5- Generate proxies at a 1Hz time resolution
proxies=struct('time',(floor(min(BP.time*3600*24)):1:max(BP.time*3600*24))'/3600/24);		% 1Hz time series starting at the first second
[~,~,idt] = unique(floor(BP.time*3600*24)); 												% find positions for each 1Hz time step
proxies.dinoflagellates = accumarray(idt,BP.min_bgrd./BP.flow,[],@(x) mean(x,'omitnan'),NaN);			% dinos proxy = min_bgrd in ph/L (simple average on 1Hz time steps)
for varname={'larvaceans','copepods','jellies'}, proxies.(varname{:})=proxies.time*NaN; end	% initialize zoo proxies vectors
for itime=1:length(proxies.time)
	iBP = BP.time>=proxies.time(itime)-window_proxies/2/3600/24 ...
		& BP.time<=proxies.time(itime)+window_proxies/2/3600/24;							% find BP time indices within the window_proxies window
	volume_sampled=mean(BP.flow(iBP),'omitnan')*sum(iBP)/60;								% volume of water sampled during the time window (assumes BP is 60Hz)
	proxies.larvaceans(itime)=sum(BP.iflash(iBP & BP.biolum-BP.med_bgrd<=flash_threshold))/volume_sampled;	% larvaceans = number of low-intensity flashes per liter
	proxies.copepods(itime)=sum(BP.iflash(iBP & BP.biolum-BP.med_bgrd>flash_threshold))/volume_sampled;	% copepods = number of high-intensity flashes per liter
	proxies.jellies(itime)=max(BP.biolum(iBP)-BP.med_bgrd(iBP));							% jellies = max flash intensity within the window (ie max biolum relative to med_bgrd)
end


% ---------------------------------------------------------------------- %


% Outputs
varargout={BP,proxies}; varargout=varargout(1:nargout);

return
