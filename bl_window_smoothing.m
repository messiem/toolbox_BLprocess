function serie_out=bl_window_smoothing(serie,time,window,smooth_function)

% BL_WINDOW_SMOOTHING: replaces each point by its averaging (or min, or median) within a window around that point.
% For fast computing, the method assumes that the data is equally spaced. When this is not the case, the result is an approximation.
% For large data gaps, a different method should be used such as a loop (much slower). 
% Note that the method builds a window * length(serie) matrix which may use a lot of RAM for a long time series.
% This function can be replaced by Matlab's native functions movmean and movmedian for the "mean" and "median" cases.
%
% serie_out = bl_window_smoothing(serie,time,window,smooth_function)
% 	e.g., serie_out = bl_window_smoothing(biolum,time,5/3600/24,'median');
%
% INPUTS: 
% 	serie: time series to be processed.
%	time: corresponding time vector. Can be of any unit (e.g., a depth vector for processing a vertical profile) or left empty (assumes regular intervals)
%	window: in the same unit as time (days if using matlab time format, number of points if time is left empty). 
%		The window will be adapted such that the corresponding number of points is an odd number (for symmetry around the data point).
%	smooth_function: function to be applied within the window. Available options are 'mean','median','min','max'.
%
% OUTPUTS: 
% 	serie_out: new, processed time series.
% 
% Monique Messié, 2018, MBARI


% Reading input data and checking that the window is short enough relative to the data
if isempty(time), time=1:length(serie); end
if max(time)-min(time)<window, serie_out=serie*NaN; disp('!!! WARNING !!! The chosen window is larger than the time, no data!'), return, end

% Calculates the number of points within a half_window, such that the processing runs within a window of (2*nb_halfwindow + 1) data points (ie odd number).
nb_pts=length(serie);
time_interval=median(time(2:end)-time(1:end-1));		% uses a median in case there are a few jumps in the time series
nb_halfwindow = floor(window/2/time_interval);

% Generates a matrix where each row contains the data within the window centered on the corresponding data point 
% (beginning & end will be approximative with less data points than a full window).
mat_data=ones(nb_pts,2*nb_halfwindow+1)*NaN;
for ipts=1:nb_pts
	iwindow=max(1,ipts-nb_halfwindow):min(nb_pts,ipts+nb_halfwindow);		% indices of data points within the window (assuming regular intervals)
	mat_data(ipts,1:length(iwindow))=serie(iwindow);					% each row in mat_data is a window of data points
end

% Computes the result for all windows at once
eval(['f_mean=@',smooth_function,';'])
switch smooth_function
	case {'min','max'}, serie_out = f_mean(mat_data,[],2); 
	case {'median','mean'}, serie_out = f_mean(mat_data,2,'omitnan');
	otherwise, error('Undefined smooth function!')
end
serie_out=reshape(serie_out,size(serie));


return




