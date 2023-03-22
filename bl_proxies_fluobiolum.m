function proxies=bl_proxies_fluobiolum(fluo,bgrd_BL,ratioAdinos,calfactor)

% BL_PROXIES_FLUOBIOLUM: calculates autotrophic dinoflagellates (adinos), heterotrophic dinoflagellates (hdinos), and other phytoplankton (aother) proxies 
% from a set of fluorescence / background bioluminescence data (such as calculated by bl_proxies_biolum).
% The proxies assume that:
%	- the phytoplankton community is dominated by dinoflagellates and other unrelated phytoplankton (often diatoms in the coastal ocean), 
%		such that adino and aother do not correlate on average
%	- the dinoflagellate community is characterized by a constant bioluminescence to fluorescence ratio, equal to ratioAdinos
%
% proxies=bl_proxies_fluobiolum(fluo,bgrd_BL,ratioAdinos,calfactor)
% 	For an example, see bl_demos.
% 
% INPUTS:
%	fluo: fluorescence (proxy for phytoplankton = adinos + aother)
%	bgrd_BL: background bioluminescence (proxy for dinoflagellates)
% 	ratioAdinos: typical bgrd_BL/fluo ratio for dinoflagellates populations, typically identified from an histogram over an entire dataset
%	calfactor: possible calibration to normalize the proxies (typically fluorescence 99th percentile). 
%		If not given, no calibration is applied (calfactor=1) and the proxies are given in fluorescence units.
%
% OUTPUTS:
%	proxies: structure containing the proxies (.aother, .adinos, .hdinos)
%
% Monique Messié, 2018, MBARI
% Reference: Messié, M., I. Shulman, S. Martini and S.D.H. Haddock (2019). 
% Using fluorescence and bioluminescence sensors to characterize auto- and heterotrophic plankton communities. 
% Progress in Oceanography, 171, 76-92, doi:10.1016/j.pocean.2018.12.010.


% Reading & checking input data
if nargin<4, calfactor=1; end
if nargin<3, error('Give fluo, bgrd_BL and ratioAdinos'), end
if ~min(size(fluo)==size(bgrd_BL)), error('fluo and bgrd_BL must have the same size'), end

% Initialize proxies output
proxies=struct('hdinos',fluo*NaN,'adinos',fluo*NaN,'aother',fluo*NaN);
fluo_dinos = bgrd_BL/ratioAdinos;		% convert bgrd_BL in fluorescence units with the ratioAdinos conversion

% Autotrophic dinoflagellates: adinos = min(fluo, fluo_dinos)
proxies.adinos=min(fluo,fluo_dinos);

% Heterotrophic dinoflagellates: hdinos = remaining fluo_dinos after removing adinos. Will be 0 if fluo_dinos<fluo.
proxies.hdinos = fluo_dinos-proxies.adinos; 

% Other autotrophic plankton (such as diatoms): aother = remaining fluo after removing adinos. Will be 0 if fluo<fluo_dinos.
proxies.aother = fluo-proxies.adinos; 

% Proxies clean up and normalization
for species={'aother','adinos','hdinos'}, species=species{:}; 
	proxies.(species)(isnan(fluo) | isnan(bgrd_BL))=NaN; 
	proxies.(species)=proxies.(species)/calfactor; 
end


return


