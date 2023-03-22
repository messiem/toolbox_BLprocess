# Toolbox BLprocess #

### Description ###

This toolbox contains the Matlab programs necessary to compute proxies from bioluminescence and fluorescence datasets as described in Messié et al. (2019), along with example datasets and demonstrations reproducing some figures from the paper. The processing assumes that the input bioluminescence time series is 60Hz and time steps regular. The toolbox requires the function peakfinder.m which can be downloaded at https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder.
The main two processing steps are:

#### bioluminescence processing ####

    bl_proxies_biolum

This function computes dinoflagellate, larvacean, copepod and small jelly proxies at a 1Hz resolution from 60Hz bioluminescence time series. In particular, it computes the background bioluminescence needed for the fluorescence/bioluminescence processing.

#### fluorescence/bioluminescence processing ####

    bl_proxies_fluobiolum

This function computes autotrophic dinoflagellate, heterotrophic dinoflagellate, and other phytoplankton (e.g., diatom) proxies from fluorescence and background bioluminescence time series. 

### Demos ###

See function: 
	
    bl_demos

This functions runs both the bl_proxies_biolum and bl_proxies_fluobiolum functions, computes proxies on example datasets and displays them, reproducting the figures bldemos_fig1.jpg (reproducing Fig. 3a in Messié et al. 2019) and bldemos_fig2.jpg (reproducing Fig. 7g,h,i in Messié et al. 2019).

* * *

### Reference ###

Please refer to this paper when using the toolbox: 

Messié, M., I. Shulman, S. Martini and S.H.D. Haddock (2019). **Using fluorescence and bioluminescence sensors to characterize auto- and heterotrophic plankton communities**. *Progress in Oceanography*, 171, 76-92, doi:10.1016/j.pocean.2018.12.010. 

(available online at https://www.sciencedirect.com/science/article/pii/S0079661118300478).

### Contact ###

monique@mbari.org



