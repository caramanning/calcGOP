function [d17O_H2O] = d17Ofromd18OH2O(d18O_H2O,excess17O)
% Calculate d17O-H2O from d18O-H2O and 17O-excess. All values are
% referenced to VSMOW.
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% [d17O_H2O] = d17Ofromd18OH2O(d18O_H2O,excess17O)
% [d17O_H2O] = d17Ofromd18OH2O(-2,-5)
%   > d17O_H2O = -1.0615
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% d18O_H2O    = delta18O of H2O (per mil vs VSMOW)
% excess17O   = 17O-excess of water (per meg vs VSMOW)
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% d17O_H2O   = delta17O of H2O (per mil vs VSMOW)
%
% Note: After calculating d17O_H2O vs VSMOW using this function, you can
% use convertVSMOWtoair to calculate d17O_H2O and d18O_H2O vs atmospheric
% O2 (air) for use in calcGOP_P11 and calcGOP_P11_NSS.
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% Calculates d17O-H2O from d18O-H2O and 17O-excess
%
% d18O_H2O  = (R18_s/R18_VSMOW - 1) * 1000 
%           where R18_s is the 18O/16O ratio of the sample and R18_VSMOW
%           is the 18O/16O ratio of VSMOW.
%
% excess17O = (log(d17O_H2O/1000+1) - lambda_w*log(d18O_H2O/1000+1))*10^6 
%           where log is the natural logarithm.
%           Note that excess17O is typically called the 17O-excess in the
%           literature.
%
% lambda_w  = the slope on a plot of log(d17O/1000+1) vs log(d18O/1000+1) 
%           for meteoric water or for seawater.
%
% NOTES:
% 1) Meijer and Li (1998), Barkan and Luz (2005), and Luz and Barkan (2010)
% find lambda_w = 0.528 for fresh and saltwater.
% 
% 2) You must provide an input for the 17O-excess. Here are some reasonable
% values:
% Luz and Barkan (2010) find for globally-distributed waters 
%       excess17O ~= 33 per meg vs VSMOW for meteoric water 
%       excess17O ~= -5 per meg vs VSMOW for seawater
% These excess17O values are the y-intercept from a plot of 
%       log(d17O_H2O/1000+1) = lambda_w*log(d18O_H2O/1000+1) + excess17O/10^6
% for meteoric water or for seawater.
%
% Li et al. find for the continental USA the average 17O-excess of tap
% water (similar to meteoric) is 17(11) per meg. Future studies will help
% to refine our understanding of variability in 17O-excess.
% 
% The delta17O-H2O and are much more sensitive to the d18O-H2O than to the
% 17O-excess, since the 17O-excess generally varies by <0.1 per mil or less
% between samples whereas d18O varies globally by 40 per mil or more (Luz
% and Barkan 2010).
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%
% Barkan, E., & Luz, B. (2005). High precision measurements of 17O/16O and
% 18O/16O ratios in H2O. Rapid Communications in Mass Spectrometry, 19(24),
% 3737-3742
%
% Li, S., Levin, N. E., & Chesson, L. A. (2015). Continental scale
% variation in 17O-excess of meteoric waters in the United States.
% Geochimica et Cosmochimica Acta, 164, 110-126.
%
% Luz, B. and Barkan, E. (2010). Variations of 17O/16O and 18O/16O in
% meteoric waters. Geochimica et Cosmochimica Acta, 74(22), 6276-6286.
%
% Meijer, H. A. J. and Li, W. J. (1998). The use of electrolysis for
% accurate d17O and d18O isotope measurements in water. Isotopes in
% Environmental and Health Studies, 34(4), 349-369.
%
% -------------------------------------------------------------------------
% DATA SOURCES FOR MEASUREMENTS OF d18O and 17O-excess
% -------------------------------------------------------------------------
%
% -------------------
% d18O and 17O-excess
% -------------------
%
% Landais A., Risi C., Bony S., Vimeux F., Descroix L., Falourd S. and
% Bouygues A. (2010) Combined measurements of 17O-excess and d-excess in
% African monsoon precipitation: Implications for evaluating convective
% parameterizations. Earth Planet Sci Lett. 298, 104-112.
%
% Landais A., Steen-Larsen H. C., Guillevic M., Masson-Delmotte V., Vinther
% B. and Winkler R. (2012a) Triple isotopic composition of oxygen in
% surface snow and water vapor at NEEM (Greenland). Geochim. Cosmochim.
% Acta 77, 304-316.
% 
% Landais A., Ekaykin A., Barkan E., Winkler R. and Luz B. (2012b) Seasonal
% variations of d-excess and 17O-excess in snow at the Vostok station (East
% Antarctica). J. Glaciol. 58, 725-733.
% 
% Li, S., Levin, N. E., & Chesson, L. A. (2015). Continental scale
% variation in 17O-excess of meteoric waters in the United States.
% Geochimica et Cosmochimica Acta, 164, 110-126.
%
% Luz, B. and Barkan, E. (2010). Variations of 17O/16O and 18O/16O in
% meteoric waters. Geochimica et Cosmochimica Acta, 74(22), 6276-6286.
%
%
% -------------------------------------------
% d18O (and salinity for non-meteoric waters)
% -------------------------------------------
%
% Bowen G. J. and Revenaugh J. (2003) Interpolating the isotopic
% composition of modern meteoric precipitation. Water Resources Research
% 39(10), 1299, doi:10.129/2003WR002086. 
%
% Bowen, G. J. (2016) The Online Isotopes in Precipitation Calculator,
% version 2.2. http://www.waterisotopes.org
% 
% LeGrande, A. N. and  Schmidt, G. A. (2006). Global gridded data set of the
% oxygen isotopic composition in seawater. Geophysical Research Letters,
% 33(12). Data set available at http://data.giss.nasa.gov/o18data/grid.html
%
% Schmidt, G.A., G. R. Bigg and E. J. Rohling. (1999). Global Seawater
% Oxygen-18 Database - v1.21. http://data.giss.nasa.gov/o18data/
% 
% IAEA GNIP and GNIR Databases (Global Networks of Isotopes in
% Precipitation and Rivers). Access through WISER (Water Isotope System for
% Data Analysis, Visualization and Electronic Retrieval). 
% http://www-naweb.iaea.org/napc/ih/IHS_resources_isohis.html
% http://nucleus.iaea.org/wiser (free signup required to access data)
%
%
% -------------------------------------------------------------------------
% SOFTWARE AUTHORSHIP, ATTRIBUTION, and URL:
% -------------------------------------------------------------------------
% Written by Cara C. Manning (cmanning@whoi.edu) and Evan M. Howard 
% Woods Hole Oceanographic Institution and Massachusetts Institute of
% Technology
%
% Cite as: CC Manning and EM Howard (2016) calcGOP: Functions for
% calculating gross oxygen production from measurements of the triple
% oxygen isotopic composition of dissolved O2. 
% http://github.com/caramanning/calcGOP/
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define slope for plot of log(d17O-H2O/1000+1) vs log(d18O-H2O/1000+1)
% --- See description note 1 in header
lambda_w = 0.528;

% Calculate the d17O-H2O
% --- See description note 2 in header
d17O_H2O = (exp(excess17O./10^6 + lambda_w*log(d18O_H2O/1000+1))-1).*1000;


end