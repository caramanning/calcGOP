function [d17O_H2O, d18O_H2O] = d17Od18OH2Ofrom2endmembers(SP,d18O_H2O_a,d18O_H2O_b,SPa,SPb,varargin)
% Estimate d17O-H2O and d18O-H2O for samples with measured salinity but
% unmeasured triple isotopic composition, using two endmembers with known
% and unique d18O-H2O and salinity. The 17O-excess for the two endmembers
% is an optional argument.
%
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% [d17O_H2O,d18O_H2O] = d17Od18OH2Ofrom2endmembers(SP,d18O_H2O_a,d18O_H2O_b,SPa,SPb)
% [d17O_H2O, d18O_H2O] = d17Od18OH2Ofrom2endmembers(10,-1.7,-9.3,30.7,0)
%   > d17O_H2O = -3.5902
%   > d18O_H2O = -6.8244
%
% [d17Oc,d18Oc] = d17Od18OH2Ofrom2endmembers(SP,d18O_H2O_a,d18O_H2O_b,SPa,SPb,excess17a,excess17b)
% [d17Oc,d18Oc] = d17Od18OH2Ofrom2endmembers(10,-1.7,-9.3,30.7,0,-10,20)
%   > d17O_H2O = -3.6006
%   > d18O_H2O = -6.8244
% 
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% SP          = practical salinity of sample with unknown d18O and d17O (PSS)
% d18O_H2O_a  = delta18O of H2O in endmember A (per mil vs VSMOW)
% d18O_H2O_b  = delta18O of H2O in endmember B (per mil vs VSMOW)
% SPa         = practical salinity of endmember A (PSS)
% SPb         = practical salinity of endmember B (PSS)
% varargin    = optional arguments
%         excess17Oa = 17O-excess of H2O in endmember A (per meg vs VSMOW)
%         excess17Ob = 17O-excess of H2O in endmember B (per meg vs VSMOW)
%
% SP can be size 1x1, 1xn or nx1.
% d18O_H2O_a, d18O_H2O_b SPa, SPb, excess17Oa, and excess17Ob must all be
% size 1x1.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% d17O_H2O   = delta17O of H2O for sample with salinity SP (per mil vs VSMOW)
% d18O_H2O   = delta18O of H2O for sample with salinity SP (per mil vs VSMOW)
%
% Note: After calculating d17O-H2O and d18O-H2O with respect to VSMOW using
% this function, you can use convertVSMOWtoair to calculate d17O-H2O and
% d18O-H2O with respect to atmospheric O2 for use in calcGOP_P11 and
% calcGOP_P11_NSS. See calcGOP_demo.m.
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% Estimates d17O and d18O for water samples as a function of salinity based
% on a linear combination of two endmembers with known d18O and salinity.
%
% d18O_H2O  = (R18_s/R18_VSMOW - 1) * 1000 
%           where R18_s is the 18O/16O ratio of the sample and R18_VSMOW
%           is the 18O/16O ratio of VSMOW.
%
% excess17O = (log(d17O_H2O/1000+1) - lambda_w*log(d18O_H2O/1000+1))*10^6 
%           where log is the natural logarithm.
%           Note that excess17O is typically called the 17O-excess (or
%           sometimes Delta17O) in the literature.
%
% lambda_w  = the slope on a plot of log(d17O_H2O/1000+1) vs log(d18O_H2O/1000+1) 
%           for meteoric water or for seawater. The slope results from
%           water isotopic fractionation from multiple environmental
%           processes and is close to the value of lambda_w for equilibrium
%           isotopic fractionation of water (0.529 +/- 0.001).
%
% NOTES:
% 1) Meijer and Li (1998), Barkan and Luz (2005), and Luz and Barkan (2010)
% find lambda_w = 0.528 for fresh and saltwater.
% 
% 2) Luz and Barkan (2010) plot
%       log(d17O_H2O/1000+1) = lambda_w*log(d18O_H2O/1000+1) + excess17O/10^6
% separately for globally-distributed meteoric water and seawater samples.
% The y-intercept on these plots is the average 17O-excess. They find 
%       excess17O = 33 per meg vs VSMOW for meteoric water 
%       excess17O = -5 per meg vs VSMOW for seawater
%
% In this function, the 17O-excess for both endmembers can be provided as
% an optional input. If not provided, it is set to 33 per meg for SP < 15
% and -5 per meg for SP >= 15 (i.e., the meteoric and seawater values from
% Luz and Barkan (2010), above).
%
% Note that Li et al. (2015) find for the continental USA the average
% 17O-excess of tap water (similar to meteoric) is 17(11) per meg, rather
% than 33 per meg. Future studies will help to refine our understanding of
% variability in 17O-excess and therefore the default 17O values in this
% function may be revised in the future.
%
% The results are much more sensitive to the choices of d18O-H2O and
% salinity for the endmembers than to the 17O-excess of the endmembers,
% since the 17O-excess generally varies by <0.1 per mil or less between
% samples whereas d18O-H2O varies globally by 40 per mil or more (Luz and
% Barkan 2010).
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
% DATA SOURCES FOR MEASUREMENTS OF d18O-H2O and 17O-excess
% -------------------------------------------------------------------------
%
% -----------------------
% d18O-H2O and 17O-excess
% -----------------------
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
% -----------------------------------------------
% d18O-H2O (and salinity for non-meteoric waters)
% -----------------------------------------------
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

% Check whether the 17O-excess for endmembers A and B was provided as an
% input. If not provided, set 17O-excess based on salinity following
% average published values for meteoric water and seawater.
% --- See description note 2 in header
if nargin > 5
    excess17Oa = varargin{1};
    excess17Ob = varargin{2};
else
    if SPa < 15
        excess17Oa = 33; % meteroric water value
    else
        excess17Oa = -5; % seawater value
    end;
    if SPb < 15
        excess17Ob = 33; % meteoric water value
    else
        excess17Ob = -5; % seawater value
    end;
end

% Calculate the d17O for endmembers A and B
% In Matlab, log gives the natural logarithm.
d17O_H2O_a = (exp(excess17Oa./10^6 + lambda_w*log(d18O_H2O_a./1000+1))-1).*1000;
d17O_H2O_b = (exp(excess17Ob./10^6 + lambda_w*log(d18O_H2O_b./1000+1))-1).*1000;

% Calculate the linear relationship for the two d18O and SP endmembers
% d18O = m*SP + b
m18 = (d18O_H2O_a - d18O_H2O_b)./(SPa - SPb);
b18 = d18O_H2O_a - m18.*SPa;

% Calculate the linear relationship for the two d17O and SP endmembers
% d17O = m*SP + b
m17 = (d17O_H2O_a - d17O_H2O_b)./(SPa - SPb);
b17 = d17O_H2O_a - m17.*SPa;

% Calculate d18O and d17O for the samples
d18O_H2O = m18.*SP + b18;
d17O_H2O = m17.*SP + b17;


end