function [d17O_air,d18O_air] = convertVSMOWtoair(d17O_VSMOW,d18O_VSMOW)
% Function to convert oxygen isotopic composition from d17O and d18O with
% respect to VSMOW to d17O and d18O with respect to air. Inputs and outputs
% are in per mil.
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% [d17O_air,d18O_air] = convertVSMOWtoair(d17O_VSMOW,d18O_VSMOW)
% [d17O_air,d18O_air] = convertVSMOWtoair(-5.884,-11.168)
%   > d17O_air = -17.6971
%   > d18O_air = -34.2315
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% d17O_VSMOW      = delta17O (per mil vs VSMOW)
% d18O_VSMOW      = delta18O (per mil vs VSMOW)
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% d17O_air      = delta17O (per mil vs atmospheric O2)
% d18O_air      = delta18O (per mil vs atmospheric O2)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% Converts oxygen isotopic composition from d17O and d18O with respect to
% VSMOW to d17O and d18O with respect to air. Inputs and outputs are in per
% mil.
%
% d18O_VSMOW     = (R18_samp/R18_VSMOW - 1) * 1000 
% d18O_air       = (R18_samp/R18_air - 1) * 1000 
% where R18_samp is the 18O/16O ratio of the sample, R18_VSMOW is the
% 18O/16O ratio of VSMOW, and R18_air is the ratio in atmospheric O2.
%
% This function uses the isotopic composition of VSMOW with respect to
% atmospheric O2 as defined in Barkan and Luz (2011) Rapid Comm. Mass
% Spectrom., Table 1 average 2007-2009.
%
% The oxygen isotopic composition of VSMOW is given in an IUPAC technical 
% report as 
%       r(18O/16O) = 0.00200520(45) 
%       r(17O/16O) = 0.0003799(8).
% See De Laeter et al. (2003), Pure Appl. Chem. 75(6), page 739
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%
% Barkan, E. and B. Luz (2011). The relationships among the three stable
% isotopes of oxygen in air, seawater and marine photosynthesis. Rapid
% communications in mass spectrometry: RCM, 25(16), 2367-2369.
%
% De Laeter, J. R., J. K. Bohlke, P. De Bievre et al. (2003). Atomic
% weights of the elements: Review 2000. Pure Appl. Chem., 75(6), 683-800.
%
% -------------------------------------------------------------------------
% SOFTWARE AUTHORSHIP, ATTRIBUTION, and URL:
% -------------------------------------------------------------------------
% Written by Cara C. Manning (cmanning@whoi.edu), Woods Hole Oceanographic
% Institution and Massachusetts Institute of Technology
%
% Cite as: CC Manning (2016) calcGOP: Functions for calculating gross
% oxygen production from measurements of the triple oxygen isotopic
% composition of dissolved O2. Zenodo. doi: 10.5281/zenodo.59268
%
% Find the latest version of the functions at 
% https://github.com/caramanning/calcGOP/
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%delta18O and delta17O of VSMOW in per mil vs air
d18O_VSMOW_air   = -23.324; 
d17O_VSMOW_air   = -11.883;

% Ratio of 18O/16O and 17O/16O of VSMOW
R18_VSMOW = 0.0020052;
R17_VSMOW = 0.0003799;

% Ratio of 18O/16O and 17O/16O of air
R18_air = R18_VSMOW./(d18O_VSMOW_air/1000+1);
R17_air = R17_VSMOW./(d17O_VSMOW_air/1000+1);

R18 = (d18O_VSMOW/1000+1)*R18_VSMOW;
R17 = (d17O_VSMOW/1000+1)*R17_VSMOW;

d18O_air = (R18./R18_air - 1) * 1000;
d17O_air = (R17./R17_air - 1) * 1000;
end