% This file demonstrates how to run the calcGOP functions to determine GOP.
% -------------------------------------------------------------------------
% SOFTWARE AUTHORSHIP, ATTRIBUTION, and URL:
% -------------------------------------------------------------------------
% Written by Cara C. Manning (cmanning@whoi.edu) and Evan M. Howard 
% Woods Hole Oceanographic Institution and Massachusetts Institute of
% Technology
%
% Cite as: CC Manning and EM Howard (2017) calcGOP: Functions for
% calculating gross oxygen production from measurements of the triple
% oxygen isotopic composition of dissolved O2. 
% http://github.com/caramanning/calcGOP/
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License, which 
% is available at http://www.apache.org/licenses/LICENSE-2.0
% -------------------------------------------------------------------------

% Initialize with sample data and endmembers
% We have the following data for one sample:
d17O_O2 = 0.1773; % sample delta17O-O2 with respect to atmospheric O2 (air)
d18O_O2 = 0.3016; % sample delta18O-O2 with respect to atmospheric O2 (air)
pt = 6.29;        % potential temperature in degC
SP=13.66;         % practical salinity in PSS
kO2 = 0.8;        % gas transfer velocity in m d^-1

% Define two endmembers for the local isotopic composition of water
% freshwater d18O-H2O with respect to VSMOW
d18O_H2O_fw_VSMOW = -9.3; 

% seawater d18O-H2O with respect to VSMOW
d18O_H2O_sw_VSMOW = -1.68; 

% freshwater and seawater salinity
SPfw = 0;
SPsw = 31.25;
% ----------

% Calculate d17O-H2O and d18O-H2O in per mil vs VSMOW for water with a
% known salinity, as a function of two d18O-H2O and salinity endmembers
[d17O_H2O_VSMOW, d18O_H2O_VSMOW] = d17Od18OH2Ofrom2endmembers(SP,d18O_H2O_fw_VSMOW,d18O_H2O_sw_VSMOW,SPfw,SPsw);

% Convert d17O-H2O and d18O-H2O vs VSMOW 
% to d17O-H2O and d18O-H2O vs atmospheric O2 (air)
[d17O_H2O_air,d18O_H2O_air] = convertVSMOWtoair(d17O_H2O_VSMOW,d18O_H2O_VSMOW);
% ----------

% Calculate GOP in mol O2 m^-2 d^-1 using the local water isotopic
% composition
GOP = calcGOP_P11(d17O_O2,d18O_O2,SP,pt,kO2,d17O_H2O_air,d18O_H2O_air);

% Now calculate GOP assuming the isotopic composition of water is
% equivalent to VSMOW
GOP_withVSMOW = calcGOP_P11(d17O_O2,d18O_O2,SP,pt,kO2);

% Calculate ratio of GOP with local water isotopic composition to GOP with
% VSMOW. If functions are working correctly GOPratio = 1.5924.
GOPratio = GOP/GOP_withVSMOW;

