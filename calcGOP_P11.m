function GOP = calcGOP_P11(d17O_O2,d18O_O2,SP,pt,kO2,varargin)
% calcGOP_P11: steady state gross oxygen production using Prokopenko 2011
% steady state equation 7
% -------------------------------------------------------------------------
% USAGE:
% -------------------------------------------------------------------------
% GOP = calcGOP_P11(d17O_O2,d18O_O2,SP,pt,kO2)
% GOP = calcGOP_P11(-0.571224,-1.27791,35,15,1.4)
%   > GOP = 0.1914
% or
% GOP = calcGOP_P11(d17O_O2,d18O_O2,SP,pt,kO2,d17O_H2O,d18O_H2O)
% GOP = calcGOP_P11(-0.571224,-1.27791,35,15,1.4,-12.1885,-23.8867) 
%   > GOP = 0.2095
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% d17O_O2      = delta17O of O2 in sample (per mil vs atmospheric O2)
% d18O_O2      = delta18O of O2 in sample (per mil vs atmospheric O2)
% SP        = practical salinity of sample (PSS)
% pt        = potential temperature of sample (degrees C)
% kO2       = gas transfer coefficient for O2 (m d^-1) 
% varargin  = optional arguments
%         d17O_H2O = delta17O of water (per mil vs atmospheric O2)
%         d18O_H2O = delta18O of water (per mil vs atmospheric O2)
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% GOP  = gross oxygen production rate in mol O2 m^2 d^-1
%
% -------------------------------------------------------------------------
% OTHER REQUIRED FUNCTIONS:
% ------------------------------------------------------------------------- 
% Use of this function requires installation of the gas_toolbox as well as
% the GSW Toolbox.
%
% CC Manning and DP Nicholson (2017). gas_toolbox.
% https://github.com/dnicholson/gas_toolbox
%
% T McDougall and P Barker (2011) Getting started with TEOS-10 and the
% Gibbs Seawater (GSW) Oceanographic Toolbox. SCOR/IAPSO WG127. Available
% at http://www.teos-10.org.
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% Calculates non-steady state gross oxygen production (GOP) in the mixed
% layer from measurements of the triple oxygen isotopic composition of
% dissolved O2, as done in: Manning et al. (2017) Impact of recently
% upwelled water on productivity quantified by in situ and incubation-based
% methods in Monterey Bay, J. Geophys. Res. Oceans doi:
% 10.1002/2016JC012306 
% and 
% Manning et al. (2017) Revising estimates of aquatic gross oxygen
% production by the triple oxygen isotope method to incorporate the local
% isotopic composition of water (submitted manuscript)
%
% This function uses equation 7 from Prokopenko et al. (2011) Geophys. Res.
% Lett.
% 
% This function includes equilibrium but NOT kinetic isotopic fractionation
% from air-sea gas exchange, nor fractionation associated with
% bubble-mediated exchange. The kinetic and bubble-mediated fractionations
% of 17O with respect to 16O currently have much larger uncertainties than
% the equilibrium fractionation. It also does not include any correction
% for entrainment/mixing of water masses containing oxygen with different
% isotopic compositons between the mixed layer and deeper waters.
%
% In this function, 
% R         = molar ratio
% X         = molar ratio divided by the molar ratio in atmospheric O2
% a         = alpha (fractionation factor)
% lambda    = (a17 - 1)/(a18 - 1), the ratio of two fractionation factors
% d         = lowercase delta in per mil
% D         = uppercase Delta in per meg (1 per meg = 0.001 per mil)
%
% e.g.,
% R18      = r(18O/16O) = the molar ratio of 18O to 16O 
% X18      = R18/R18air 
% a18      = R18a/R18b where a and b are two different substances
%            that can exchange isotopes (e.g., reactant and product)
% lambda   = (a17 - 1)/(a18 - 1) = ratio of fractionation factors
% d18O     = (R18samp/R18std - 1) * 1000 
% D17      = (log(d17O_O2/1000+1) - lambda*log(d18O/1000+1))*10^6 where log is
%           the natural logarithm
%
% 
% STEADY STATE: At steady state as defined in Prokopenko et al. (2011),
% there is no change in D17 with time. This occurs when changes in D17 due
% to photosynthesis and air-sea gas exchange exactly balance each other.
%
% Prokopenko et al. (2011) equation 7 is written in terms of molar ratios,
% R18 = r(18O/16O). However, each ratio in equation 7 can be replaced by X,
% the ratio normalized to air, because R18air and R17air cancel out of the
% equation as written below (Hamme et al. (2012)).
%
% Steady state GOP is calculated as
%
%                      (X17 - X17air*a17eq)/X17 - lambda*(X18 - X18air*a18eq)/X18
% GOP = kO2 * [O2]eq * ------------------------------------------------------------
%                      (X17w*a17p - X17)/X17 - lambda*(X18w*a18p - X18)/X18
%
% where X17 is the ratio of 17O/16O in O2 dissolved in water divided by the
% ratio in air, X17w is the ratio of of 17O/16O in H2O divided by the ratio
% in air, X17air is the ratio of 17O/16O in air divided by the ratio in air
% (and therefore equals 1), a17eq is the equilibrium fractionation factor
% for air-sea gas exchange, and a17p is the fractionation factor for
% photosynthesis, and lambda is the ratio of 17O and 18O fractionation
% factors for respiratory O2 consumption.
%
% NOTES (see function below):
%
% 1) This function uses the isotopic composition of VSMOW with respect to
% atmospheric O2 as defined in Barkan and Luz (2011) Rapid Comm. Mass
% Spectrom., Table 1 average 2007-2009.
%
% 2) Uses lambda = 0.5179(0.0006) where lambda is the ratio of kinetic
% respiratory fractionation factors for 17O and 18O, i.e., when respiration
% is the only process affecting oxygen isotopic composition
%       lambda = 17epsilon/18epsilon = (a17 -1)/(a18 - 1)
% Here, lambda = 0.5179 is from Luz and Barkan (2005) Geochim. Cosmochim.
% Acta, Table 1 weighted average. 
%
% 3) The oxygen isotopic composition of VSMOW is given in an IUPAC
% technical report as 
%       r(18O/16O) = 0.00200520(45) 
%       r(17O/16O) = 0.0003799(8).
% See De Laeter et al. (2003), Pure Appl. Chem. 75(6), page 739
%
% 4) Photosynthetic fractionation factors are calculated from Luz and
% Barkan (2011) Geophys. Res. Lett., Table 1, bottom row. They estimate
% that the isotopic composition of O2 produced by photosynthesis from VSMOW
% is d18O_O2_p = -20.014 per mil and d17O_O2_p = -10.126 per mil with
% respect to amospheric O2. These delta values are for "average
% phytoplankton", and are the average of four different isotopic
% compositions of O2, measured independently for cyanobacteria, green
% algae, diatoms, and coccolithophores.
%   We calculate a18p, the photosynthetic fractionation factor, as
%       a18p = r(18O/16O)p / r(18O/16O)VSMOW
%            = (d18O_O2_p+1)*r(18O/16O)air / r(18O/16O)VSMOW.
%
% 5) Equilibrium isotopic fractionation during air-sea gas exchange is
% calculated as
%      a18eq = r(18O/16O)eq / r(18O/16O)air 
% where r(18O/16O)eq is the ratio of dissolved O2 in water at equilibrium
% with the atmosphere. a18eq is calculated using results of Benson and
% Krause (1980) Limnol. Oceanogr. and Benson and Krause (1984) Limnol.
% Oceanogr., BK80 and BK84 respectively. BK80 find that for freshwater
% between 0 to 60 deg Celsius, the equation
%      delta = -0.730 + (427/Tk) 
% with Tk the temperature in Kelvin, quantifies delta, the per mil
% difference between the abundance of 18O-O2 (m/z=34) and O2 (m/z=32) in
% the dissolved gas in water relative to the abundance in air (BK80 eqn
% 29). BK84 conclude that "the effect of salinity on the oxygen
% fractionation should be undetectable in the range 0 < S < 40" and
% therefore recommend using the freshwater equation from BK80. (See
% Isotopic fractionation of atmospheric oxygen during solution, page 631 of
% BK84.)
%
% a17eq is from Stanley et al. (2010), Appendix A who finds that 17Delta =
% 7.7(3.0) per meg for water at room temperature equilibrated with air,
% using lambda = 0.518. This result is similar to Reuer et al. (2007) who
% found 17Delta = 8 per meg at 11 degC and 25 degC using lambda = 0.516,
% which is approximately equal to 17Delta = 6.5 per meg with lambda =
% 0.518. Using the definition of a18eq from BK84 and noting that
%      (delta17Oeq + 1) = a17eq 
% we can solve the equation for 17D to yield 
%      a17eq = exp(8E-6+0.518*log(a18eq))
%
% -------------------------------------------------------------------------
% REFERENCES:
% -------------------------------------------------------------------------
%
% Barkan, E. and B. Luz (2011). The relationships among the three stable
% isotopes of oxygen in air, seawater and marine photosynthesis. Rapid
% communications in mass spectrometry: RCM, 25(16), 2367-2369.
%
% Benson, B. B. and D. Krause (1980). The concentration and isotopic
% fractionation of gases dissolved in freshwater in equilibrium with the
% atmosphere. 1. Oxygen. Limnology and Oceanography, 25(4), 662-671.
%
% Benson, B. B. and D. Krause (1984). The concentration and isotopic
% fractionation of oxygen dissolved in freshwater and seawater in
% equilibrium with the atmosphere. Limnology and oceanography, 29(3),
% 620-632.
%
% De Laeter, J. R., J. K. Bohlke, P. De Bievre et al. (2003). Atomic
% weights of the elements: Review 2000. Pure Appl. Chem., 75(6), 683-800.
%
% Hamme, R. C., N. Cassar, V.P. Lance et al. (2012), Dissolved O2/Ar and
% other methods reveal rapid changes in productivity during a Lagrangian
% experiment in the Southern Ocean, J. Geophys. Res., 117, C00F12,
% doi:10.1029/2011JC007046.
%
% Luz, B., & Barkan, E. (2005). The isotopic ratios 17O/16O and 18O/16O
% in molecular oxygen and their significance in biogeochemistry. Geochimica
% et Cosmochimica Acta, 69(5), 1099-1110.
%
% Luz, B., and E. Barkan (2011), Proper estimation of marine gross O2
% production with 17O/16O and 18O/16O ratios of dissolved O2, Geophys. Res.
% Lett., 38, L19606, doi:10.1029/2011GL049138.
%
% Prokopenko, M. G., O. M. Pauluis, J. Granger, and L. Y. Yeung (2011),
% Exact evaluation of gross photosynthetic production from the oxygen
% triple-isotope composition of O2: Implications for the net-to-gross
% primary production ratios, Geophys. Res. Lett., 38, L14603,
% doi:10.1029/2011GL047652.
%
% Reuer, M. K., B.A. Barnett, M.L. Bender, P. G. Falkowski and M. B.
% Hendricks (2007). New estimates of Southern Ocean biological production
% rates from O2/Ar ratios and the triple isotope composition of O2. Deep
% Sea Research I, 54(6), 951-974.
%
% Stanley, R. H. R., J. B. Kirkpatrick, N. Cassar, B. A. Barnett, and M. L.
% Bender (2010), Net community production and gross primary production
% rates in the western equatorial Pacific, Global Biogeochem. Cycles, 24,
% GB4001, doi:10.1029/2009GB003651.
%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Define constants and fractionation factors
% -------------------------------------------------------------------------

% delta18O and delta17O of VSMOW in per mil vs air O2
% --- See description note 1 in header
d18O_VSMOW_air = -23.324; 
d17O_VSMOW_air = -11.883;

% Check whether the isotopic compositon of water (d17O_H2O and d18O_H2O) was
% provided as an input. If not provided, set d17O_H2O and d18O_H2O to the
% isotopic composition of VSMOW
if nargin > 5
    d17O_H2O = varargin{1};
    d18O_H2O = varargin{2};
else
    % set to VSMOW
    d17O_H2O = d17O_VSMOW_air; 
    d18O_H2O = d18O_VSMOW_air; 
end

% convert deltas to X, the isotope ratio divided by the ratio of
% atmospheric O2
X17     = d17O_O2/1000+1; 
X18     = d18O_O2/1000+1; 
X17w    = d17O_H2O/1000+1;
X18w    = d18O_H2O/1000+1;
X17air  = 1;
X18air  = 1;

% Ratio of respiration fractionation factors for 17O and 18O
% --- See description note 2 in header
lambda  = 0.5179; 

% Ratio of 18O/16O and 17O/16O of VSMOW
% --- See description note 3 in header
R18_VSMOW = 0.0020052;
R17_VSMOW = 0.0003799;

% Ratio of 18O/16O and 17O/16O of atmospheric O2
R18_air = R18_VSMOW./(d18O_VSMOW_air/1000+1);
R17_air = R17_VSMOW./(d17O_VSMOW_air/1000+1);


% Isotopic composition of photosynthetic O2 when the substrate is water
% with the isotopic composition of VSMOW
% --- See description note 4 in header
d18O_O2_p = -20.014;
d17O_O2_p = -10.126;

% Photosynthetic fractionation factors 
a18p    = (d18O_O2_p/1000+1).*R18_air./R18_VSMOW;
a17p    = (d17O_O2_p/1000+1).*R17_air./R17_VSMOW;

% Equilibrium isotopic fractionation during air-sea gas exchange
% as a function of temperature in Kelvin
% --- See description note 5 in header
% In Matlab, the function log yields the natural logarithm
lambdaS10 = 0.518; % lambda used in Stanley et al. (2010)
Tk      = pt+273.15;
a18eq   = (-0.730+427./Tk)/1000+1;
a17eq   = exp(7.7E-6+lambdaS10*log(a18eq));

% -------------------------------------------------------------------------
% Calculate GOP using equation 7 in Prokopenko et al. (2011) Geophys. Res.
% Lett. (steady state)
% -------------------------------------------------------------------------

% GOP/(kO2 * [O2]eq) = (N1 + N2) / (D1 + D2); 
% GOPr = (N1 + N2) / (D1 + D2);
% GOP = GOPr * k * [O2]eq
N1 = (X17 - X17air.*a17eq)./X17;
N2 = -1*lambda*(X18-X18air.*a18eq)./X18;

D1 = (X17w.*a17p - X17)./X17;
D2 = -1*lambda*(X18w.*a18p - X18)./X18;

GOPr = (N1+N2)./(D1+D2);

% GOP in mol O2 m^-2 d^-1
GOP = GOPr .* kO2 .* gasmoleq(SP,pt,'O2');

end

