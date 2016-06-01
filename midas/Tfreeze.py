import numpy as np
"""
==============================

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit   
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.

===============================
"""


def calculate_TFreeze_Millero(S,p):
  """
    This subroutine computes the freezing point potential temparature
  (in deg C) from salinity (in psu), and pressure (in Pa) using the expression
  from Millero (1978) (and in appendix A of Gill 1982), but with the of the
  pressure dependence changed from 7.53e-8 to 7.75e-8 to make this an
  expression for potential temperature (not in situ temperature), using a
  value that is correct at the freezing point at 35 PSU and 5e6 Pa (500 dbar).

 Arguments: S - salinity in PSU.
         pres - pressure in Pa.
  
 ***********************************************************************
 
 """

  cS1= -0.0575
  cS3_2 = 1.710523e-3
  cS2 = -2.154996e-4
  dTFr_dp = -7.75e-8

  T_Fr = S*(cS1+ (cS3_2 * np.sqrt(np.max(S,0.0)) + cS2 * S)) + dTFr_dp*p

  return T_Fr
