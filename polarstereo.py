import numpy as np

pi_180 = np.pi/180.0

def polarstereo_fwd(phi=None,lam=None,a=6.378137e6,e=0.08181919,phi_c=70.,lamb_0=0.0):
#POLARSTEREO_FWD transforms lat/lon data to map coordinates for a polar stereographic system
#   [X,Y]=POLARSTEREO_FWD(LAT,LONG,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY) 
#   X and Y are the map coordinates (scalars, vectors, or matrices of equal size).
#   LAT and LON are in decimal degrees with negative numbers (-) for S and W. 
#   EARTHRADIUS is the radius of the earth defined in the projection
#       (default is 6378137.0 m, WGS84)
#   ECCENTRICITY is the earth's misshapenness 
#       (default is 0.08181919)
#   LAT_TRUE is the latitude of true scale in degrees, aka standard parallel 
#       (default is -70). Note that some NSIDC data use 70 and some use 71.
#   LON_POSY is the meridian in degrees along the positive Y axis of the map 
#       (default is 0)
#   

#   The National Snow and Ice Data Center (NSIDC) and Scientific Committee
#   on Antarctic Research (SCAR) use a version of the polar stereographic
#   projection that Matlab does not have. This file does transformations to
#   map coordinates from geographic coordinates to facilitate
#   comparisons with other datasets.
#
#   Equations from: Map Projections - A Working manual - by J.P. Snyder. 1987 
#   http://kartoweb.itc.nl/geometrics/Publications/Map#20Projections#20-#20A#20Working#20manual#20-#20by#20J.P.#20Snyder.pdf
#   See the section on Polar Stereographic, with a south polar aspect and
#   known phi_c not at the pole.
#
#   See also: MINVTRAN, MFWDTRAN.
#
#   Written by Andy Bliss, 9/12/2011
#   Translated to Python by M Harrison, 10/8/2012

#   Changes since version 01:
#       1. Split into two functions and vectorized code.

#############
#some standard info
#############
#WGS84 - radius: 6378137.0 eccentricity: 0.08181919
#   in Matlab: axes2ecc(6378137.0, 6356752.3142)
#Hughes ellipsoid - radius: 6378.273 km eccentricity: 0.081816153
#   Used for SSM/I  http://nsidc.org/data/polar_stereo/ps_grids.html
#International ellipsoid (following Snyder) - radius: 6378388.0 eccentricity: 0.0819919 

#{
#check the code using Snyder's example. Should get x=-1540033.6; y=-560526.4;
#phi=-75; lam=150;
#x,y=polarstereo_fwd(phi,lam,6378388.0,0.0819919,-71,-100);
#print x,y

############
#check with AntDEM
############
#http://nsidc.org/data/docs/daac/nsidc0304_0305_glas_dems.gd.html
# Center Point of Corner Grid Cell
#x	y	Latitude	Longitude
#test=[-2812000.0  2299500.0   -57.3452815 -50.7255753
#    2863500.0   2299500.0   -57.0043684 51.2342036
#    -2812000.0  -2384000.0  -56.8847122 -130.2911169
#    2863500.0   -2384000.0  -56.5495152  129.7789915];
#x,y=polarstereo_fwd(test[:,2],test[:,3],6378137.0,axes2ecc(6378137.0, 6356752.3),-70,0);
#plt.figure,plt.plot(test[:,0],test[:,1],'.'),plt.plot(x,y,'r+')
#print test[:,1] test[:,0]-x,test[:,1] test[:,1]-y
#error is less than half a meter (probably just round-off error).



############
#input checking
############
    if a is None:
        a=6378137.0; #radius of ellipsoid, WGS84

    if e is None:
        e=0.08181919; #eccentricity, WGS84

    if phi_c is None:
        phi_c=-70; #standard parallel, latitude of true scale

    if lamb_0 is None:
        lamb_0=0; #meridian along positive Y axis

#convert to radians
    phi=phi*pi_180
    phi_c=phi_c*pi_180
    lam=lam*pi_180
    lamb_0=lamb_0*pi_180

#if the standard parallel is in S.Hemi., switch signs.
    if phi_c < 0 :
        pm=-1; #plus or minus, north lat. or south
        phi=-phi
        phi_c=-phi_c
        lam=-lam
        lamb_0=-lamb_0
    else:
        pm=1


#this is not commented very well. See Snyder for details.
    t=np.tan(np.pi/4.0-phi/2.0)/((1.0-e*np.sin(phi))/(1.0+e*np.sin(phi)))**(e/2.0)
#    t =np.sqrt((1.0-np.sin(phi))/(1.0+np.sin(phi))*((1.0+e*np.sin(phi))/(1.0-e*np.sin(phi)))**e)
    t_c=np.tan(np.pi/4.0 - phi_c/2.0)/((1.0-e*np.sin(phi_c))/(1+e*np.sin(phi_c)))**(e/2.0)
#    t_c =np.sqrt((1.0-np.sin(phi_c))/(1.0+np.sin(phi_c))*((1.0+e*np.sin(phi_c))/(1.0-e*np.sin(phi_c)))**e)
    m_c=np.cos(phi_c)/np.sqrt(1.0-(e**2.0)*(np.sin(phi_c))**2.0)
    rho=a*m_c*t/t_c #true scale at lat phi_c

    m=np.cos(phi)/np.sqrt(1.0-(e**(2.0))*(np.sin(phi))**2.0)
    x=pm*rho*np.sin(lam-lamb_0)
    y=-pm*rho*np.cos(lam - lamb_0)
    k=rho/(a*m)

    return x,y

def polarstereo_inv(x=None,y=None,a=6378137.0,e=0.08181919,phi_c=-70.,lamb_0=0.0):
#POLARSTEREO_INV transforms map coordinates to lat/lon data for a polar stereographic system
#   [LAT,LON]=POLARSTEREO_INV(X,Y,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY) 
#   LAT and LON are in decimal degrees with negative numbers (-) for S and W. 
#   X and Y are the map coordinates (scalars, vectors, or matrices of equal size).
#   EARTHRADIUS is the radius of the earth defined in the projection
#       (default is 6378137.0 m, WGS84)
#   ECCENTRICITY is the earth's misshapenness 
#       (default is 0.08181919)
#   LAT_TRUE is the latitude of true scale in degrees, aka standard parallel 
#       (default is -70). Note that some NSIDC data use 70 and some use 71.
#   LON_POSY is the meridian in degrees along the positive Y axis of the map 
#       (default is 0)
#   
#   The National Snow and Ice Data Center (NSIDC) and Scientific Committee
#   on Antarctic Research (SCAR) use a version of the polar stereographic
#   projection that Matlab does not have. This file does transformations to
#   geographic coordinates from map coordinates to facilitate
#   comparisons with other datasets.
#
#   Equations from: Map Projections - A Working manual - by J.P. Snyder. 1987 
#   http://kartoweb.itc.nl/geometrics/Publications/Map#20Projections#20-#20A#20Working#20manual#20-#20by#20J.P.#20Snyder.pdf
#   See the section on Polar Stereographic, with a south polar aspect and
#   known phi_c not at the pole.
#
#   See also: MINVTRAN, MFWDTRAN.
#
#   Written by Andy Bliss, 9/12/2011
#   Translated to Python by M Harrison, 10/8/2012

#   Changes since version 01:
#       1. Split into two functions and vectorized code.

#############
#some standard info
#############
#WGS84 - radius: 6378137.0 eccentricity: 0.08181919
#   in Matlab: axes2ecc(6378137.0, 6356752.3142)
#Hughes ellipsoid - radius: 6378.273 km eccentricity: 0.081816153
#   Used for SSM/I  http://nsidc.org/data/polar_stereo/ps_grids.html
#International ellipsoid (following Snyder) - radius: 6378388.0 eccentricity: 0.0819919 

#{
#check the code using Snyder's example. Should be about -75 (S) and 150 (E).
#x=-1540033.6;
#y=-560526.4;
#[phi,lambda]=polarstereo_inv(x,y,6378388.0,0.0819919,-71,-100);
#print phi,lambda

############
#check with AntDEM
############
#http://nsidc.org/data/docs/daac/nsidc0304_0305_glas_dems.gd.html
# Center Point of Corner Grid Cell
#x	y	Latitude	Longitude
#test=[-2812000.0  2299500.0   -57.3452815 -50.7255753
#    2863500.0   2299500.0   -57.0043684 51.2342036
#    -2812000.0  -2384000.0  -56.8847122 -130.2911169
#    2863500.0   -2384000.0  -56.5495152  129.7789915];
#[phi,lamb]=polarstereo_inv(test[:,0],test[:,1],6378137.0,axes2ecc(6378137.0, 6356752.3),-70,0);
#plt.figure;plt.plot(test[:,3],test[:,2],'.');plt.plot(lamb,phi,'r+')
#print test[:,2],test[:,2]-phi,test[:,3],test[:,3]-lamb
#error is less than half a meter (probably just round-off error).

############
#check with Greenland
############
#projected from the WGS 84 Ellipsoid, with 70 N as the latitude of true scale and a rotation of 45.
#test=[-890000.0 -629000.0 79.9641229 -99.7495626 #center point of cell
#    1720000.0 -629000.0 73.2101234 24.9126514
#    -890000.0 -3410000.0 58.2706251 -59.6277136
#    1720000.0 -3410000.0 55.7592932 -18.2336765];
#[phi,lamb]=polarstereo_inv(test(:,1),test(:,2),6378273,0.081816153,70,-45); #slightly off
#[phi2,lamb2]=polarstereo_inv(test(:,1),test(:,2),6378137.0,0.08181919,70,-45); #correct, FINALLY!
#figure,hold on,plot(test(:,4),test(:,3),'.'),plot(lamb,phi,'r+'),plot(lamb2,phi2,'gx')
#[test(:,3) test(:,3)-phi test(:,3)-phi2],[test(:,4) test(:,4)-lamb test(:,4)-lamb2]
#error is less than half a meter (probably just round-off error).
#}

############
#input checking
############
    if a is None:
        a=6378137.0 #radius of ellipsoid, WGS84

    if e is None:
        e=0.08181919 #eccentricity, WGS84

    if phi_c is None:
        phi_c=-70 #standard parallel, latitude of true scale

    if lamb_0 is None:
        lamb_0=0 #meridian along positive Y axis

#convert to radians
    phi_c=phi_c*pi_180
    lamb_0=lamb_0*pi_180

#if the standard parallel is in S.Hemi., switch signs.
    if phi_c < 0:
        pm=-1 #plus or minus, north lat. or south
        phi_c=-phi_c
        lamb_0=-lamb_0
        x=-x
        y=-y
    else:
        pm=1


#this is not commented very well. See Snyder for details.
    t_c=np.tan(np.pi/4.0 - phi_c/2.0)/((1.0-e*np.sin(phi_c))/(1.0+e*np.sin(phi_c)))**(e/2.0);
    m_c=np.cos(phi_c)/np.sqrt(1.0-e**2.0*(np.sin(phi_c))**2.0);
    rho=np.sqrt(x**2.0+y**2.0); 
    t=rho*t_c/(a*m_c)

#iterate to find phi, with a threshold of pi*1e-8
# phi_alt=pi/2 - 2 * atan(t); #guess for phi
# phiold=phi_alt+10; #+10 to make sure it executes the while loop
# while any(abs(phi_alt(:)-phiold(:)) > pi*1e-8)
#     phiold=phi_alt;
#     phi_alt=pi/2 - 2*atan(t.*((1-e*sin(phi_alt))./(1+e*sin(phi_alt))).^(e/2));
#     #add a break in case it doesn't converge?
# end

#find phi with a series instead of iterating.
    chi=np.pi/2.0 - 2.0 * np.arctan(t)
    phi=chi+((e**2.0)/2.0 + 5.0*(e**4.0)/24.0 + (e**6.0)/12.0 + 13.0*e**(8.0)/360.)*np.sin(2.0*chi) + (7.0*(e**4.0)/48.0 + 29.0*(e**6.0)/240.0 + 811.0*(e**8.0)/11520.)*np.sin(4.0*chi)+ (7.0*(e**6.0)/120.0+81.0*(e**8.0)/1120.)*np.sin(6.0*chi) + (4279.0*(e**8.0)/161280.0)*np.sin(8.0*chi)
# phicheck=phi-phi_alt #difference is usually < 1e-10

    lamb=lamb_0 + np.arctan2(y,x);

#correct the signs and phasing
    phi=pm*phi;
# phi_alt=pm*phi_alt;
    lamb=pm*lamb;
    lamb=np.mod(lamb+np.pi,2.0*np.pi)-np.pi; #want longitude in the range -pi to pi

#convert back to degrees
    phi=phi/pi_180
    lamb=lamb/pi_180

    return lamb,phi
