MIDAS
=====

Python for Oceanographic and Atmospheric Data Analysis . Modular Isosurface Data Analysis Software is 
designed to operate using the MOM6 API (Hallberg, Adcroft, Griffies, in preparation). The Pythonic
MIDAS interface rests at the level of the Python and interfaces with C/Fortran95 APIs to MOM6 using the
FMS infrastructure.  The hierarchy of (class or subclass) objects is described here: 

<- API (Pythonic User Interface)->
                .state  (global ocean) [(t),(s),y,x]              <--- Top level class object for MIDAS
                    .supergrid[y,x]                                 <---  Top level grid (FMS supergrid)
                        .grid[y,x]                                    <---  Derived from supergrid (or externally)
                          .metrics[(t),y,x]                             <---  (Optional) metrics for cell lengths, areas and orientation)
                          .domain2d[y,x]                                <---  Property of the grid (FMS or ESMF based object)
                          .bathymetry[y,x]                              <---  Static information about the location of bedrock (relative to a geopotential)
                            .thin_walls[(s),y,x]                          <--- MOM6 specific
                        .potential_temperature(stagger='11,units='degC',...)[(t),s,y,x]
                        .interfaces[(t),s,y,x]                        <---  Interfaces reside on a static grid but vary in time in the vertical
                        .salinity                                     
                          .interfaces                                   <--- Use pointers to interfaces
                          .transport[(t),(z),y,x]                       <--- Use keys to associate transports with tracers
                        .mass transport[(t),(z),y,x]                    
                          .porous_barriers[(t),(s),y,x]                 <--- Volume or mass transports through layers 
                        .[tracers]                                    <--- Additional tracers    
                          .[bling,cobalt,topaz]
                        .surface[(t),y,x]                             <--- Surface information for coupler
                          .state
                            .SSH,SST,SSS,SSU,SSV[(t),y,x]
                            .[tracers][(t),y,x]
                              .[bling,cobalt,topaz]
                            .wave_model[(t),y,x]
                          .fluxes[(t),y,x]
                        .sea_ice[(t),(s),y,x]
                        .icebergs[(t),(s),y,x]
                        .land_ice[(t),(s),y,x]
                        
                        
The current version of MIDAS (as of 11/29/2012) has the following structure fully or partially implemented:

<- MIDAS Snapshot 11/29/2012 ->
                .state  (global ocean) [(t),(s),y,x]              <--- Top level class object for MIDAS
                    .supergrid[y,x]                                 <---  Top level grid (FMS supergrid)
                        .grid[y,x]                                    <---  Derived from supergrid (or externally)
                                                                      <--- Grids are MOM6 specific but generic grids are incorporated as a sub-class
                          .metrics[(t),y,x]                             <---  (Optional) metrics for cell lengths, areas and orientation)
                          .bathymetry[y,x]                              <---  Static information about the location of bedrock (relative to a geopotential)
                        .interfaces[(t),s,y,x]                        <---  Interfaces reside on a static grid but vary in time in the vertical                  
                        .[tracers]                                    <--- Generic tracers    

In summary, the current definition of the ocean "state" is limited and we only have generic variables.  However,
experience indicates that it is easy to extend class attributes. 
    
An object can be initialized using a set of pre-defined methods triggered by input parameters. These parameters are defined
by the INPUT(and OVERRIDE) files or through Python arguments, e.g.:

sgrid=supergrid(path=None,lon=None,lat=None,rotate_sp=False,tripolar_n=False)
grid=sgrid.mom6_grid

A MOM6 state is instantiated at runtime using F95 APIs. In Python, we will call the same APIs but through a different
interface

S=MOM6_state(path=None,grid=grid,fields=[None],...)

A newer option could be:

S=MOM6_state(npes=1,task_names=['ocean','sea_ice','icebergs','land_ice'],pe_list=[None],...)



