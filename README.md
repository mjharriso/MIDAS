MIDAS
=====

Modular Isosurface Data Analysis Software is designed for working 
across multiple environments:

Offline mode:   Read a restart or history file and use to initialize a state.
                State is limited to what is provided (e.g. tracers without
                dynamics). Any number of generic tracers lying on a simple 
                grid with chosen stagger which may or may not include metrics.
                
                Example:
                
                from midas import *
                sgrid=supergrid(MOM6_override)
                grid=sgrid.mom6_rectgrid(refine=2)
                State=generic_state(grid=grid,path='/tmp/foo.nc',fields=['temp','salt'],stagger=['11','11'])
                State.add_field(['taux','/tmp/foo_tau.nc',stagger=['21'])
                
                
                
Online mode:    Initialize a GOLD or MOM6 object from "scratch" or from a restart
                or history file.  Uses "GOLD_input" and "GOLD_override". 
                
                State=mom6_state(grid=grid,params="MOM6_override")
                State.time_step_tracers()
                State.time_step_dynamics()
                State.write_diagnostics()
                

                





using the proposed MOM6 API (Hallberg, Adcroft, Griffies, in preparation). The
MIDAS interface will interface via a Pywrapper with C/Fortran95 APIs to MOM6 using the
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
                        
                        
The current version of MIDAS (as of 11/29/2012) has the following, implemented:

<- MIDAS Snapshot 11/29/2012 ->
                .state  (global ocean) [(t),(s),y,x]              <--- Top level class object for MIDAS
                    .supergrid[y,x]                                 <---  Top level grid (FMS supergrid)
                        .grid[y,x]                                    <---  Derived from supergrid (or externally)
                                                                      <--- Grids are MOM6 specific but generic grids are incorporated as a sub-class
                          .metrics[(t),y,x]                             <---  (Optional) metrics for cell lengths, areas and orientation)
                          .bathymetry[y,x]                              <---  Static information about the location of bedrock (relative to a geopotential)
                        .interfaces[(t),s,y,x]                        <---  Interfaces reside on a static grid but vary in time in the vertical                  
                        .[tracers]                                    <--- Generic tracers    

    
An object can be initialized using a set of pre-defined methods triggered by input parameters. These parameters are provided by 
INPUT(and OVERRIDE) files or through Python arguments, e.g.:

sgrid=supergrid(path=None,lon=None,lat=None,rotate_sp=False,tripolar_n=False,cyclic_x=False)
grid=sgrid.mom6_grid

A MOM6 state is instantiated at runtime using F95 APIs. In Python, run interactively or potentially at runtime,
we will reference the same APIs but through a different interface, likely involving a mediator F90 layer in order to
pass contiguous address spaces of native float,double,string or logical arrays. Currently, we can only instantiate a 
state through an external file. The file can be a regular path or OpenDAP or Python-NetCDF aggregation virtual file:

S=MOM6_state(path=None,grid=grid,fields=[None],...)

A newer option could be:

S=MOM6_state(npes=1,task_names=['ocean','sea_ice','icebergs','land_ice'],pe_list=[None],...)




