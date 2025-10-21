  
 --------------------------------------------
 Physics Parameters:
 -------------------
    gravity:   9.8100000000000005     
    density water:   1025.0000000000000     
    density air:   1.1499999999999999     
    ambient pressure:   101300.00000000000     
    earth_radius:   6367500.0000000000     
    coordinate_system:           2
    sea_level:   0.0000000000000000     
  
    coriolis_forcing: F
    theta_0:   0.0000000000000000     
    friction_forcing: T
    manning_coefficient:   2.5000000000000001E-002
    friction_depth:   100.00000000000000     
  
    dry_tolerance:   1.0000000000000000E-003
  
 --------------------------------------------
 Refinement Control Parameters:
 ------------------------------
    wave_tolerance:   1.0000000000000000E-002
    speed_tolerance:   0.0000000000000000     
    Variable dt Refinement Ratios: T
 
  
 --------------------------------------------
 SETDTOPO:
 -------------
    num dtopo files =            1
    fname:/home/aila/clawpack_src/sulawesi_source_filter/geoclaw_output/fine_B0/dtopo.tt3                                                                       
    topo type:           3
  
 --------------------------------------------
 SETTOPO:
 ---------
    mtopofiles =            2
    
    /home/aila/clawpack_src/sulawesi_source_filter/geoclaw_driver/DataFiles/etopo1_-130_-124_38_45_1min.asc                                               
   itopotype =            3
   mx =          391   x = (  -129.99166666666650      ,  -123.49166666653650      )
   my =          421   y = (   38.008333333333503      ,   45.008333333473502      )
   dx, dy (meters/degrees) =    1.6666666667000000E-002   1.6666666667000000E-002
    
    /home/aila/clawpack_src/sulawesi_source_filter/geoclaw_driver/DataFiles/cc-1sec.asc                                                                   
   itopotype =            3
   mx =         3601   x = (  -124.99977000000000      ,  -123.99977280000000      )
   my =         2449   y = (   41.419950000000000      ,   42.099948095999999      )
   dx, dy (meters/degrees) =    2.7777699999999999E-004   2.7777699999999999E-004
  
   Ranking of topography files  finest to coarsest:            2           3           1
  
  
 --------------------------------------------
 SETQINIT:
 -------------
   qinit_type = 0, no perturbation
  
 --------------------------------------------
 Multilayer Parameters:
 ----------------------
    check_richardson: T
    richardson_tolerance:  0.94999999999999996     
    eigen_method:           4
    inundation_method:           2
    dry_tolerance:   1.0000000000000000E-003
