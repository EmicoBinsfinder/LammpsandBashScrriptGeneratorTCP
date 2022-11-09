"""
Author: Egheosa Ogbomo


- Sliding speed
- Temperature
- Normal pressure
- System (TCP or Hexadecane etc)

Saves the scripts with the defined parameters as a text file in the specified directory.

"""

SlidingSpeeds = ['0.0002', '0.0003', '0.0004', '0.0005']
Temperatures = ['700']
Pressures = ['1000', '2000', '3000', '4000', '5000']

for SlidingSpeed in SlidingSpeeds:
    for Temperature in Temperatures:
        for Pressure in Pressures:
            text = """
echo both

units real 
atom_style charge
dimension 3 
boundary p p f

#----------------------initial variables-------------------------

variable Temp_eq equal 300 # Equilibration temperature in K
variable       Press_eq  equal 10*9.86923                  # Equilibration pressure in atm (10 MPa)

variable Temp equal {Temperature} # Temperature in K
variable Press_MPa equal {Pressure} # Applied pressure in MPa
variable       Press     equal ${{Press_MPa}}*9.86923        # Applied pressure in atm

variable       Nequil    equal 400000                      # Number of timesteps (fs) to equilibrate (300 K, 1 atm)
variable       Nheat     equal 2000000                     # Number of timesteps (fs) to heat
variable       Ncomp     equal 4000000                     # Number of timesteps (fs) to compress
#variable       Nshear    equal 2000000                     # Number of timesteps (fs) to shear

variable       Wall_v    equal {SlidingSpeed}                      # Wall velocity (A/fs) - equal to 10 m/s

variable       res       equal 0.5                         # Resolution for density profiles etc (A)

timestep       0.25                                        # in fs

#----------------------initial position of atoms------------------------- 

read_data 48_TCP_Fe.data

#----------------------calculate variables-------------------------

variable        Sz       equal lx*ly
variable        z_top    equal (bound(all,zmax))
variable        z_bot    equal (bound(all,zmin))
variable        volume   equal ${{Sz}}*(${{z_top}}-${{z_bot}})

variable        bot_wall_z equal bound(bot_wall,zmax)
variable        top_wall_z equal bound(top_wall,zmin)
variable        delta_h equal (v_top_wall_z-v_bot_wall_z)               # inner wall-wall
variable        delta_out equal bound(all,zmax)-bound(all,zmin)         # outer wall-wall

#-----------------------Wall and region definition ----------------------------------- 

group O type 4
group Fe type 1
group C type 3
group H type 2
group P type 5

variable z_bot equal (bound(all,zmin))
variable z_top equal (bound(all,zmax))
variable wall_z equal 7.1
variable z_top_l equal ${{z_top}}-${{wall_z}}
variable z_bot_h equal ${{z_bot}}+${{wall_z}}

region top_wall block INF INF INF INF ${{z_top_l}} ${{z_top}}
region bot_wall block INF INF INF INF ${{z_bot}} ${{z_bot_h}}

group top_wall region top_wall
group bot_wall region bot_wall
group walls union top_wall bot_wall
group molecules subtract all walls

variable fix_z equal 1.2
variable thermo_z equal 2.4
variable z_top_fix equal ${{z_top}}-${{fix_z}}
variable z_top_thermo equal ${{z_top_fix}}-${{thermo_z}}

variable z_bot_fix equal ${{z_bot}}+${{fix_z}}
variable z_bot_thermo equal ${{z_bot_fix}}+${{thermo_z}}

region top_fixed block INF INF INF INF ${{z_top_fix}} ${{z_top}}
region top_thermo block INF INF INF INF ${{z_top_thermo}} ${{z_top_fix}} 

region bot_fixed block INF INF INF INF ${{z_bot}} ${{z_bot_fix}}
region bot_thermo block INF INF INF INF ${{z_bot_fix}} ${{z_bot_thermo}}

group top_fixed region top_fixed
group top_thermo region top_thermo
group top_free subtract top_wall top_thermo top_fixed

group bot_fixed region bot_fixed
group bot_thermo region bot_thermo
group bot_free subtract bot_wall bot_thermo bot_fixed

group total_check union molecules top_fixed top_thermo top_free bot_fixed bot_thermo bot_free

#----------------------Estimate density-------------------------- 

variable        dens_conv   equal 1.0E24/6.02214857E23                              #g/mol/A^3 to g/cm^3 
variable        dens_liq    equal (mass(molecules)/(lx*ly*v_delta_h))*v_dens_conv

#----------------------Defining the interactions-------------------------- 

pair_style reax/c NULL checkqeq yes safezone 8.0 mincap 800 
pair_coeff * * reaxFF_CHOFeP300 Fe H C O P

#------------------------------------------------------------------------- 

neighbor 2.0 bin #2 because the default parameter for skin is 2 in real units 
neigh_modify every 10 delay 0 check no

fix charge all qeq/reax 1 0.0 10.0 1e-6 reax/c 

#----------------------ReaxFF Energies-------------------------

compute reax all pair reax/c
variable eb equal c_reax[1] # bond energy
variable ea equal c_reax[2] # atom energy
variable elp equal c_reax[3] # lone-pair energy
variable emol equal c_reax[4] # molecule energy (0)
variable ev equal c_reax[5] # valence angle energy
variable epen equal c_reax[6] # double-bond valence angle penalty
variable ecoa equal c_reax[7] # valence angle conjugation energy
variable ehb equal c_reax[8] # hydrogen bond energy
variable et equal c_reax[9] # torsion energy
variable eco equal c_reax[10] # conjugation energy
variable ew equal c_reax[11] # van der Waals energy
variable ep equal c_reax[12] # Coulomb energy
variable efi equal c_reax[13] # electric field energy (0)
variable eqeq equal c_reax[14] # charge equilibration energy

######################################################################
#---------------------ReaxFF Minimization----------------------------#
###################################################################### 

fix freeze_top top_fixed setforce 0.0 0.0 0.0
fix freeze_bot bot_fixed setforce 0.0 0.0 0.0 

thermo 10 
thermo_style custom step temp pe ke press pxx pyy pzz v_delta_h v_dens_liq

#fix minconn all reax/c/bonds 100 emin_bondsu.reaxc 

dump ovito_min all custom 10 dump_min.lammpstrj id type x y z 

minimize 1.0e-7 1.0e-7 40000 40000 
min_style cg

undump ovito_min

#unfix  minconn

#reset_timestep 0

##############################################
#-------------Equilibrate--------------------#
##############################################

#-------------Top wall piston-------------------------

unfix           freeze_top

fix             piston_top top_fixed setforce 0.0 0.0 NULL

# ----------- Heat and Thermostat ----------------------------

group           movable subtract all top_fixed bot_fixed
compute         temp_free movable temp
compute         temp_surf walls temp

velocity        movable create ${{Temp_eq}} 45722345

compute         temp_y_top top_thermo temp/partial 0 1 0
fix             lang_top   top_thermo langevin ${{Temp_eq}} ${{Temp_eq}} $(100.0*dt) 699483 zero yes
fix_modify      lang_top   temp temp_y_top
compute         temp_y_bot bot_thermo temp/partial 0 1 0
fix             lang_bot   bot_thermo langevin ${{Temp_eq}} ${{Temp_eq}} $(100.0*dt) 2847563 zero yes
fix_modify      lang_bot   temp temp_y_bot

###########################################################################
#------------Move top wall down--------------------
#
#variable        move_dist  equal 30              # in A
#variable        move_speed equal 1e-4            # 10 m/s
#variable        move_time  equal 4*(${{move_dist}}/${{move_speed}})
#
#group           most subtract all top_fixed
#
#fix             nve_most most nve
#fix             top_down top_fixed move linear NULL NULL -${{move_speed}}  
#
#thermo          1000
#
#run             ${{move_time}}
#
#unfix           top_down
#unfix           nve_most
#
############################################################################

#-------------Compress to equil pressure--------------

fix             nve_all all nve

#1GPa = 0.14393 kcal mol-1 A-3, 1atm = 0.0000145837 kcal mol-1 A-3, kcal mol-1 A-3 x A2 = kcal mol-1 A-1
variable        applied_press_eq equal ${{Press_eq}}                                                    # 10 MPa
variable        applied_force_eq equal (${{applied_press_eq}}*0.0000145837*${{Sz}})                       # in kcal mol-1 A-1

group top_fixed_fe subtract top_fixed O
variable        z_force_top_fe_eq equal (-${{applied_force_eq}})/(count(top_fixed_fe))

fix             comp_top_fe_eq top_fixed_fe aveforce NULL NULL v_z_force_top_fe_eq

#-------------Bond info------------------------------

fix             bonds_equil all reax/c/bonds 4000 bonds_equil.txt
#fix             species_equil all reax/c/species 1 100 40000 species_equil.txt element Fe O C H

dump            ovito_equil all custom 4000 dump_equil.lammpstrj id type x y z

#-------------Normal and shear forces-----------------------

variable        s_bot equal -f_freeze_bot[1]/0.000014393/${{Sz}}
variable        p_bot equal -f_freeze_bot[3]/0.000014393/${{Sz}}

#------------Thermo Equil/Heat/Comp/Shear-----------------

thermo          4000 
thermo_style    custom step pe etotal c_temp_free c_temp_surf press pxx pyy pzz v_delta_h v_dens_liq v_s_bot v_p_bot
thermo_modify   lost ignore flush yes

#-------------Run Equilibration-------------------

run             ${{Nequil}}

unfix           lang_top
unfix           lang_bot

unfix           bonds_equil
#unfix           species_equil
undump          ovito_equil

###################################################
# ----------------Heating-------------------------#
###################################################

#-------------Heat to set temperature--------------

fix             lang_top_2 top_thermo langevin ${{Temp}} ${{Temp}} $(100*dt) 699483 zero yes
fix_modify      lang_top_2 temp temp_y_top

fix             lang_bot_2 bot_thermo langevin ${{Temp}} ${{Temp}} $(100*dt) 284756 zero yes
fix_modify      lang_bot_2 temp temp_y_bot

#-------------Bond info------------------------------

#fix             bonds_heat all reax/c/bonds 4000 bonds_heat.txt
#fix             species_heat all reax/c/species 1 100 40000 species_heat.txt element Fe O C H

#dump            ovito_heat all custom 4000 dump_heat.lammpstrj id type x y z

#-----------Run Heating---------------------------

#variable        Nheat_10 equal ${{Nheat}}/10
#restart         ${{Nheat_10}} Heating.restart

#run             ${{Nheat}}

#unfix           comp_top_fe_eq

#unfix           bonds_heat
#unfix           species_heat
#undump          ovito_heat

###################################################
# ----------Compression and Shear-----------------#
###################################################

#-------------Compress to set pressure--------------

#1GPa = 0.14393 kcal mol-1 A-3, 1atm = 0.0000145837 kcal mol-1 A-3, kcal mol-1 A-3 x A2 = kcal mol-1 A-1
variable        applied_press equal ${{Press}}
variable        applied_force equal (${{applied_press}}*0.0000145837*${{Sz}})                       #in kcal mol-1 A-1 (/2 if applied to both surfaces)
print           "Applied Force = ${{applied_force}} kcal mol^-1 A^-1"                            
print           "Surface Area = ${{Sz}} A^2"

#variable        z_force_top   equal (-${{applied_force}})/(count(top_fixed))
#print           "z_force_top = ${{z_force_top}} kcal mol-1 A-1"
#fix             comp_top top_fixed aveforce NULL NULL ${{z_force_top}}

group top_fixed_fe subtract top_fixed O
variable        z_force_top_fe equal (-${{applied_force}})/(count(top_fixed_fe))
print           "z_force_top_fe = ${{z_force_top_fe}} kcal mol-1 A-1"

fix             comp_top_fe top_fixed_fe aveforce NULL NULL v_z_force_top_fe

#---------Add Velocity-----------------------------

variable        vel_top equal  ${{Wall_v}}/2     
variable        vel_bot equal -${{Wall_v}}/2 
velocity        top_fixed set ${{vel_top}} 0.0 0.0 units box
velocity        top_fixed set ${{vel_bot}} 0.0 0.0 units box

#---------Calculate Friction-----------------------------

variable        div_s    equal 100
variable        ts_div_s equal ${{Ncomp}}/${{div_s}}

fix             fc_ave all ave/time 1 ${{ts_div_s}} ${{ts_div_s}} v_s_bot v_p_bot file fc_ave.dump

#-------------Bond info-------------------------------

fix             bonds_comp all reax/c/bonds 4000 bonds_comp.txt
#fix             species_comp all reax/c/species 1 100 40000 species_comp.txt element Fe O P C H

dump            ovito_comp all custom 4000 dump_comp.lammpstrj id type x y z

#-----------Run Compression and Shear-----------------------

variable        Ncomp_10 equal ${{Ncomp}}/10
restart         ${{Ncomp_10}} comp.restart

run             ${{Ncomp}}

unfix           bonds_comp
#unfix           species_comp
undump          ovito_comp

""".format(SlidingSpeed=SlidingSpeed, Temperature=Temperature, Pressure=Pressure)
            file = open(
                f"D:/PhD/TCPDecompositionExperiments/Templates/DifferentSlidingSpeeds/{SlidingSpeed[-1]}0ms/{Temperature}K/{Pressure[0]}GPa/TCP.lammps",
                "w")
            file.write(text)
            file.close()