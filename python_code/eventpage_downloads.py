# -*- coding: utf-8 -*-
"""
The routines here create the remaining files in the Downloads section of the NEIC eventpages
"""

from obspy.imaging.scripts import mopad
from numpy import arange, linspace, meshgrid, array, c_, deg2rad, sin, cos, zeros, ones, shape
from glob import glob
import pyproj
import os
from okada_wrapper import dc3dwrapper
import matplotlib.pyplot as plt

def write_Okada_displacements(directory=None):
    print('Writing Okada Displacement File...')

    if directory==None:
        directory = os.getcwd()

    ##########################
    #### FAULT INFORMATION ###
    ##########################

    fsp = open(directory+'/fsp_sol_file','r')
    for line in fsp:
        if line.startswith('% Loc'):
            hypo_lat = float(line.split()[5])
            hypo_lon = float(line.split()[8])
    fsp.close()
    sol = open(directory+'/Solucion.txt','r')
    for line in sol:
        if line.startswith('#Fault_segment'):
            sf_length_as = float(line.split()[7].split('km')[0])
            sf_length_ad = float(line.split()[12].split('km')[0])
    sol.close()
    alpha = 2./3
    #########################################
    ### COORDINATES OF OBSERVATION POINTS ###
    #########################################
    x = linspace(-400,400,40) #in km
    y = linspace(-400,400,40) #in km
    g = pyproj.Geod(ellps='WGS84') # Use WGS84 Ellipsoid
    # Convert to grid of lon/lat around epicenter
    _,xLats,_ = pyproj.Geod.fwd(g, hypo_lon*ones(len(x)), hypo_lat*ones(len(x)), 0*ones(len(x)), x*1000)
    xLons,_,_ = pyproj.Geod.fwd(g, hypo_lon*ones(len(y)), hypo_lat*ones(len(y)), 90*ones(len(y)), y*1000)
    # Make Lon 0 to 360 (not -180 to 180) to avoid plotting issues
    for kLon in range(len(xLons)):
        if xLons[kLon] < 0:
            xLons[kLon] = 360 + xLons[kLon]

    gridLon, gridLat = meshgrid(xLons,xLats)

    LAT=[]; LON=[]; DEP=[]; SLIP=[]; RAKE=[]; STRIKE=[]; DIP=[]; T_RUP=[]; T_RIS=[]; T_FAL=[]; MO=[]

    sol = open(directory+'/Solucion.txt')
    for line in sol:
        if '#' in line: #HEADER LINES
            continue
        if len(array(line.split())) < 4: #FAULT BOUNDARY LINES
            continue
        else: #ACTUAL SUBFAULT DETAILS
            lat, lon, dep, slip, rake, strike, dip, t_rup, t_ris, t_fal, mo = line.split()
            # Make rake be -180 to 180 (not 0-360)
            if float(rake) > 180:
                rake = float(rake) - 360
            # Make Lon 0 to 360 (not -180 to 180)
            if float(lon) < 0:
                lon = float(lon) + 360
            LAT.append(float(lat)); LON.append(float(lon)); DEP.append(float(dep)); SLIP.append(float(slip));
            RAKE.append(float(rake)); STRIKE.append(float(strike)); DIP.append(float(dip));
    sol.close()

    ###############################
    ### CALCULATE DISPLACEMENTS ###
    ###############################
    Nfaults = len(LAT)
    ### Initialize northing, easting, and vertical displacements ###
    ux = zeros(len(xLats)*len(xLons))
    uy = zeros(len(xLats)*len(xLons))
    uz = zeros(len(xLats)*len(xLons))

    #for ksub in range(1):
    for ksub in range(Nfaults):
        if ksub % 100 == 0:
            print('...Subfault: '+str(ksub))
        sf_lon = LON[ksub]
        sf_lat = LAT[ksub]
        strike = STRIKE[ksub]
        dip = DIP[ksub]
        slip = SLIP[ksub]
        rake = RAKE[ksub]
        depth = DEP[ksub]
        ss_slip = (slip/100.) * cos(deg2rad(rake)) 
        ds_slip = (slip/100.) * sin(deg2rad(rake))

        ### Make new x,y vectors, grid of distances in km from current subfault to lon/lat grid defined earlier ###
        fwd_az,b_az,distance_m = pyproj.Geod.inv(g,sf_lon*ones(shape(gridLon)),sf_lat*ones(shape(gridLat)),gridLon,gridLat)
        distance_km = distance_m/1000.
        x = distance_km*sin(deg2rad(fwd_az))
        y = distance_km*cos(deg2rad(fwd_az))
        ### Rotation matrices (code assumes N-S fault, must rotate positions by strike angle and then rotate back) ###
        theta = strike - 90
        theta = deg2rad(theta)
        R_fwd = array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])
        R_inv = array([[cos(-theta),-sin(-theta)],[sin(-theta),cos(-theta)]])

        k=0                            
        for kxy in range(len(x)*len(y)):
            #Calculate on rotated position
            xy=R_fwd.dot(array([[x.flatten()[kxy]], [y.flatten()[kxy]]]))
            success, u, grad_u = dc3dwrapper(alpha, [xy[0], xy[1], 0.0],
                             depth, dip, [-sf_length_as/2, sf_length_as/2], [-sf_length_ad/2, sf_length_ad/2],
                             [ss_slip, ds_slip, 0.0])
            # "un-rotate"
            urot=R_inv.dot(array([[u[0]], [u[1]]]))
            u[0]=urot[0]
            u[1]=urot[1]
        
            #populate output vector
            ux[k]=ux[k]+u[0]
            uy[k]=uy[k]+u[1]
            uz[k]=uz[k]+u[2]
            k+=1

    DISPout = open(directory+'/surface_deformation.disp','w')
    DISPout.write('#Longitude, Latitude, Elevation, Easting Displacement (m), Northing Displacement (m), Vertical Displacement (m)\n')
    for ko in range(len(x)):
        DISPout.write('%10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \t %10.4f \n' % (xLons[ko], xLats[ko], 0, ux[ko], uy[ko], uz[ko]))
    ############
    ### PLOT ###
    ###########
    plt.figure(figsize=(8,6))
    horizontal = (ux**2 + uy**2) **0.5
    horizontal_grid = horizontal.reshape((shape(gridLon)))
    horiz = abs(horizontal)
    horiz.sort()
    minmax = horiz[-2]
    plt.scatter(gridLon,gridLat,marker='s',c=horizontal*100, lw=0, s=100, vmin=0.0, vmax=minmax*100, cmap='rainbow')
    cb = plt.colorbar()
    cb.set_label('Horizontal Surface Displacement (cm)')
    plt.xlim([xLons.min(),xLons.max()])
    plt.ylim([xLats.min(),xLats.max()])
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.savefig(directory+'/Horizontal_Surface_Displacement.png', dpi=300)
    plt.close()

    plt.figure(figsize=(8,6))
    vertical_grid = uz.reshape((shape(gridLon)))
    abs_uz = abs(uz)
    abs_uz.sort()
    minmax = abs_uz[-2]
    plt.scatter(gridLon,gridLat, marker ='s', c=uz*100,lw=0, s=100, vmin=-minmax*100, vmax=minmax*100,cmap='coolwarm')
    cb = plt.colorbar()
    cb.set_label('Vertical Surface Displacement (cm)')
    plt.xlim([xLons.min(),xLons.max()])
    plt.ylim([xLats.min(),xLats.max()])
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.savefig(directory+'/Vertical_Surface_Displacement.png', dpi=300)
    plt.close()

def write_CMTSOLUTION_file(pdefile=None,directory=None):
    print('Writing CMTSOLUTION to file...')

    if directory==None:
        directory = os.getcwd()
    if pdefile==None:
        pdefile = glob('../../../info/*_cmt_CMT')[0]
    pdef = open(pdefile,'r')
    lines = pdef.readlines()
    pde = lines[0]
    eve = lines[1]
    pdef.close()

    LAT=[]; LON=[]; DEP=[]; SLIP=[]; RAKE=[]; STRIKE=[]; DIP=[]; T_RUP=[]; T_RIS=[]; T_FAL=[]; MO=[]

    sol = open(directory+'/Solucion.txt')
    for line in sol:
        if '#' in line: #HEADER LINES
            continue
        if len(array(line.split())) < 4: #FAULT BOUNDARY LINES
            continue
        else: #ACTUAL SUBFAULT DETAILS
            lat, lon, dep, slip, rake, strike, dip, t_rup, t_ris, t_fal, mo = line.split()
            # Make rake be -180 to 180 (not 0-360)
            if float(rake) > 180:
                rake = float(rake) - 360
            LAT.append(float(lat)); LON.append(float(lon)); DEP.append(float(dep)); SLIP.append(float(slip)); 
            RAKE.append(float(rake)); STRIKE.append(float(strike)); DIP.append(float(dip));
            T_RUP.append(float(t_rup)); T_RIS.append(float(t_ris)); T_FAL.append(float(t_fal)); MO.append(float(mo))
    sol.close()

    SUBFAULTS = c_[LAT,LON,DEP,SLIP,RAKE,STRIKE,DIP,T_RUP,T_RIS,T_FAL,MO]
    SUBFAULTS = SUBFAULTS[SUBFAULTS[:,7].argsort()]
    LAT = SUBFAULTS[:,0]; LON = SUBFAULTS[:,1]; DEP = SUBFAULTS[:,2]; SLIP = SUBFAULTS[:,3];
    RAKE = SUBFAULTS[:,4]; STRIKE = SUBFAULTS[:,5]; DIP = SUBFAULTS[:,6]; T_RUP = SUBFAULTS[:,7];
    T_RIS = SUBFAULTS[:,8]; T_FAL = SUBFAULTS[:,9]; MO = SUBFAULTS[:,10]

    CMTout = open(directory+'/CMTSOLUTION','w')        
    for ksubf in range(len(SUBFAULTS)):
        MT = mopad.MomentTensor([STRIKE[ksubf],DIP[ksubf],RAKE[ksubf]])
        MT = array(MT.get_M())
        MT = MT*MO[ksubf]
        
        Mrr = MT[0,0]
        Mtt = MT[1,1]
        Mpp = MT[2,2]
        Mrt = MT[0,1]
        Mrp = MT[0,2]
        Mtp = MT[1,2]
        
        CMTout.write(pde)
        CMTout.write(eve)
        CMTout.write('time shift:{:>13} \n'.format('{:.4f}'.format(T_RUP[ksubf])))
        CMTout.write('half duration:{:>10} \n'.format('{:.4f}'.format(T_RIS[ksubf])))
        CMTout.write('latitude:{:>15} \n'.format('{:.4f}'.format(LAT[ksubf])))
        CMTout.write('longitude:{:>14} \n'.format('{:.4f}'.format(LON[ksubf])))
        CMTout.write('depth:{:>18} \n'.format('{:.4f}'.format(DEP[ksubf])))
        CMTout.write('Mrr: {:>19} \n'.format('{:.6e}'.format(Mrr)))
        CMTout.write('Mtt: {:>19} \n'.format('{:.6e}'.format(Mtt)))
        CMTout.write('Mpp: {:>19} \n'.format('{:.6e}'.format(Mpp)))
        CMTout.write('Mrt: {:>19} \n'.format('{:.6e}'.format(Mrt)))
        CMTout.write('Mrp: {:>19} \n'.format('{:.6e}'.format(Mrp)))
        CMTout.write('Mtp: {:>19} \n'.format('{:.6e}'.format(Mtp)))
    CMTout.close()

def write_Coulomb_file(directory=None,eventID=None):
    print('Writing Coulomb.inp to file...')

    if directory==None:
        directory = os.getcwd()
    if eventID==None:
        eventID='EventX'

    ##########################################################################
    ### GRAB RELEVANT SUBFAULT PARAMETER INFO FROM EVENT_MULT.IN, AND .FSP ###
    ##########################################################################
    
    fsp = open(directory+'/fsp_sol_file','r')
    for line in fsp:
        if line.startswith('% Loc'):
            hypo_lat = float(line.split()[5])
            hypo_lon = float(line.split()[8])
            hypo_dep = float(line.split()[11])
        if line.startswith('% Size'):
            Mw = float(line.split()[13])
            M0 = float(line.split()[16])
        if line.startswith('% Invs : Nx '):
            n_sf_strike = int(line.split()[5]) # Number of subfaults along strike
            n_sf_dip = int(line.split()[8])    # Number of subfaults along dip
        if line.startswith('% Invs : Dx '):
            sf_length_as = float(line.split()[5]) # Subfault length along strike
            sf_length_ad = float(line.split()[9]) # Subfault length along dip
    nft = n_sf_strike * n_sf_dip

    ##########################################################################
    ########################### WRITE HEADER INFO ############################
    ##########################################################################

    coulomb = open(directory+'/Coulomb.inp','w')    
    coulomb.write(eventID+'_coulomb.inp ' + 'inverted by Dara Goldberg (USGS/NEIC) with Mo = ' + str('{:.2e}'.format(M0)) + 'N-m, and Mw = ' + str('{:.2f}'.format(Mw))+'\n')
    coulomb.write('See Hayes(2017), The finite, kinematic rupture properties of great-sized earthquakes since 1990, EPSL 468, 94-100\n')
    coulomb.write('#reg1=  0  #reg2=  0  #fixed= ' + str(nft) +'  sym=  1\n')
    coulomb.write(' PR1=       0.250     PR2=       0.250   DEPTH=       10.000\n')
    coulomb.write('  E1=      8.000e+05   E2=      8.000e+05\n')
    coulomb.write('XSYM=       .000     YSYM=       .000\n')
    coulomb.write('FRIC=          0.400\n')
    coulomb.write('S1DR=         19.000 S1DP=         -0.010 S1IN=        100.000 S1GD=          0.000\n')
    coulomb.write('S2DR=         89.990 S2DP=         89.990 S2IN=         30.000 S2GD=          0.000\n')
    coulomb.write('S3DR=        109.000 S3DP=         -0.010 S3IN=          0.000 S3GD=          0.000\n')
    coulomb.write('\n')
    coulomb.write('  #   X-start    Y-start     X-fin      Y-fin   Kode  rake     netslip   dip angle     top      bot\n')
    coulomb.write('xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n')
    
    ##########################################################################
    ########################## SUBFAULT BY SUBFAULT ##########################
    ##########################################################################

    LAT=[]
    LON=[]
    DEP=[]
    SLIP=[]
    RAKE=[]
    STRIKE=[]
    DIP=[]
    T_RUP=[]
    T_RIS=[]
    T_FAL=[]
    MO=[]
    SS_SLIP = []
    DS_SLIP = []
    
    param = open(directory+'/Solucion.txt')
    for line in param:
        if '#' in line: # HEADER LINES
            continue
        elif len(array(line.split())) < 4: # FAULT BOUNDARY LINES
            continue
        else:
            lat, lon, dep, slip, rake, strike, dip, t_rup, t_ris, t_fal, mo = line.split()
            LAT.append(float(lat)); LON.append(float(lon)); DEP.append(float(dep)); SLIP.append(float(slip)); 
            RAKE.append(float(rake)); STRIKE.append(float(strike)); DIP.append(float(dip)); T_RUP.append(float(t_rup)); 
            T_RIS.append(float(t_ris)); T_FAL.append(float(t_fal)); MO.append(float(mo))
            SS_SLIP.append(float(slip)*cos(deg2rad(float(rake))))
            DS_SLIP.append(float(slip)*sin(deg2rad(float(rake))))
    param.close()
    
    ##########################################################################
    ######### CONVERT COORDINATE SUBFAULT CENTERS TO LOCAL CARTESIAN #########
    ##########################################################################

    Nfaults = len(LAT)
    x = zeros(Nfaults)
    y = zeros(Nfaults)
    g = pyproj.Geod(ellps='WGS84') # Use WGS84 Ellipsoid
    for nsubfault in range(Nfaults):
        backaz, az, dist_m = pyproj.Geod.inv(g, LON[nsubfault], LAT[nsubfault], hypo_lon, hypo_lat)
        x[nsubfault] = (dist_m/1000.) * sin(deg2rad(az))
        y[nsubfault] = (dist_m/1000.) * cos(deg2rad(az))
        
    top_mid_x = zeros(Nfaults)
    top_mid_y = zeros(Nfaults)
    xstart = zeros(Nfaults)
    ystart = zeros(Nfaults)
    xfin = zeros(Nfaults)
    yfin = zeros(Nfaults)
    ztop = zeros(Nfaults)
    zbot = zeros(Nfaults)
    
    for ksubfault in range(Nfaults):
        #top_mid_x[ksubfault] = x[ksubfault]+((sf_length_as/2)*cos(deg2rad(DIP[ksubfault])))*sin(deg2rad(STRIKE[ksubfault]-90))
        #top_mid_y[ksubfault] = y[ksubfault]+((sf_length_as/2)*cos(deg2rad(DIP[ksubfault])))*cos(deg2rad(STRIKE[ksubfault]-90))
        #xstart[ksubfault] = top_mid_x[ksubfault]+(sf_length_as/2)*sin(deg2rad(STRIKE[ksubfault]-180))
        #ystart[ksubfault] = top_mid_y[ksubfault]+(sf_length_as/2)*cos(deg2rad(STRIKE[ksubfault]-180))
        #xfin[ksubfault] = top_mid_x[ksubfault]+(sf_length_as/2)*sin(deg2rad(STRIKE[ksubfault]))
        #yfin[ksubfault] = top_mid_y[ksubfault]+(sf_length_as/2)*cos(deg2rad(STRIKE[ksubfault]))
        #ztop[ksubfault] = DEP[ksubfault]-(sf_length_ad/2)*sin(deg2rad(DIP[ksubfault]))
        #zbot[ksubfault] = DEP[ksubfault]+(sf_length_ad/2)*sin(deg2rad(DIP[ksubfault]))
        xstart[ksubfault] = x[ksubfault]
        ystart[ksubfault] = y[ksubfault]
        xfin[ksubfault] = x[ksubfault] + sf_length_as * sin(deg2rad(STRIKE[ksubfault]))
        yfin[ksubfault] = y[ksubfault] + sf_length_as * cos(deg2rad(STRIKE[ksubfault]))
        ztop[ksubfault] = DEP[ksubfault]
        zbot[ksubfault] = DEP[ksubfault] + sf_length_ad * sin(deg2rad(DIP[ksubfault]))

        out='  1 %10.4f %10.4f %10.4f %10.4f 100 %10.4f %10.4f %10.4f %10.4f %10.4f FFM%d\n' % (xstart[ksubfault],ystart[ksubfault],xfin[ksubfault],yfin[ksubfault],RAKE[ksubfault],SLIP[ksubfault]/100.,DIP[ksubfault],ztop[ksubfault],zbot[ksubfault],ksubfault)
        coulomb.write(out)
    
    ##########################################################################
    ########################### WRITE FOOTER INFO ############################
    ##########################################################################

    deg2km = 111.19
    gsize=6
    
    sx = -(gsize/2.)*deg2km*cos(deg2rad(hypo_lat))
    ex = (gsize/2.)*deg2km*cos(deg2rad(hypo_lat))
    dx = (ex-sx)/100.
    sy = -(gsize/2.)*deg2km
    ey = gsize/2.*deg2km
    dy = (ey-sy)/100.
    min_lat = hypo_lat - (gsize/2.)
    max_lat = hypo_lat + (gsize/2.)
    min_lon = hypo_lon - (gsize/2.)
    max_lon = hypo_lon + (gsize/2.)
    dz = 1.0
    eqdep = -hypo_dep
    
    coulomb.write('\n')
    coulomb.write('     Grid Parameters\n')
    coulomb.write('1 ----------------------------  Start-x =' + '{:16.7f}'.format(sx) + '\n')
    coulomb.write('2 ----------------------------  Start-y =' + '{:16.7f}'.format(sy) + '\n')
    coulomb.write('3 --------------------------   Finish-x =' + '{:16.7f}'.format(ex) + '\n')
    coulomb.write('4 --------------------------   Finish-y =' + '{:16.7f}'.format(ey) + '\n')
    coulomb.write('5 -----------------------   x-increment =' + '{:16.7f}'.format(dx) + '\n')
    coulomb.write('6 -----------------------   y-increment =' + '{:16.7f}'.format(dy) + '\n')
    coulomb.write('\n')
    
    coulomb.write('     Size Parameters\n')
    coulomb.write('1 --------------------------  Plot size =' + '{:16.7f}'.format(gsize) + '\n')
    coulomb.write('2 --------------  Shade/Color increment =' + '{:16.7f}'.format(float(1.0)) + '\n')
    coulomb.write('3 ------  Exaggeration for disp.& dist. =' + '{:16.7f}'.format(float(10000)) + '\n')
    coulomb.write('\n')
    
    coulomb.write('     Cross section default\n')
    coulomb.write('1 ----------------------------  Start-x =' + '{:16.7f}'.format(min_lon) + '\n')
    coulomb.write('2 ----------------------------  Start-y =' + '{:16.7f}'.format(min_lat) + '\n')
    coulomb.write('3 ---------------------------  Finish-x =' + '{:16.7f}'.format(max_lon) + '\n')
    coulomb.write('4 ---------------------------  Finish-y =' + '{:16.7f}'.format(max_lat) + '\n')
    coulomb.write('5 ------------------  Distant-increment =' + '{:16.7f}'.format(dx) + '\n')
    coulomb.write('6 ----------------------------- Z-depth =' + '{:16.7f}'.format(eqdep) + '\n')
    coulomb.write('7 ------------------------  Z-increment =' + '{:16.7f}'.format(dz) + '\n')
    coulomb.write('\n')
    
    coulomb.write('     Map info\n')
    coulomb.write('1 ---------------------------  min. lon =' + '{:16.7f}'.format(min_lon) + '\n')
    coulomb.write('2 ---------------------------  max. lon =' + '{:16.7f}'.format(max_lon) + '\n')
    coulomb.write('3 ---------------------------  zero lon =' + '{:16.7f}'.format(hypo_lon) + '\n')
    coulomb.write('4 ---------------------------  min. lat =' + '{:16.7f}'.format(min_lat) + '\n')
    coulomb.write('5 ---------------------------  max. lat =' + '{:16.7f}'.format(max_lat) + '\n')
    coulomb.write('6 ---------------------------  zero lat =' + '{:16.7f}'.format(hypo_lat) + '\n')
    coulomb.write('7 ------------------------  Z-increment =' + '{:16.7f}'.format(dz) + '\n')
    
    coulomb.close()
