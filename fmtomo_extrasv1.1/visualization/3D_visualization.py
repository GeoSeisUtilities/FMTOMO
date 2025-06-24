#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  8 17:42:27 2025

@author: donato
"""
## Read output of FMTOMO and produce x,y,z,v text files for easier visulization

# Setting input
folder_path = '/home/donato/Scaricati/observed/'
sez_strike = 45 # angle from North
sez_spacing = 20 # km

# Importing modules
import os
import pandas as pd
import numpy as np
import math
from geographiclib.geodesic import Geodesic
import pygmt

# Defining functions
def get_grid(grid3dgi):
    with open(grid3dgi) as o:
        _ = [o.readline() for i in range(16)]
        Z = o.readline()
        LAT = o.readline()
        LON = o.readline()
    # depth range
    Z = Z.split(' ')
    Zlen = Z.index('c:')
    Z1 = float(Z[0])
    for i in range(Zlen - 1):
        if Z[i] != '':
            j = i
    Z2 = float(Z[j])
    # lat range
    LAT = LAT.split(' ')
    LATlen = LAT.index('c:')
    LAT1 = float(LAT[0])
    for i in range(LATlen - 1):
        if LAT[i] != '':
            j = i
    LAT2 = float(LAT[j])
    # lon range
    LON = LON.split(' ')
    LONlen = LON.index('c:')
    LON1 = float(LON[0])
    for i in range(LONlen - 1):
        if LON[i] != '':
            j = i
    LON2 = float(LON[j])
    if LAT1 < LAT2:
        LATs = LAT1
        LATe = LAT2
    else:
        LATs = LAT2
        LATe = LAT1
    if LON1 < LON2:
        LONs = LON1
        LONe = LON2
    else:
        LONs = LON2
        LONe = LON1
    if Z1 > Z2:
        Zs = Z1
        Ze = Z2
    else:
        Zs = Z2
        Ze = Z1
    return LATs, LATe, LONs, LONe, Zs, Ze

def get_velocity(vgrids):
    with open(vgrids) as o:
        _ = o.readline()
        vgrid = o.readline()
        _ = o.readline()
        _ = o.readline()
        l = o.readlines()
    if len(l[0].split(' '))-(l[0].split(' ')).count('') > 1:
        for j,v in enumerate(l):
            v = v.split(' ')
            while '' in v:
                v.remove('')
            l[j] = v[0]
    while '\n' in l:
        l.remove('\n')
    vgrid = vgrid.split(' ')
    while '' in vgrid:
        vgrid.remove('')
    n_rad = int(vgrid[0])
    n_lat = int(vgrid[1])
    n_lon = int(vgrid[2])
    spac_lat = (LATe - LATs) / (n_lat - 1)
    spac_lon = (LONe - LONs) / (n_lon - 1)
    spac_rad = (Ze - Zs) / (n_rad - 1)
    lon = LONs
    lat = LATs
    depth = Ze
    w = 0
    velocity_grid = pd.DataFrame(columns=["Lon(°)", "Lat(°)", "Dep(km)", 
                    "Dep(m)", "Vel(km/s)", "LonSpacing(°)" , "LatSpacing(°)", 
                    "DepSpacing(km)", "DepSpacing(m)", "n_lon", "n_lat", "n_dep"])
    for i in range(n_rad):
        if i == 0:
            depth = Ze
        else:
            depth = depth - spac_rad
        for j in range(n_lat):
            if j == 0:
                lat = LATs
            else:
                lat = lat + spac_lat
            for k in range(n_lon):
                if k == 0:
                    lon = LONs
                else:
                    lon = lon + spac_lon
                velocity_grid.loc[len(velocity_grid)] = [lon, lat, depth,
                                    depth*1000, float(l[w]), spac_lon, spac_lat,
                                    spac_rad, spac_rad*1000, n_lon, n_lat, n_rad]
                w += 1
    return velocity_grid

def get_sources(source):
    with open(source) as o:
        ll = o.readlines()
    source_loc = pd.DataFrame(columns=["Dep(km)", "Dep(m)", "Lat(°)", "Lon(°)",
                              "ErrDep(km)", "ErrDep(m)", "ErrLat(°)", "ErrLon(°)"])
    for l in ll:
        l = l.split(' ')
        if (len(l)-l.count('')) >= 3:
            while '' in l:
                l.remove('')
            if len(l) > 4:
                source_loc.loc[len(source_loc)] = [-float(l[0]), -float(l[0])*1000,
                            float(l[1]), float(l[2]), float(l[3]), float(l[3])*1000,
                            float(l[4]), float(l[5])]
            else:
                source_loc.loc[len(source_loc)] = [-float(l[0]), -float(l[0])*1000,
                            float(l[1]), float(l[2]), 0, 0 ,0, 0]
    return source_loc

def haversine(lat1, lon1, lat2, lon2):
    '''
    Calculate distance between two points.

    Parameters
    ----------
    lat1 : float
        latitude of the first point (bottom left corner).
    lon1 : float
        longitude of the first point (bottom left corner).
    lat2 : float
        latitude of the second point (top right corner).
    lon2 : float
        longitude of the second point (top right corner).

    Returns
    -------
    distance : float
        distance between the two points.
    '''
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    r = 6371
    distance = r * c
    return distance

def move_point(lat, lon, direction_degrees, distance_km):
    """
    Move a point (lat, lon) along a given bearing (degrees from North, clockwise)
    by the specified distance in kilometers, returning (lat2, lon2).
    Uses WGS-84 ellipsoid via GeographicLib for robustness in all hemispheres.

    Parameters
    ----------
    lat : float
        latitude of the point.
    lon : float
        longitude of the point.
    direction_degrees : float
        direction for moving the point in degrees from North (clockwise).
    distance_km : float
        distance to move along the bearing in kilometers.

    Returns
    -------
    lat_new_deg : float
        latitude of the moved point.
    lon_new_deg : float
        longitude of the moved point.
    """
    distance_m = distance_km * 1000.0 # Convert kilometers to meters
    moved = Geodesic.WGS84.Direct(lat, lon, direction_degrees, distance_m)
    lat_new_deg = moved['lat2']
    lon_new_deg = moved['lon2']
    lon_new_deg = (lon_new_deg + 180) % 360 - 180 # Ensure longitude is in [–180, +180]
    return lat_new_deg, lon_new_deg


def line_intersection(p1, p2, q1, q2):
    """
    Returns the intersection point of two lines (p1-p2 and q1-q2) if it exists.
    Each point is a tuple (x, y).
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = q1
    x4, y4 = q2
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if denom == 0:
        return None
    px = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) / denom
    py = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) / denom
    return (px, py)

# Formatting data in new format
LATs, LATe, LONs, LONe, Zs, Ze = get_grid(f'{folder_path}/invert_p/mkmodel/grid3dgi.in')
Pvelocity = get_velocity(f'{folder_path}/invert_p/vgrids.in')
try:
    Pvelocity_ref = get_velocity(f'{folder_path}/invert_p/mkmodel/reference_velocity.txt')
except:
    Pvelocity_ref = get_velocity(f'{folder_path}/invert_p/vgridsref.in')
Svelocity = get_velocity(f'{folder_path}/invert_s/vgrids.in')
try:
    Svelocity_ref = get_velocity(f'{folder_path}/invert_s/mkmodel/reference_velocity.txt')
except:
    Svelocity_ref = get_velocity(f'{folder_path}/invert_s/vgridsref.in')
source = get_sources(f'{folder_path}/invert_p/sources.in')
source_ref = get_sources(f'{folder_path}/invert_p/sourcesref.in')
receivers = get_sources(f'{folder_path}/invert_p/receivers.in')
receivers_ref = get_sources(f'{folder_path}/invert_p/receiversref.in')
try:
    source_reloc = get_sources(f'{folder_path}/serror.dat')
except:
    print(f'{folder_path}/serror.dat not found')

# Exporting results
path = os.getcwd()+'/'
if not os.path.exists(path+'output_files/'):
    os.makedirs(path+'output_files/')
path = path + 'output_files/'
Pvelocity.to_csv(path+'Pvelocity.csv', index=False)
Pvelocity_ref.to_csv(path+'Pvelocity_ref.csv', index=False)
Svelocity.to_csv(path+'Svelocity.csv', index=False)
Svelocity_ref.to_csv(path+'Svelocity_ref.csv', index=False)
receivers.to_csv(path+'Receivers.csv', index=False)
receivers_ref.to_csv(path+'Receivers_ref.csv', index=False)
source.to_csv(path+'Source.csv', index=False)
source_ref.to_csv(path+'Source_ref.csv', index=False)
try:
    source_reloc.to_csv(path+'Source_reloc.csv', index=False)
except:
    print(f'{folder_path}/serror.dat not found')

# Preparing input for GMT scripts
path = os.getcwd()+'/gmt_files/'
if not os.path.exists(path):
    os.makedirs(path)
with open(path+'data_path', 'w') as o:
    o.write(folder_path)
sez_length = haversine(LATs, LONs, LATe, LONe)
lat_mean = LATs + (LATe - LATs)/2
lon_mean = LONs + (LONe - LONs)/2
if sez_strike >= 0 and sez_strike <= 90:
    y_c, x_c = LATe, LONs
    lim1 = LATs
    lim2 = LONe
    conf1 = '>='
    conf2 = '<='
    cond1 = f'<={LATe}'
    cond2 = f'>={LONs}'
elif sez_strike > 90 and sez_strike <= 180:
    y_c, x_c = LATe, LONe
    lim1 = LATs
    lim2 = LONs
    conf1 = '>='
    conf2 = '>='
    cond1 = f'<={LATe}'
    cond2 = f'<={LONe}'
elif sez_strike > 180 and sez_strike <= 270:
    y_c, x_c = LATs, LONe
    lim1 = LATe
    lim2 = LONs
    conf1 = '<='
    conf2 = '>='
    cond1 = f'>={LATs}'
    cond2 = f'<={LONe}'
elif sez_strike > 270 and sez_strike <= 360:
    y_c, x_c = LATs, LONs
    lim1 = LATe
    lim2 = LONe
    conf1 = '<='
    conf2 = '<='
    cond1 = f'>={LATs}'
    cond2 = f'>={LONs}'
y_c, x_c = move_point(y_c, x_c, (sez_strike + 90) % 360, sez_spacing)
y_st, x_st = move_point(y_c, x_c,(sez_strike+270)%360,sez_length/2)
coor_dict = {}
i=0
while eval(f'{y_c}{conf1}{lim1}') and eval(f'{x_c}{conf2}{lim2}') is True:
    i+=1
    y_end, x_end = move_point(y_st, x_st, sez_strike, sez_length)
    if lim1 == LATs and lim2 == LONe:
        x_1, y_1 = line_intersection((LONs, LATs), (LONs, LATe), (x_st, y_st), (x_end, y_end))
        if y_1 < LATs:
            x_1, y_1 = line_intersection((LONs, LATs), (LONe, LATs), (x_1, y_1), (x_end, y_end))
        x_2, y_2 = line_intersection((LONs, LATe), (LONe, LATe), (x_1, y_1), (x_end, y_end))
        if x_2 > LONe:
            x_2, y_2 = line_intersection((LONe, LATs), (LONe, LATe), (x_1, y_1), (x_2, y_2))
    elif lim1 == LATs and lim2 == LONs:
        x_1, y_1 = line_intersection((LONs, LATe), (LONe, LATe), (x_st, y_st), (x_end, y_end))
        if x_1 < LONs:
            x_1, y_1 = line_intersection((LONs, LATs), (LONs, LATe), (x_1, y_1), (x_end, y_end))
        x_2, y_2 = line_intersection((LONe, LATs), (LONe, LATe), (x_1, y_1), (x_end, y_end))
        if y_2 < LATs:
            x_2, y_2 = line_intersection((LONs, LATs), (LONe, LATs), (x_1, y_1), (x_2, y_2))
    elif lim1 == LATe and lim2 == LONs:
        x_1, y_1 = line_intersection((LONe, LATs), (LONe, LATe), (x_st, y_st), (x_end, y_end))
        if y_1 < LATe:
            x_1, y_1 = line_intersection((LONs, LATe), (LONe, LATe), (x_1, y_1), (x_end, y_end))
        x_2, y_2 = line_intersection((LONs, LATs), (LONe, LATs), (x_1, y_1), (x_end, y_end))
        if x_2 < LONs:
            x_2, y_2 = line_intersection((LONs, LATs), (LONs, LATe), (x_1, y_1), (x_2, y_2))
    elif lim1 == LATe and lim2 == LONe:
        x_1, y_1 = line_intersection((LONs, LATs), (LONe, LATs), (x_st, y_st), (x_end, y_end))
        if x_1 > LONe:
            x_1, y_1 = line_intersection((LONe, LATs), (LONe, LATe), (x_1, y_1), (x_end, y_end))
        x_2, y_2 = line_intersection((LONs, LATs), (LONs, LATe), (x_1, y_1), (x_end, y_end))
        if y_2 > LATe:
            x_2, y_2 = line_intersection((LONs, LATe), (LONe, LATe), (x_1, y_1), (x_2, y_2))
    if eval(f'{y_1}{cond1}') and eval(f'{x_1}{cond2}') is True:
        y_c, x_c = y_1 + (y_2-y_1)/2, x_1 + (x_2-x_1)/2
        sez_tomo = haversine(y_1, x_1, y_2, x_2)
        key = f"Section_{i}"
        coor_dict[key] = {
            "y_start": round(y_1, 3),
            "x_start": round(x_1, 3),
            "y_center": round(y_c, 3),
            "x_center": round(x_c, 3),
            "y_end": round(y_2, 3),
            "x_end": round(x_2, 3),
            "sez_length": sez_tomo
        }
        y_st, x_st = move_point(y_st, x_st, (sez_strike + 90) % 360, sez_spacing)
    else:
        i-=1
        y_c, x_c = y_1 + (y_2-y_1)/2, x_1 + (x_2-x_1)/2
        y_st, x_st = move_point(y_st, x_st, (sez_strike + 90) % 360, sez_spacing)

# Plotting sections
region = [LONs-1, LONe+1, LATs-1, LATe+1]
proj_x=abs(((LONe+1)-(LONs-1))*4)
proj_y=abs(((LATe+1)-(LATs-1))*4)
fig = pygmt.Figure()
pygmt.config(MAP_FRAME_TYPE="plain")
fig.basemap(region=region, projection=f"M{proj_x}/{proj_y}c", frame=True)
fig.coast(land="lightgrey", water="lightblue", shorelines="0.5p,black")
fig.plot(data=np.array([[LONs,LATs,LONe,LATe]]), style='r+s', pen="3p,blue")
for sec in coor_dict.keys():
    start_x = coor_dict[sec]["x_start"]
    start_y = coor_dict[sec]["y_start"]
    end_x = coor_dict[sec]["x_end"]
    end_y = coor_dict[sec]["y_end"]
    fig.plot(x=[start_x, end_x], y=[start_y, end_y], pen="1.5p,red")
    fig.plot(x=[start_x, end_x], y=[start_y, end_y], style="s0.3c", fill="white", pen="red")
    angle_degrees = -sez_strike % 360 + 90
    fig.text(x=start_x, y=start_y, text=f"{sec}", angle=angle_degrees,
    font="10p,Helvetica-Bold,white", fill="red",)
fig.basemap(map_scale="jBR+w100k+l'100 km'") # ScaleBar for 100 km
fig.basemap(compass="jTR+w0.7i") # North arrow
fig.show()

with open(path+'Section_coordinates', 'w') as o:
    for i in coor_dict.keys():
        line = coor_dict[i]
        o.write(f'{line["y_start"]} {line["x_start"]} {line["y_end"]} {line["x_end"]} {line["sez_length"]} {i} \n')
xp = np.unique(Pvelocity["Lon(°)"].values)
xp.tofile(path+'LONGY', sep='\n')
yp = np.unique(Pvelocity["Lat(°)"].values)
yp.tofile(path+'LATY', sep='\n')
zp = np.unique(Pvelocity["Dep(km)"].values)
zp.tofile(path+'DEPTH', sep='\n')
with open(path+'Area_boundaries', 'w') as o:
    o.write(f'{LONs}\n{LONe}\n{LATs}\n{LATe}\n')
with open(path+'Spacing', 'w') as o:
    o.write(f'{sez_spacing/2}')
with open(path+'Project', 'w') as o:
    o.write(f'-Jl{lon_mean}/{LATs}/{LATe}/{lat_mean}/{LONe - LONs}')
