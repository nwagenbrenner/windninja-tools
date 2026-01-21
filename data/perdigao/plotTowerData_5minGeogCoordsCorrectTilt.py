#!/usr/bin/env python


import netCDF4
import numpy as np
from matplotlib import pyplot as plt


### sudo apt install python3-scipy
#import scipy
import scipy.stats


import datetime
from netCDF4 import num2date, date2num, date2index


import os   # blegh, need this to check if folder/file/path exists and delete or make folders accordingly
import shutil  # how dumb, need this to delete a NON empty folder, os only works for empty folder, and deleting files, so dumb need a new import for this



# need to set the desired netcdf file to read in
perdigao_ncFile = "../data/2017-05-04/5-min_GeographicCoords-tiltCorrected/isfs_qc_tiltcor_20170504.nc"


# will be using the exact vars of interest, no need to set a list of which ones to decide on
# so will be using "ldiag", "u", "v", "w", "u_u_", "v_v_", and "w_w_", as inputs are ux, uy, uz, k, epsilon, p
# so the plan is to calculate "spd" and "dir" as needed, and to calculate "tke" as needed as well




#### now set the VenkatramanEtAl Jan 2023 paper towers and plot axes lims as the data and plot style to use

plotCasename = "VenkatramanEtAl2023"

# define the list of towers of interest to be used from the netcdf file and from the line sample post process results
towerList = ["rsw06","tse04","rsw03",  "tse06","tse09","tnw07","tse11",  "tnw10","tse13"]
# note that rsw06 is tower 37, tse04 is tower 20, rsw03 is tower 34, these are the towers from the southwest ridge
# note that tse06 is tower 22, tse09 is tower 25, tnw07 is tower 7, tse11 is tower 27, these are the towers from the valley
# note that tnw10 is tower 10, tse13 is tower 29, these are the towers from the northeast ridge

# looks like we are going to try a crazy subplot with all of them together, with specific grouping
towerGroups = [ 3, 4, 2 ]
towerGroupNames = [ "1stRidgeline", "valley", "2ndRidgeline" ]


# define the set of heights to use for each tower, this can vary tower by tower
# make sure this is the same size and in the same order as the tower list
# order per tower should probably be smallest to greatest to make plot ordering easier
towerHeights = [
                [                      "10m",        "20m", "30m", "40m", "60m"                       ],    # tower 37, rsw06
                #[                      "10m",        "20m", "30m", "40m", "60m", "78m", "80m", "100m" ],    # tower 20, tse04
                [                      "10m",        "20m", "30m", "40m", "60m",        "80m", "100m" ],    # tower 20, tse04
                [                      "10m",        "20m", "30m", "40m", "60m"                       ],    # tower 34, rsw03
                
                [                      "10m",        "20m", "30m", "40m", "60m"                       ],    # tower 22, tse06
                [                      "10m",        "20m", "30m", "40m", "60m",        "80m", "100m" ],    # tower 25, tse09
                [  "4m",  "6m",  "8m", "10m", "12m", "20m", "30m", "40m", "60m"                       ],    # tower  7, tnw07
                [                      "10m",        "20m", "30m", "40m", "60m"                       ],    # tower 27, tse11
                
                [                      "10m",        "20m", "30m", "40m", "60m"                       ],    # tower 10, tnw10
                [                      "10m",        "20m", "30m", "40m", "60m",        "80m", "100m" ]     # tower 29, tse13
               ]
#print(towerHeights)   # for debugging


# specify the desired time range to use out of the full list of times
# do so using datetime()
desiredTimeRangeMin = datetime.datetime(2017,5,4,hour=22,minute=0)
desiredTimeRangeMax = datetime.datetime(2017,5,4,hour=22,minute=30)
#print(desiredTimeRangeMin)
#print(desiredTimeRangeMax)



### for each variable, need to pick and set xLim, yLim, xTicks, yTicks
### if value is "" or [] (lims vs ticks), then skip it and use default standard values for a given lim or set of tick marks
### seems easiest to set the values to the proper size as empty, then to fill if needed by uncommenting the filling of them
spd_xLim = []
spd_yLim = []
spd_xTicks = []
spd_yTicks = []

dir_xLim = []
dir_yLim = []
dir_xTicks = []
dir_yTicks = []

tke_xLim = []
tke_yLim = []
tke_xTicks = []
tke_yTicks = []

## notice that no epsilon or p is done, not plotting those ones in these cases. Technically ux, uy, uz are also side things of interest to plot as well


## spd
spd_yLim = [ 0, 125 ]
spd_yTicks = np.arange( 0, 125+1, 25 )
spd_xLim = [ 0, 7 ]
spd_xTicks = np.arange( 0, 7+1, 1 )

## dir
dir_yLim = [ 0, 125 ]
dir_yTicks = np.arange( 0, 125+1, 25 )
dir_xLim = [ 0, 360 ]
dir_xTicks = np.arange( 0, 360+1, 60 )

## tke
tke_yLim = [ 0, 125 ]
tke_yTicks = np.arange( 0, 125+1, 25 )
tke_xLim = [ 0, 5 ]
tke_xTicks = np.arange( 0, 5+1, 1 )





### function to be edited a lot if needed to get the file paths right and to do general size and type checks on many inputs
### expects the script to be run in the directory where plot outputs are desired
def setScriptDirs( towerList, towerHeights,  plotCasename ):
    
    
    if len(np.shape(towerList)) != 1:
        raise RuntimeError("\n\n!!! input towerList is not a 1D array of values!!!\n\n")
    
    #if len(np.shape(towerHeights)) != 2:
    #    raise RuntimeError("\n\n!!! input towerHeights are not a 2D list of values!!!\n\n")
    
    
    nTowers = len(towerList)
    
    if len(towerHeights) != nTowers:
        raise RuntimeError("\n\n!!! input towerHeights 1st dimension does not equal the number of towers in input towerList !!!\n\n len(towerList) = "+str(len(towerList))+", len(towerHeights) = "+str(len(towerHeights))+"\n\n")
    
    
    
    plotOutputDir = "./" + plotCasename
    if os.path.exists(plotOutputDir):
        # folder exists, need to delete it, this will delete the intermediate output folder if it exists as well, plus all files in the folders
        shutil.rmtree(plotOutputDir)
    ## now the folders are deleted if they existed, can just make them as if they never were there, no need for a check just make them
    os.mkdir(plotOutputDir)
    
    
    return ( nTowers,  plotOutputDir )



### putting the formulas for this into one place
### notice that dir by itself can't be used because it is a python variable, uMag and spd are equivalent, ang and dir are equivalent
def calcSpdAndDirFromComponents( ux, uy, uz ):
    
    
    ## now should be able to calculate and set the desired values
    ## https://www.eol.ucar.edu/content/wind-direction-quick-reference
    ## https://en.wikipedia.org/wiki/Turbulence_kinetic_energy
    
    ##spd = np.sqrt( ux*ux + uy*uy + uz*uz )   # This is the form that I'm used to using
                                               # but we decided that we should use this form for full 3D data, like simulation data,
                                               # and use the 2D form for tower data, and line sample plots
    spd = np.sqrt( ux*ux + uy*uy )   ## Use this 2D form for tower data, and line sample plots, use the above 3D form for full 3D data, like simulation data
    
    
    # test atan2 to see the type
    # https://www.eol.ucar.edu/content/wind-direction-quick-reference
    # if result of atan2(1,-1) is 2.36 radians (135 degrees), then software is y,x. If result is -0.79 radians (-45 degrees), then software is x,y.
    # looks like the result was the 1,-1 form, so this means that I should use the y,x style not the x,y style that I'm used to.
    #print(np.arctan2(1,-1))
    
    
    dirVal = 270.0 - (np.arctan2(uy,ux)*180.0/np.pi)  # rad to degrees
    
    ## had a small problem where the values were mostly correct, but off for some values. Turns out that it just needed one last correction
    ## interestingly, only needed to correct for going OVER 360 degrees, not under 0 degrees, for the above calculation
    if dirVal > 360:
        dirVal = dirVal - 360
    
    
    return ( spd, dirVal )


### notice that dir by itself can't be used because it is a python variable, uMag and spd are equivalent, ang and dir are equivalent
def calcHozComponentsFromSpdAndDir( spd, dirVal ):
    
    
    ### NaN values should just calculate to be NaN values
    
    ### these formulas are what I expected it to be, but ended up being wrong
    ##ux = spd*math.cos(math.radians(dirVal))
    ##uy = spd*math.sin(math.radians(dirVal))
    
    ### wind is in REVERSE direction, formula is different, this ended up being right
    ux = -1*spd*math.sin(math.radians(dirVal))
    uy = -1*spd*math.cos(math.radians(dirVal))
    
    
    return ( ux, uy )


### k and tke are equivalent
def calcTkeFromComponents( uu, vv, ww ):
    
    ## input data should already be averaged covariances of velocity fluctuations
    ## https://en.wikipedia.org/wiki/Turbulence_kinetic_energy
    tke = 0.5*( uu + vv + ww )
    
    return ( tke )



def calc_towerDataRangeAndMeans( tower_var, nTowerHeights,  calcType ):
    
    ### prepare the output storage containers to be filled
    tower_var_min = []
    tower_var_max = []
    tower_var_mean = []
    tower_var_minToMeanDist = []
    tower_var_meanToMaxDist = []
    
    for heightIdx in range(nTowerHeights):
        
        current_data = tower_var[heightIdx]#[:]   # [:] shouldn't be needed now, is the time dimension
        ##print(current_data)
        
        if calcType != "angle":
        
            current_min = np.min(current_data)
            current_max = np.max(current_data)
            current_mean = np.mean(current_data)
            
            current_minToMeanDist = current_mean - current_min
            current_meanToMaxDist = current_max - current_mean
        
        else:
            
            ### use the scipy circular mean calculation, rather than manually doing a polar coordinate transformation myself
            current_min = np.min(current_data)
            current_max = np.max(current_data)
            ##print("min = "+str(current_min)+", max = "+str(current_max))
            
            current_mean = scipy.stats.circmean(current_data, high=360, low=0)
            current_minToMeanDist = current_mean - current_min
            current_meanToMaxDist = current_max - current_mean
            ##print("scipy")
            ##print("mean = "+str(current_mean)+", minToMeanDist = "+str(current_minToMeanDist)+", meanToMaxDist = "+str(current_meanToMaxDist))
        
        tower_var_min.append( current_min )
        tower_var_max.append( current_max )
        tower_var_mean.append( current_mean )
        tower_var_minToMeanDist.append( current_minToMeanDist )
        tower_var_meanToMaxDist.append( current_meanToMaxDist )
    
    #print(tower_var_mean)
    #print(np.asarray(tower_var_mean).shape)
    
    return ( tower_var_min, tower_var_max, tower_var_mean, tower_var_minToMeanDist, tower_var_meanToMaxDist )



### function for reading in perdigao netcdf data
def read_perdigaoTowerData( perdigao_ncFile,  towerList, towerHeights,  desiredTimeRangeMin, desiredTimeRangeMax,   plotOutputDir ):
    
    
    nTowers = len(towerList)
    
    nTowerHeights = [ len(towerHeights[x]) for x in range(nTowers)] # a weird one, is a list of nVals for each tower
    #print(towerHeights)
    #print(nTowerHeights)
    
    
    # now read in the netcdf file
    rootgrp = netCDF4.Dataset( perdigao_ncFile, "r" )
    #print(rootgrp)
    #print(rootgrp.data_model)
    #print(rootgrp.ncattrs())
    #print(rootgrp.dimensions)
    #print(rootgrp.variables)
    
    
    ### example reading in and accessing a variable and its dimensions and nvals
    #temp = rootgrp.variables['u_10m_tse04'] # tower 20, as netcdf variable object
    #print(temp)
    #print(temp.dimensions)
    #print(temp.shape)
    #print(temp._FillValue)
    #print(temp[:])
    #temp = rootgrp.variables['u_10m_tse04'][:] # tower 20, as variable values, so the object .dimensions and .shape are no longer accessible
    #print(temp)
    
    
    ### now to read in and process time variables
    #baseTime = rootgrp.variables['base_time'][:]   # hrm, the date time conversions somehow don't end up needing this after all
    rawTimes = rootgrp.variables['time']   # the full time list, as netcdf variable object
    #print(rawTimes)
    #print(rawTimes.units)
    #print(rawTimes[:])
    #print(type(rawTimes[:]))
    #print(type(rawTimes[0]))
    nRawTimes = len(rawTimes[:])
    #print(nRawTimes)
    raw_dt = rawTimes[1] - rawTimes[0]
    #print(raw_dt)
    #print(type(raw_dt))
    #raw_dt = datetime.timedelta(seconds=raw_dt)
    #print(type(raw_dt))
    #print(raw_dt.total_seconds())
    
    #print('units = %s, values = %s' % (rawTimes.units, rawTimes[:]))
    rawTimesAsDates = num2date(rawTimes[:], rawTimes.units)
    #print([date.strftime('%Y-%m-%d %H:%M:%S') for date in rawTimesAsDates[:]])
    #rawDatesAsTimes = date2num(rawTimesAsDates, rawTimes.units)
    #print('units = %s, values = %s' % (rawTimes.units, rawDatesAsTimes[:]))
    #print(type(rawDatesAsTimes[:]))
    #print(type(rawDatesAsTimes[0]))
    
    # now get the indices of the time range from the time list
    # could do a for loop, but date2index() appears to already exist
    # 
    # hrm, it seems to get the times after the first one, and before the last one, 
    # so I guess need to be more clever after all
    # for now I'll just manually do the before and after as I desire, 
    # cause probably going to need to modify it again later anyways
    # looks like the indexing dropped one when doing the slicing for the second index, 
    # is fixed by adjusting at the slice indexing not here, it is working fine here
    # 
    # caught an additional problem, the times are at centers of the averaged times, so like a cell center vs face center problem
    # so they are stated at centers but are FOR a face center range. So the indexing actually needs to be handling the face center range
    # for its choice. So I guess let's make a new set of time variables representing the range the times are for, the face values list
    # then search for the desired times using that list
    # this means subtracting 1/2*dt to the list and adding a time at the final time + 1/2*dt to get the list of time ranges represented by the times list
    # 
    # hrm, when I went to do this, I had a ROUGH TIME. Turns out there are differences between datetime objects, timedelta objects, and netcdf objects
    # and the date2index() function requires netcdf object as input, NOT a datetime object. So while I was able to figure out how to manipulate stuff
    # as datetime objects, was not effective. The only working fix was to make a new diskless in memory netcdf file in memory to create a full new
    # list of times as a netcdf object, trying to make a copy of the existing time netcdf object adding in space gives a trying to write to a read only
    # file error, so a new diskless object is the only easy way while still being able to use date2index().
    # when I made the diskless netcdf dataset to make the new netcdf list of times object for date2index(), I made the mistake of closing the netcdf
    # file before using the netcdf object time list, and it also gave the nasty errors of before, just be aware.
    
    #desiredTimeRangeMinIndex = date2index(desiredTimeRangeMin,rawTimes,select='before')#'before')#'nearest')   # for this case, 'nearest' is returning 22:07:30 which is the same as 'before'
    #desiredTimeRangeMaxIndex = date2index(desiredTimeRangeMax,rawTimes,select='after')#'after')#'nearest')    # for this case, 'nearest' is returning 22:37:30 which is 1 less than 'after'
    ## now get the desired times from the full time list
    #times = rawTimes[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex+1]   # darn python slicing is so messed up, had to add 1 to get the right end spot
    ##print('units = %s, values = %s' % (rawTimes.units, times[:]))
    #timesAsDates = num2date(times[:], rawTimes.units)
    ##print([date.strftime('%Y-%m-%d %H:%M:%S') for date in timesAsDates[:]])
    #nTimes = len(times)
    ##print(nTimes)
    
    nc = netCDF4.Dataset('diskless_example.nc','w',diskless=True,persist=False) #'a' was causing trouble, no need to 'append' when it is supposed to be a diskless dataset
    d = nc.createDimension('time',None)
    v = nc.createVariable('time',np.float64,('time',))
    v.units = rawTimes.units
    timesRange_list = rawTimes[:] - raw_dt/2        # phew, didn't need to convert raw_dt into a timedelta object with: datetime.timedelta(seconds=raw_dt/2)
    timesRange_list = np.append( timesRange_list, rawTimes[-1] + raw_dt/2 )
    #print(timesRange_list)
    v[:] = timesRange_list[:]
    #print(v)
    #print(v.units)
    #print(v[:])
    ##nc_rawTimesRange = v[:]   # either works fine, turns out that just using v directly into date2index() was also doable
    nc_rawTimesRange = nc.variables['time']
    #print(nc_rawTimesRange)
    #print(nc_rawTimesRange.units)
    #print(nc_rawTimesRange[:])
    rawTimesRange = nc_rawTimesRange[:]
    #print(rawTimesRange)
    rawTimesRangeAsDates = num2date(rawTimesRange[:], rawTimes.units)
    #print([date.strftime('%Y-%m-%d %H:%M:%S') for date in rawTimesRangeAsDates[:]])
    desiredTimeRangeMinIndex = date2index(desiredTimeRangeMin,nc_rawTimesRange,select='nearest')#'before')#'nearest')   # for this case, 'nearest' is returning 22:10:00 for the range, and 22:12:30 for the center, which is the desired set of times
    desiredTimeRangeMaxIndex = date2index(desiredTimeRangeMax,nc_rawTimesRange,select='nearest')#'after')#'nearest')    # for this case, 'nearest' is correctly returning 22:40:00 for the range, incorrectly returning 22:42:30 for the center if keeping the +1 on that index and correctly returning 22:37:30  for the center if dropping the +1 on that index
    nc.close()  # used to have this before the date2index() function, and it was causing me an error claiming that the netcdf variable was unusable or wrong somehow
    #print(desiredTimeRangeMinIndex)
    #print(desiredTimeRangeMaxIndex)
    
    # now get the desired times from the full time list
    times = rawTimes[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]   # had to drop the added +1 to get the right end spot
    timesRange = rawTimesRange[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex+1]     # darn python slicing is so messed up, had to add 1 to get the right end spot. Interesting, had to drop it for the center time to get that right.
    #print('units = %s, values = %s' % (rawTimes.units, times[:]))
    #print('units = %s, values = %s' % (rawTimes.units, timesRange[:]))
    timesAsDates = num2date(times[:], rawTimes.units)
    timesRangeAsDates = num2date(timesRange[:], rawTimes.units)
    #print([date.strftime('%Y-%m-%d %H:%M:%S') for date in timesAsDates[:]])
    #print([date.strftime('%Y-%m-%d %H:%M:%S') for date in timesRangeAsDates[:]])
    nTimes = len(times)
    #print(nTimes)
    nTimesRange = len(timesRange)
    #print(nTimesRange)
    
    
    
    ### output the time info to a text file, a bit easier than just having it in the plot title
    ### use this file to debug the selected time vs what time you wanted
    ### https://www.geeksforgeeks.org/reading-writing-text-files-python/#
    outputFilename = plotOutputDir + "/zz_timeInfo.txt"
    outputFile = open( outputFilename, "w" )
    outputFile.write( 'desiredTimeRange: ' + desiredTimeRangeMin.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + desiredTimeRangeMax.strftime('%Y-%m-%d %H:%M:%S') + '\n' )
    outputFile.write( 'actualTimeRange: ' + timesRangeAsDates[0].strftime('%Y-%m-%d %H:%M:%S') + ' to ' + timesRangeAsDates[-1].strftime('%Y-%m-%d %H:%M:%S') + '\n' )
    outputFile.write( 'actualTimes: ' + timesAsDates[0].strftime('%Y-%m-%d %H:%M:%S') + ' to ' + timesAsDates[-1].strftime('%Y-%m-%d %H:%M:%S') + '\n' )
    outputFile.close()
    
    
    
    ## now read in all desired tower data that is independent of the desired times but is still on a per tower basis
    towerLats = -99999.999*np.ones(nTowers)
    towerLons = -99999.999*np.ones(nTowers)
    for towerIdx in range(nTowers):
        lat_data_name = "latitude_" + towerList[towerIdx]
        towerLats[towerIdx] = rootgrp.variables[lat_data_name][:]
        lon_data_name = "longitude_" + towerList[towerIdx]
        towerLons[towerIdx] = rootgrp.variables[lon_data_name][:]
    #print(towerLats)
    #print(towerLons)
    
    
    
    ## now read in all the desired data for the desired times into a data structure
    ## this is a complex data structure, the method for creating/reading it in 
    ## involves starting with an empty 1D list and making it 4D as I go
    ## https://www.adamsmith.haus/python/answers/how-to-create-a-3d-array-in-python
    # overall lists, to be filled with lists, which are lists of lists. eww
    ### WATCH OUT FOR [idx1][idx2] vs [idx1,idx2] list vs vector indexing/slicing!!! If it is a list, use list indexing methods, if it is a vector part in a list, use vector slicing, or you get some weird quirky undefined behavior types of data structures. It is almost easiest to always pull out the given dimension into the next size lower dimension with a single [:].
    tower_ux = []
    tower_uy = []
    tower_uz = []
    tower_tke = []
    tower_spd = []
    tower_dir = []
    for towerIdx in range(nTowers):
        tower_ux.append([])
        tower_uy.append([])
        tower_uz.append([])
        tower_tke.append([])
        tower_spd.append([])
        tower_dir.append([])
        for heightIdx in range(nTowerHeights[towerIdx]):
            #print(heightIdx)
            ##tower_ux[towerIdx].append([])   ## oops, I apparently was throwing in one extra layer to the data, somehow it worked fine despite it though
            ##tower_uy[towerIdx].append([])
            ##tower_uz[towerIdx].append([])
            ##tower_tke[towerIdx].append([])
            ##tower_spd[towerIdx].append([])
            ##tower_dir[towerIdx].append([])
            
            
            towerDataString = "_" + towerHeights[towerIdx][heightIdx] + "_" + towerList[towerIdx]
            current_nc_ux = rootgrp.variables[ "u"+towerDataString ]
            current_nc_uy = rootgrp.variables[ "v"+towerDataString ]
            current_nc_uy = rootgrp.variables[ "v"+towerDataString ]
            current_nc_uz = rootgrp.variables[ "w"+towerDataString ]
            current_nc_uu = rootgrp.variables[ "u_u_"+towerDataString ]
            current_nc_vv = rootgrp.variables[ "v_v_"+towerDataString ]
            current_nc_ww = rootgrp.variables[ "w_w_"+towerDataString ]
            ##current_nc_spd = rootgrp.variables[ "spd"+towerDataString ]  ## uncomment out if want to read in spd and dir instead of calculating them from the components
            ##current_nc_dir = rootgrp.variables[ "dir"+towerDataString ]
            #print(current_nc_ux)
            
            current_allTimes_ux = current_nc_ux[:]
            current_allTimes_uy = current_nc_uy[:]
            current_allTimes_uz = current_nc_uz[:]
            current_allTimes_uu = current_nc_uu[:]
            current_allTimes_vv = current_nc_vv[:]
            current_allTimes_ww = current_nc_ww[:]
            ##current_allTimes_spd = current_nc_spd[:]
            ##current_allTimes_dir = current_nc_dir[:]
            #print(current_allTimes_ux)
            
            current_ux = current_allTimes_ux[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]#+1]   # darn python slicing is so messed up, at one point in time, I had to add 1 to get the right end spot
            current_uy = current_allTimes_uy[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]
            current_uz = current_allTimes_uz[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]
            current_uu = current_allTimes_uu[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]
            current_vv = current_allTimes_vv[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]
            current_ww = current_allTimes_ww[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]
            ##current_spd_raw = current_allTimes_spd[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]
            ##current_dir_raw = current_allTimes_dir[desiredTimeRangeMinIndex:desiredTimeRangeMaxIndex]
            #print(current_ux)
            
            
            ### always assumes output units are "m/s", so no need to set a variable for it
            #print("current_nc_ux.units")
            #print(current_nc_ux.units)
            #print("current_nc_uy.units")
            #print(current_nc_uy.units)
            #print("current_nc_uz.units")
            #print(current_nc_uz.units)
            
            ### always assumes output units are (m/s)^2
            #print("current_nc_uu.units")
            #print(current_nc_uu.units)
            #print("current_nc_vv.units")
            #print(current_nc_vv.units)
            #print("current_nc_ww.units")
            #print(current_nc_ww.units)
            
            
            
            ## now should be able to calculate and set the desired values
            
            current_spd = []
            current_dir = []
            current_tke = []
            for timeIdx in range(len(current_ux)):
                
                ##current_spd_pt = current_spd_raw[timeIdx]
                ##current_dir_pt = current_dir_raw[timeIdx]
                current_spd_pt, current_dir_pt = calcSpdAndDirFromComponents( current_ux[timeIdx], current_uy[timeIdx], current_uz[timeIdx] )  ## comment out if reading in vals directly
                current_spd.append( current_spd_pt )
                current_dir.append( current_dir_pt )
                
                current_tke_pt = calcTkeFromComponents( current_uu[timeIdx], current_vv[timeIdx], current_ww[timeIdx] )
                current_tke.append( current_tke_pt )
            
            
            ## now append the values to the datasets, this is an all desired times at once per time append
            ##tower_ux[towerIdx][heightIdx].append( current_ux )   ## oops, I apparently was throwing in one extra layer to the data, somehow it worked fine despite it though
            ##tower_uy[towerIdx][heightIdx].append( current_uy )
            ##tower_uz[towerIdx][heightIdx].append( current_uz )
            ##tower_tke[towerIdx][heightIdx].append( current_tke )
            ##tower_spd[towerIdx][heightIdx].append( current_spd )
            ##tower_dir[towerIdx][heightIdx].append( current_dir )
            tower_ux[towerIdx].append( current_ux )
            tower_uy[towerIdx].append( current_uy )
            tower_uz[towerIdx].append( current_uz )
            tower_tke[towerIdx].append( current_tke )
            tower_spd[towerIdx].append( current_spd )
            tower_dir[towerIdx].append( current_dir )
            
        
    #print(len(tower_ux))
    #print(len(tower_ux[0]))
    #print(len(tower_ux[0][0]))
    #print(tower_ux)
    #print(tower_uy)
    #print(tower_uz)
    #print(tower_tke)
    #print(tower_spd)
    #print(tower_dir)
    
    # now finished with the netcdf file, close it
    rootgrp.close()
    
    
    return ( nTowers, nTowerHeights,  times, timesRange, timesAsDates, timesRangeAsDates, nTimes, nTimesRange,  towerLats, towerLons,  tower_ux, tower_uy, tower_uz, tower_tke, tower_spd, tower_dir )



def plotTowerVsLineSampleData( var_name, var_units,  axesLimType, var_xLim, var_yLim, var_xTicks, var_yTicks,  towerGroups, towerGroupNames,  var_towerData,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename ):
    
    
    nTowerGroups = len(towerGroups)
    correcttowerGroups = False
    ## correct for if no tower groups
    ## easiest to just construct it as a single tower group, using nTowers as the tower group
    if nTowerGroups == 0:
        towerGroups = [ nTowers ]
        nTowerGroups = 1
        correcttowerGroups = True
    
    
    
    max_nTowersPerGroup = np.max(towerGroups)
    plt.rc('figure', figsize=(18, 10))  # set the figure size (in inches)
    fig, axs = plt.subplots(nTowerGroups,max_nTowersPerGroup, constrained_layout=True)
    #print(axs)
    
    #suptitleString = "perdigao data for " + var_name + " (" + var_units + ")"
    #fig.suptitle(suptitleString, fontsize=20)
    
    # set the towerIdx counter for the new loop
    towerIdx = 0
    for towerGroupListIdx in range(nTowerGroups):
        
        current_nTowers = towerGroups[towerGroupListIdx]
        current_towerGroupName = towerGroupNames[towerGroupListIdx]
        
        ## keep using towerIdx for references to full tower data, replace towerIdx with towerGroupIdx for all references to plot towers, 
        ## so no change for almost anything except for axes[towerGroupListIdx, towerIdx] becomes axes[towerGroupListIdx, towerGroupIdx]
        #for towerIdx in range(nTowers):
        for towerGroupIdx in range(current_nTowers):
            
            #print('towerGroupListIdx = ' + str(towerGroupListIdx) + ', current_nTowers = ' + str(current_nTowers) + ', towerGroupIdx = ' + str(towerGroupIdx) + ', towerIdx = ' + str(towerIdx))
            
            ### plot the tower data stuff first
            
            ## cut off the "m" before the conversion to float, needed to properly put distance between the values as the axis labels
            ## no need to correct this for z_ground (z[0]), these are already in terms of AGL
            current_z = [ float(x[0:-1]) for x in towerHeights[towerIdx] ]
            #print(current_z)
            #print(np.asarray(current_z).shape)
            
            current_towerData = var_towerData[towerIdx]
            current_nTowerHeights = nTowerHeights[towerIdx]
            
            if var_name == "dir":
                
                calcType = "angle"
                current_var_min, current_var_max, current_var_mean, current_var_minToMeanDist, current_var_meanToMaxDist = calc_towerDataRangeAndMeans( current_towerData, current_nTowerHeights,  calcType )
                
            else:
                
                calcType = "standard"
                current_var_min, current_var_max, current_var_mean, current_var_minToMeanDist, current_var_meanToMaxDist = calc_towerDataRangeAndMeans( current_towerData, current_nTowerHeights,  calcType )
            
            for zIdx in range(len(current_z)):
                current_towerHeightData = current_towerData[zIdx]
                current_z_val = current_z[zIdx]
                current_z_vals = np.repeat(current_z_val,len(current_towerHeightData))
                if zIdx == 0:
                    axs[towerGroupListIdx, towerGroupIdx].plot( current_towerHeightData, current_z_vals, linestyle='None', marker='o', markersize=7.5, color="grey", label="Observed 5-min mean" )
                else:
                    axs[towerGroupListIdx, towerGroupIdx].plot( current_towerHeightData, current_z_vals, linestyle='None', marker='o', markersize=7.5, color="grey" )
            
            axs[towerGroupListIdx, towerGroupIdx].plot( current_var_mean, current_z, 'k^', markersize=7.5, label="Observed 30-min mean" )
            
            ##axs[towerGroupListIdx, towerGroupIdx].plot( current_var_min, current_z, 'k|', markersize=12.5, markeredgewidth=2.5)
            ##axs[towerGroupListIdx, towerGroupIdx].plot( current_var_max, current_z, 'k|', markersize=12.5, markeredgewidth=2.5)
            
            
            ### now set the axis stuff that is always the same
            
            axes_title = str(towerList[towerIdx])
            axs[towerGroupListIdx, towerGroupIdx].set_title(axes_title, fontsize=20)
            axs[towerGroupListIdx, towerGroupIdx].tick_params(axis='x', labelsize=18)  ## maybe make size 20
            axs[towerGroupListIdx, towerGroupIdx].tick_params(axis='y', labelsize=18)  ## maybe make size 20
            axs[towerGroupListIdx, towerGroupIdx].tick_params(axis="x", length=4.25, width=1.3)
            axs[towerGroupListIdx, towerGroupIdx].tick_params(axis="y", length=4.25, width=1.3)
            axs[towerGroupListIdx, towerGroupIdx].grid(True)
            axs[towerGroupListIdx, towerGroupIdx].minorticks_on()
            axs[towerGroupListIdx, towerGroupIdx].tick_params(which='minor', length=3, width=0.75)
            if ( towerGroupListIdx == nTowerGroups-1 ):
                ### trying to get rid of italics while still going latex https://stackoverflow.com/questions/19671659/remove-italics-in-latex-subscript-in-matplotlib
                ### challenge is wanting to still keep the original font? cause superscript still retains new font
                if var_name == "spd":
                    xLabelString = "wind speed ($\mathrm{m}$ $\mathrm{s^{-1}}$)"
                elif var_name == "dir":
                    #xLabelString = "direction ($\mathrm{°}$)"
                    xLabelString = "direction ($\mathrm{^{o}}$)"
                elif var_name == "tke":
                    xLabelString = "TKE ($\mathrm{m^{2}}$ $\mathrm{s^{-2}}$)"
                else:
                    xLabelString = var_name + " (" + var_units + ")"
                axs[towerGroupListIdx, towerGroupIdx].set_xlabel(xLabelString, fontsize=20)
            if ( towerGroupIdx == 0 ):
                axs[towerGroupListIdx, towerGroupIdx].set_ylabel("z (m AGL)", fontsize=20)
            
            
            ### now set the axis stuff that is data specific
            
            if len(var_xLim) != 0:
                axs[towerGroupListIdx, towerGroupIdx].set_xlim(var_xLim)
            if len(var_yLim) != 0:
                axs[towerGroupListIdx, towerGroupIdx].set_ylim(var_yLim)
            if len(var_xTicks) != 0:
                axs[towerGroupListIdx, towerGroupIdx].xaxis.set_ticks(var_xTicks)
            if len(var_yTicks) != 0:
                axs[towerGroupListIdx, towerGroupIdx].yaxis.set_ticks(var_yTicks)
            
            
            # increment the towerIdx counter for the next loop
            towerIdx = towerIdx + 1
        
        
        # remove unused axes
        for ax in axs[towerGroupListIdx,current_nTowers:max_nTowersPerGroup]:
            ax.remove()
    
    
    handles, labels = axs[0,0].get_legend_handles_labels()   ## using the first subplot axes for the full figure legend
    fig.legend(handles, labels, loc="upper right", shadow=True, fontsize=16)   ## might need to adjust the location of this legend
    
    plotDPI = 100
    #plotDPI = 300
    #fig.tight_layout(rect=[0.001,0.001,0.999, 0.999])
    plotFilename = plotOutputDir + "/" + axesLimType + "_" + var_name + ".png"
    plt.savefig( plotFilename, bbox_inches='tight', pad_inches=0.1, dpi=plotDPI )
    
    ## method to just show, manually close each time
    ##plt.show()
    
    ## method to just show for a second, then automatically closes each time
    ##plt.show(block=False)
    ##plt.pause(1)
    ##plt.clf()
    ##plt.cla()
    ##plt.close()
    
    ## probably safest to just cleanup after the plot either way
    plt.clf()
    plt.cla()
    plt.close()
    
    return (  )




def main():
    
    
    ### call the function to set the various directories for use in the script
    nTowers,  plotOutputDir = setScriptDirs( towerList, towerHeights,  plotCasename )
    
    
    ### read in the netcdf tower data
    nTowers, nTowerHeights,  times, timesRange, timesAsDates, timesRangeAsDates, nTimes, nTimesRange,  towerLats, towerLons,  tower_ux, tower_uy, tower_uz, tower_tke, tower_spd, tower_dir = read_perdigaoTowerData( perdigao_ncFile,  towerList, towerHeights,  desiredTimeRangeMin, desiredTimeRangeMax,   plotOutputDir )
    
    
    ### now plot the tower data
    
    ### gonna do 3 sets of plots, repeat with possibly 3 sets of lims, first no lims, 2nd zLimOnly, 3rd the particular axis lims that were input
    ### so need to prepare a set of empty lims to use, and a set of zLimOnly lims as well, the particular lims are already input and set
    
    empty_xLim = []
    empty_yLim = []
    empty_xTicks = []
    empty_yTicks = []
    
    
    ##axesLimType = "noLims"
    ##
    ##plotTowerVsLineSampleData( "spd", "m/s",  axesLimType, empty_xLim, empty_yLim, empty_xTicks, empty_yTicks,  towerGroups, towerGroupNames,  tower_spd,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    ##
    ##plotTowerVsLineSampleData( "dir", "deg",  axesLimType, empty_xLim, empty_yLim, empty_xTicks, empty_yTicks,  towerGroups, towerGroupNames,  tower_dir,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    ##
    ##plotTowerVsLineSampleData( "tke", "m²/s²",  axesLimType, empty_xLim, empty_yLim, empty_xTicks, empty_yTicks,  towerGroups, towerGroupNames,  tower_tke,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    
    
    axesLimType = "zLimOnly"
    
    plotTowerVsLineSampleData( "spd", "m/s",  axesLimType, empty_xLim, spd_yLim, empty_xTicks, spd_yTicks,  towerGroups, towerGroupNames,  tower_spd,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    
    plotTowerVsLineSampleData( "dir", "deg",  axesLimType, empty_xLim, dir_yLim, empty_xTicks, dir_yTicks,  towerGroups, towerGroupNames,  tower_dir,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    
    plotTowerVsLineSampleData( "tke", "m²/s²",  axesLimType, empty_xLim, tke_yLim, empty_xTicks, tke_yTicks,  towerGroups, towerGroupNames,  tower_tke,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    
    
    axesLimType = "inputLims"
    
    plotTowerVsLineSampleData( "spd", "m/s",  axesLimType, spd_xLim, spd_yLim, spd_xTicks, spd_yTicks,  towerGroups, towerGroupNames,  tower_spd,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    
    plotTowerVsLineSampleData( "dir", "deg",  axesLimType, dir_xLim, dir_yLim, dir_xTicks, dir_yTicks,  towerGroups, towerGroupNames,  tower_dir,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    
    plotTowerVsLineSampleData( "tke", "m²/s²",  axesLimType, tke_xLim, tke_yLim, tke_xTicks, tke_yTicks,  towerGroups, towerGroupNames,  tower_tke,  nTowers, nTowerHeights, towerList, towerLats, towerLons, towerHeights,  plotOutputDir, plotCasename )
    
    
    
    
    
    
    
    
    
    
    

if __name__ == '__main__':
    main()



