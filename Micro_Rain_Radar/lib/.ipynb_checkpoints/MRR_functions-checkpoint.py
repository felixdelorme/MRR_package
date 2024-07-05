#######################################################################################################
#######################################################################################################
## 
##raw2snow: Process raw files using Maahn and Kollias method (2012)
##
##
#######################################################################################################
#######################################################################################################

def raw2snow(file_in,file_out, TRES = 60, Descr = "MRR data", author = "APRES3 project",ncForm="NETCDF3_CLASSIC",removeOverlap=False,varsToSave="all"): #convert raw files in MRR Doppler moments
    import os, time, warnings, sys
    import matplotlib.pyplot as plt
    path_IMProToo="/IMProToo/" #See more in https://github.com/maahn/IMProToo
    sys.path.append(path_IMProToo)
    import core3 as IMProToo
    warnings.filterwarnings("ignore") #Ignore IMProToo warmings
    rawData = IMProToo.mrrRawData(file_in,removeOverlap=removeOverlap) # Read RawData
    processedSpec = IMProToo.MrrZe(rawData) #create the IMProToo object and load rawData
    processedSpec.averageSpectra(TRES)# integrates the data in 60 seconds by default.
    processedSpec.co["ncCreator"] = author
    processedSpec.co["ncDescription"] = Descr
    processedSpec.co["dealiaseSpectrum"] = True 
    processedSpec.rawToSnow() # Converts RawData into Radar moments
    processedSpec.writeNetCDF(file_out,ncForm=ncForm,varsToSave=varsToSave) # Saves the processed data in a nc file    

def MRRQlooks(file_in, file_out, Ze_ranges = [-10, 30], W_ranges = [-1, 3], SW_ranges = [0, 2],S_ranges = [0, 3],format='png',dpi=300,name_station = '',max_h=3, height_plot = False, cmap = 'jet', Zecorrected=True, IncludeS=True):
    import numpy as np
    from netCDF4 import Dataset
    import matplotlib.dates as mdate
    import matplotlib.pyplot as plt
    import time

    font = {'family'    :   'serif',
            'weight'    :   'normal',
            'size'      :   22}

    plt.rc('font', **font)        

    ds = Dataset(file_in)
    
    if Zecorrected == True:
        Ze = ds.variables["ZeX"][:]#.T #Ze X band in dB
    else:
        Ze = ds.variables["Ze"][:]#.T #Ze in dB
            
    W = ds.variables["W"][:]#.T #Radial velocity in m/s
    SW = ds.variables["spectralWidth"][:]#.T #spectral width in m/s
    if IncludeS == True: S = ds.variables["SnowfallRate"][:]#.T #S mm/h
    Time = ds.variables["time"][:] # Epoch time seconds since 1970/1/1
    H = ds.variables["height"][:]/1000.#.T #height in km 

    #single vector of heights
    h = np.linspace(np.nanmin(H),np.nanmax(H),31)

    #Mask of time steps without profiles
    H = np.ma.array(H)
    NoDataMask =np.ma.masked_where(H.filled(-1) > 0, np.ones(shape=np.shape(Ze)))
    N_NoData = np.sum(H[:,0].filled(-1)==-1)
    #print('N profiles without Data = ',N_NoData)

    #Get year, month and day
    time_structure = time.gmtime(np.ma.mean(Time)) #using the center of the file as referece
    year,month,day = time_structure.tm_year,time_structure.tm_mon,time_structure.tm_mday

    #Define time sticks for the plots format = 'HH'
    time_ticks_epoch = np.linspace(Time[0],Time[-1],13) 

    time_ticks = [str(time.gmtime(j).tm_hour).zfill(2) for j in time_ticks_epoch]

    #Make figure
    fig = plt.figure(figsize=(18,18))
    
    if IncludeS == True:
        n = 4
        mat = [Ze, W, SW, S]

    else:
        n = 3
        mat = [Ze, W, SW]

    for i in range(1,n+1):
        ax = fig.add_subplot(n, 1, i)
        if (i == 1) & (N_NoData > 0):
            plt.plot([-1,1],[-1,1],'-',color='black',label='No Data', linewidth=2)    
            legend = plt.legend(loc='upper right')
            legend.get_frame().set_edgecolor('1.0')

        if i == 1: vmMn = [Ze_ranges[0], Ze_ranges[1],r'$Z_e$']  
        if i == 2: vmMn = [W_ranges[0], W_ranges[1],r'$W$']    
        if i == 3: vmMn = [SW_ranges[0], SW_ranges[1],r'$\sigma$']    
        if i == 4: vmMn = [S_ranges[0], S_ranges[1],'Snowfall rate']    
        
        plt.pcolormesh(Time,h,np.transpose(NoDataMask), shading='auto', cmap='gray')
        plt.pcolormesh(Time,h,np.transpose(mat[i-1]), vmin = vmMn[0],vmax = vmMn[1],cmap = cmap, shading='auto')

        if i == 1:
            plt.colorbar(label = r'$Z_e$'+" [dB"+r'$Z_e$'+"]")    
            plt.title(name_station+" - MRR, "+str(year).zfill(4)+"-"+str(month).zfill(2)+"-"+str(day).zfill(2))
        if i == 2:
            plt.colorbar(label = r'$W$'+" [m s"+r'$^{-1}$'+"]")    
        if i == 3: 
            plt.colorbar(label = r'$\sigma$'+" [m s"+r'$^{-1}$'+"]")
        if i == 4:
            plt.colorbar(label = 'Snowfall rate'+" [mm h"+r'$^{-1}$'+"]")    

        if i == n: 
            plt.xlabel('Time [UTC]') 
        else: 
            plt.xlabel('')

        plt.ylabel('Height [km a.g.l.]')

        plt.axis([Time[0],Time[-1],0,max_h])
        plt.xticks(time_ticks_epoch,time_ticks)

    plt.savefig(file_out+".png",bbox_inches='tight',format='png',dpi=300)
     
    ds.close()

# Saves Multiple netCDF files
def MFDataset_save(fileout,ds_merged,format="NETCDF3_CLASSIC"):
    from netCDF4 import Dataset # Read and write ncCDF files
    #Create new database
    ds_merged2 = Dataset(fileout, 'w', format=format) 
    ds_merged2.description = ds_merged.description+' Time merged'
    ds_merged2.history = ds_merged.history
    ds_merged2.author = ds_merged.author
    ds_merged2.source = ds_merged.source + "CDA processing"
    ds_merged2.properties = ds_merged.properties

    # dimensions
    ds_merged2.createDimension('time', None)
    ds_merged2.createDimension('range', 31)
    ds_merged2.createDimension('velocity', 192)

    # variables

    time = ds_merged2.createVariable('time', 'i', ('time',), fill_value=-9999)
    range = ds_merged2.createVariable('range', 'i', ('range',),fill_value=-9999)
    velocity = ds_merged2.createVariable('velocity', 'f', ('velocity'),fill_value=-9999.) 
    height = ds_merged2.createVariable('height', 'f', ('time','range'),fill_value=-9999.) 
    Ze = ds_merged2.createVariable('Ze', 'f', ('time','range'), fill_value=-9999.) 
    W = ds_merged2.createVariable('W', 'f', ('time','range'), fill_value=-9999.) 

    #Input Variables

    time[:] = ds_merged.variables['time'][:]
    range[:] = ds_merged.variables['range'][:]
    velocity[:] = ds_merged.variables['velocity'][:]
    height[:] = ds_merged.variables['height'][:]
    Ze[:] = ds_merged.variables['Ze'][:]
    W[:] = ds_merged.variables['W'][:]


    # Variable Attributes

    time.description = ds_merged.variables['time'].description
    range.description = ds_merged.variables['range'].description
    velocity.description = ds_merged.variables['velocity'].description
    height.description = ds_merged.variables['height'].description
    Ze.description = ds_merged.variables['Ze'].description
    W.description = ds_merged.variables['W'].description


    time.units = ds_merged.variables['time'].units
    range.units = ds_merged.variables['range'].units
    velocity.units = ds_merged.variables['velocity'].units
    height.units = ds_merged.variables['height'].units
    Ze.units = ds_merged.variables['Ze'].units
    W.units = ds_merged.variables['W'].units

    time.timezone = ds_merged.variables['time'].timezone    

    varkeys = list(ds_merged.variables.keys())

    if 'quality' in varkeys: 
        quality = ds_merged2.createVariable('quality', 'i', ('time', 'range',),fill_value=-9999)
        quality[:] = ds_merged.variables['quality'][:]
        quality.description = ds_merged.variables['quality'].description
        quality.units = ds_merged.variables['quality'].units

    if 'etaMask' in varkeys: 
        etaMask = ds_merged2.createVariable('etaMask', 'i', ('time', 'range','velocity'),fill_value=-9999)
        etaMask[:] = ds_merged.variables['etaMask'][:]
        etaMask.description = ds_merged.variables['etaMask'].description
        etaMask.units = ds_merged.variables['etaMask'].units

    if 'eta' in varkeys: 
        eta = ds_merged2.createVariable('eta', 'f', ('time','range','velocity'),fill_value=-9999.)  
        eta[:] = ds_merged.variables['eta'][:]
        eta.description = ds_merged.variables['eta'].description
        eta.units = ds_merged.variables['eta'].units
    if 'TF' in varkeys: 
        TF = ds_merged2.createVariable('TF', 'f', ('time','range'), fill_value=-9999.) 
        TF[:] = ds_merged.variables['TF'][:]
        TF.description = ds_merged.variables['TF'].description
        TF.units = ds_merged.variables['TF'].units
    if 'spectralWidth' in varkeys: 
        spectralWidth = ds_merged2.createVariable('spectralWidth', 'f', ('time','range'), fill_value=-9999.) 
        spectralWidth[:] = ds_merged.variables['spectralWidth'][:]
        spectralWidth.description = ds_merged.variables['spectralWidth'].description
        spectralWidth.units = ds_merged.variables['spectralWidth'].units
    if 'skewness' in varkeys: 
        skewness = ds_merged2.createVariable('skewness', 'f', ('time','range'), fill_value=-9999.) 
        skewness[:] = ds_merged.variables['skewness'][:]
        skewness.description = ds_merged.variables['skewness'].description
        skewness.units = ds_merged.variables['skewness'].units

    if 'kurtosis' in varkeys: 
        kurtosis = ds_merged2.createVariable('kurtosis', 'f', ('time','range'), fill_value=-9999.) 
        kurtosis[:] = ds_merged.variables['kurtosis'][:]
        kurtosis.description = ds_merged.variables['kurtosis'].description
        kurtosis.units = ds_merged.variables['kurtosis'].units        
    if 'peakVelLeftBorder' in varkeys: 
        peakVelLeftBorder = ds_merged2.createVariable('peakVelLeftBorder', 'f', ('time', 'range'), fill_value=-9999.) 
        peakVelLeftBorder[:] = ds_merged.variables['peakVelLeftBorder'][:]
        peakVelLeftBorder.description = ds_merged.variables['peakVelLeftBorder'].description
        peakVelLeftBorder.units = ds_merged.variables['peakVelLeftBorder'].units

    if 'peakVelRightBorder' in varkeys: 
        peakVelRightBorder = ds_merged2.createVariable('peakVelRightBorder', 'f', ('time','range'), fill_value=-9999.) 
        peakVelRightBorder[:] = ds_merged.variables['peakVelRightBorder'][:]
        peakVelRightBorder.description = ds_merged.variables['peakVelRightBorder'].description
        peakVelRightBorder.units = ds_merged.variables['peakVelRightBorder'].units
    if 'leftSlope' in varkeys: 
        leftSlope = ds_merged2.createVariable('leftSlope', 'f', ('time','range'), fill_value=-9999.) 
        leftSlope[:] = ds_merged.variables['leftSlope'][:]
        leftSlope.description = ds_merged.variables['leftSlope'].description
        leftSlope.units = ds_merged.variables['leftSlope'].units
    if 'rightSlope' in varkeys: 
        rightSlope = ds_merged2.createVariable('rightSlope', 'f', ('time','range'), fill_value=-9999.) 
        rightSlope[:] = ds_merged.variables['rightSlope'][:]
        rightSlope.description = ds_merged.variables['rightSlope'].description
        rightSlope.units = ds_merged.variables['rightSlope'].units
    if 'etaNoiseAve' in varkeys: 
        etaNoiseAve = ds_merged2.createVariable('etaNoiseAve', 'f', ('time','range'), fill_value=-9999.) 
        etaNoiseAve[:] = ds_merged.variables['etaNoiseAve'][:]
        etaNoiseAve.description = ds_merged.variables['etaNoiseAve'].description
        etaNoiseAve.units = ds_merged.variables['etaNoiseAve'].units
    if 'etaNoiseStd' in varkeys: 
        etaNoiseStd = ds_merged2.createVariable('etaNoiseStd', 'f', ('time','range'), fill_value=-9999.) 
        etaNoiseStd[:] = ds_merged.variables['etaNoiseStd'][:]
        etaNoiseStd.description = ds_merged.variables['etaNoiseStd'].description
        etaNoiseStd.units = ds_merged.variables['etaNoiseStd'].units
    if 'SNR' in varkeys: 
        SNR = ds_merged2.createVariable('SNR', 'f', ('time','range'), fill_value=-9999.) 
        SNR[:] = ds_merged.variables['SNR'][:]
        SNR.description = ds_merged.variables['SNR'].description
        SNR.units = ds_merged.variables['SNR'].units
    if 'ZeX' in varkeys: 
        ZeX = ds_merged2.createVariable('ZeX', 'f', ('time','range'), fill_value=-9999.) 
        ZeX[:] = ds_merged.variables['ZeX'][:]
        ZeX.description = ds_merged.variables['ZeX'].description
        ZeX.units = ds_merged.variables['ZeX'].units
    if 'SnowfallRate' in varkeys: 
        S = ds_merged2.createVariable('SnowfallRate', 'f', ('time','range'), fill_value=-9999.) 
        S[:] = ds_merged.variables['SnowfallRate'][:]
        S.description = ds_merged.variables['SnowfallRate'].description
        S.units = ds_merged.variables['SnowfallRate'].units

    ds_merged2.close()

def merge_files(fnames,fname_out,name_station):
    #from MRR_functions import MRRQlooks # Make the Quicklooks  
    from netCDF4 import MFDataset
    from MRR_functions import MFDataset_save

    import glob

    print("Merging files...")
    ds_merged = MFDataset(fnames)

    print("Saving merged file...")
    MFDataset_save(fname_out,ds_merged,format="NETCDF3_CLASSIC")

    #NameFile_fig = fname_out[:-3]
    #print("Preparing Quicklooks...") #time stamps -> to be improved
    #MRRQlooks(fname_out, NameFile_fig, name_station=name_station, Zecorrected=False, IncludeS=False)
    print("Done.")
