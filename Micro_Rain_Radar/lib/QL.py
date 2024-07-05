def MRRQlooks(file_in, file_out, Ze_ranges = [-20, 20], W_ranges = [-1, 3], SW_ranges = [0, 2],S_ranges = [0, 3],format='png',dpi=300,name_station = '',max_h=3, height_plot = False, cmap = 'jet'):
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
    
    Ze = ds.variables["ZeX"][:] #Ze in dB
    W = ds.variables["W"][:] #Radial velocity in m/s
    SW = ds.variables["spectralWidth"][:] #spectral width in m/s
    S = ds.variables["SnowfallRate"][:] #S mm/h
    Time = ds.variables["time"][:] # Epoch time seconds since 1970/1/1
    H = ds.variables["height"][:]/1000. #height in km 

    mat = [Ze, W, SW, S]

    #single vector of heights
    h = np.linspace(np.nanmin(H),np.nanmax(H),31)

    #Mask of time steps without profiles
    H = np.ma.array(H)
    NoDataMask =np.ma.masked_where(H.filled(-1) > 0, np.ones(shape=np.shape(Ze)))
    N_NoData = np.sum(H[:,0].filled(-1)==-1)
    print('N profiles without Data = ',N_NoData)

    #Get year, month and day
    time_structure = time.gmtime(np.ma.mean(Time)) #using the center of the file as referece
    year,month,day = time_structure.tm_year,time_structure.tm_mon,time_structure.tm_mday

    #Define time sticks for the plots format = 'HH'
    time_ticks_epoch = np.linspace(Time[0],Time[-1],13) 

    time_ticks = [str(time.gmtime(j).tm_hour).zfill(2) for j in time_ticks_epoch]

    #Make figure
    fig = plt.figure(figsize=(18,18))
    
    n = 4
    for i in range(1,5):
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
            plt.xlabel('Time [UTC]') 
            plt.colorbar(label = 'Snowfall rate'+" [mm h"+r'$^{-1}$'+"]")    
            
        else: plt.xlabel('')

        plt.ylabel('Height [km a.g.l.]')

        plt.axis([Time[0],Time[-1],0,max_h])
        plt.xticks(time_ticks_epoch,time_ticks)

    plt.savefig(file_out+".png",bbox_inches='tight',format='png',dpi=300)
     
    ds.close()


