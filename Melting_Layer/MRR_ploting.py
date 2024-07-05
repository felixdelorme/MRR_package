import netCDF4 as nc
import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.dates as mdate
import matplotlib.pyplot as plt
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import astropy.convolution.convolve
from astropy.convolution import Gaussian2DKernel, convolve
""" Félix DELORME on Monday 24 June 2024 """



"""                             REFLECTIVITY                                """

def read_data(file_path, event_start, event_stop):
    ds = nc.Dataset(file_path, mode='r')
    time = ds.variables['time']
    ind_start = np.where(event_start == pd.to_datetime(time[:].data, unit='s'))[0][0]
    ind_stop = np.where(event_stop == pd.to_datetime(time[:].data, unit='s'))[0][0]
    time = time[ind_start:ind_stop]
    Ze = ds.variables['Ze'][ind_start:ind_stop]
    W = ds.variables['W'][ind_start:ind_stop]
    height = ds.variables['height'][ind_start:ind_stop]
    date_event = pd.to_datetime(time[0], unit='s').strftime("%m-%d-%Y")

    return time, Ze, W, height, date_event

def filter_reflectivity(Ze):
    for i in range(len(Ze)):
        if ma.count_masked(Ze[i, :]) / len(Ze[i, :]) > 0.2:
            Ze[i, :] = ma.masked_where(True, Ze[i, :])
    k = Gaussian2DKernel(0.5)
    Ze_filtered = convolve(Ze, k, boundary='fill')
    return ma.array(data=Ze_filtered, mask=Ze.mask, fill_value=-9999)

def plot_maxZe(time, maxZe, date_event):
    plt.plot_date(pd.to_datetime(time[:], unit='s'), maxZe, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('max Ze value')
    plt.xlabel('Time [UTC]')
    plt.title('MRR' + date_event + ' \n Maximum reflectivity')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_maxZe_elevation(time, hmaxZe, date_event):
    plt.plot_date(pd.to_datetime(time[:], unit='s'), hmaxZe, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('max Ze elevation [m]')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Maximum reflectivity elevation')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_filtered_maxZe_elevation(time, hmaxZe, date_event):
    hmaxZe = ma.masked_where(hmaxZe < 400, hmaxZe)
    std_Ze = 1.5
    z_score = np.abs(stats.zscore(ma.filled(hmaxZe, fill_value=np.nan), nan_policy='omit'))
    hmaxZe = ma.masked_where(z_score > std_Ze, hmaxZe)
    plt.plot_date(pd.to_datetime(time[:], unit='s'), hmaxZe, '+', markersize=2, markeredgewidth=1)
    plt.ylim(0, 3000)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Maximum reflectivity position')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_smoothed_maxZe_elevation(time, hmaxZe, date_event):
    hmaxZe_interpol = pd.DataFrame(hmaxZe).interpolate(method='linear', limit=90, limit_area='inside')
    hmaxZe_smooth = gaussian_filter1d(hmaxZe_interpol.values[:, 0], sigma=20)
    hmaxZe_smooth = ma.masked_values(hmaxZe_smooth, np.nan)
    hmaxZe_smooth = ma.masked_array(data=hmaxZe_smooth.data, mask=np.isnan(hmaxZe_smooth.data))
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmaxZe_smooth, '-', label=r'max(Ze)', markersize=3, markeredgewidth=1)
    plt.ylim(0, 3000)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Maximum reflectivity elevation')
    plt.legend(fontsize=10)
    plt.grid(which='major')
    plt.grid(which='minor', alpha=0.5)
    plt.minorticks_on()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.show()

def plot_min_gradZe(time, mingradZe, date_event):
    plt.plot_date(pd.to_datetime(time[:], unit='s'), mingradZe, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('min grad Ze')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Minimum reflectivity gradient value')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_min_gradZe_elevation(time, hmingradZe, date_event):
    plt.plot_date(pd.to_datetime(time[:], unit='s'), hmingradZe, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('height')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Minimum reflectivity gradient position')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_filtered_min_gradZe_elevation(time, hmingradZe, hmaxZe_smooth, date_event):
    hmingradZe = ma.masked_where(hmingradZe < 400, hmingradZe)
    hmingradZe = ma.masked_where(hmingradZe > (hmaxZe_smooth.mean() + 600), hmingradZe)
    std_gradZe = 1
    z_score = np.abs(stats.zscore(ma.filled(hmingradZe, fill_value=np.nan), nan_policy='omit'))
    hmingradZe = ma.masked_where(z_score > std_gradZe, hmingradZe)
    plt.plot_date(pd.to_datetime(time[:], unit='s'), hmingradZe, '+', markersize=2, markeredgewidth=1)
    plt.ylim(0, 3000)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Minimum reflectivity gradient position')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_hmingradZe(time, hmingradZe, hmaxZe_smooth, date_event):

    hmingradZe_interpol = pd.DataFrame(hmingradZe).interpolate(method='linear', limit=90, limit_area='inside')
    hmingradZe_smooth = gaussian_filter1d(hmingradZe_interpol.values[:, 0], sigma=20)
    hmingradZe_smooth = ma.masked_values(hmingradZe_smooth, np.nan)
    hmingradZe_smooth = ma.masked_array(data=hmingradZe_smooth.data, mask=np.isnan(hmingradZe_smooth.data))

    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmaxZe_smooth, '-', label=r'max(Ze)', markersize=3, markeredgewidth=1)
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmingradZe_smooth, '-', label=r'min(gradZe)', markersize=2, markeredgewidth=1)
    plt.ylim(0, 3000)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Min reflectivity elevation')
    plt.legend(fontsize=10)
    plt.grid(which='major')
    plt.grid(which='minor', alpha=0.5)
    plt.minorticks_on()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.rc('xtick', labelsize=12)
    plt.show()

def plot_max_gradZe(time, gradZe, imaxZe, height, date_event):
    gradZe_formax = gradZe.copy()
    for i in range(len(time)):
        gradZe_formax[i, :] = ma.masked_where(np.arange(0, 31) >= imaxZe[i], gradZe[i, :])

    maxgradZe = ma.MaskedArray.max(gradZe_formax, axis=1)
    imaxgradZe = ma.masked_array(ma.MaskedArray.argmax(gradZe_formax, axis=1), fill_value=0, mask=maxgradZe.mask)
    hmaxgradZe = ma.masked_array(height[:, imaxgradZe.data][0], mask=imaxgradZe.mask)

    plt.plot_date(pd.to_datetime(time[:], unit='s'), maxgradZe, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('min grad Ze')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Minimum reflectivity gradient value')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

    plt.plot_date(pd.to_datetime(time[:], unit='s'), hmaxgradZe, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('height')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Minimum reflectivity gradient position')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_max_gradZe(time, gradZe, imaxZe, height, date_event):
    gradZe_formax = gradZe.copy()
    for i in range(len(time)):
        gradZe_formax[i, :] = ma.masked_where(np.arange(0, 31) >= imaxZe[i], gradZe[i, :])

    maxgradZe = ma.MaskedArray.max(gradZe_formax, axis=1)
    imaxgradZe = ma.masked_array(ma.MaskedArray.argmax(gradZe_formax, axis=1), fill_value=0, mask=maxgradZe.mask)
    hmaxgradZe = ma.masked_array(height[:, imaxgradZe.data][0], mask=imaxgradZe.mask)

    plt.plot_date(pd.to_datetime(time[:], unit='s'), maxgradZe, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('min grad Ze')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Minimum reflectivity gradient value')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

    plt.plot_date(pd.to_datetime(time[:], unit='s'), hmaxgradZe, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('height')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Minimum reflectivity gradient position')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_filtered_hmaxgradZe(time, hmaxgradZe, hmaxZe_smooth, date_event):
    hmaxgradZe = ma.masked_where(hmaxgradZe < 400, hmaxgradZe)
    hmaxgradZe = ma.masked_where(hmaxgradZe < (hmaxZe_smooth.mean() - 600), hmaxgradZe)
    hmaxgradZe = ma.masked_where(hmaxgradZe > hmaxZe_smooth, hmaxgradZe)

    from scipy import stats
    std_gradZe = 2
    z_score = np.abs(stats.zscore(ma.filled(hmaxgradZe, fill_value=np.nan), nan_policy='omit'))
    print('remove all values above ' + str(hmaxgradZe.mean() + std_gradZe * hmaxgradZe.std()) +
          ' and below ' + str(hmaxgradZe.mean() - std_gradZe * hmaxgradZe.std()))

    hmaxgradZe = ma.masked_where(z_score > std_gradZe, hmaxgradZe)

    plt.plot_date(pd.to_datetime(time[:], unit='s'), hmaxgradZe, '+', markersize=2, markeredgewidth=1)
    plt.ylim(0, 3000)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n Minimum reflectivity gradient position')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_smoothed_hmaxgradZe(time, hmaxgradZe, hmaxZe_interpol, hmaxZe_smooth, hmingradZe_smooth, date_event):
    hmaxgradZe_interpol = pd.DataFrame(hmaxgradZe).interpolate(method='linear', limit=90, limit_area='inside')
    hmaxgradZe_interpol[hmaxgradZe_interpol >= hmaxZe_interpol] = np.nan
    hmaxgradZe_smooth = gaussian_filter1d(hmaxgradZe_interpol.values[:, 0], sigma=20)
    hmaxgradZe_smooth = ma.masked_values(hmaxgradZe_smooth, np.nan)
    hmaxgradZe_smooth = ma.masked_array(data=hmaxgradZe_smooth.data, mask=np.isnan(hmaxgradZe_smooth.data))

    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmaxZe_smooth, '-', label=r'max(Ze)', markersize=3, markeredgewidth=1)
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmingradZe_smooth, '-', label=r'min(gradZe)', markersize=2, markeredgewidth=1)
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmaxgradZe_smooth, '-', label=r'max(gradZe)', markersize=2, markeredgewidth=1)
    plt.ylim(0, 3000)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title('MRR  ' + date_event + '\n ML elevation from reflectivity')
    plt.legend(fontsize=10)
    plt.grid(which='major')
    plt.grid(which='minor', alpha=0.5)
    plt.minorticks_on()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%d/%m %H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.rc('xtick', labelsize=12)
    plt.show()

def plot_REFspe_time(time, Ze_filtered, gradZe, height, hmaxZe, hmaxgradZe, specific_time):
    id_sat_time = np.argwhere(pd.to_datetime(time[:].data, unit='s') == specific_time)[0][0]
    Ze_avg = Ze_filtered[id_sat_time-5:id_sat_time+5,:].mean(axis=0)

    fig = plt.figure(facecolor='white', figsize=(5,4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_ylim(0, 2000)
    ax1.set_xlim(0, 35)
    ax2.set_xlim(-4, 6)
    p1 = ax1.plot(Ze_avg, height[id_sat_time,:], '-+', c='b', label='Ze')
    p2 = ax2.plot(gradZe[id_sat_time-5:id_sat_time+5,:].mean(axis=0), height[id_sat_time,:], '-+', c='r', label=r'$\partial_{z}{Ze}$')
    ax2.set_xlabel(r"$10^{-2}.dBZ.m^{-1}$", color='r')
    ax2.tick_params(axis='x', labelcolor='r')
    ax1.set_xlabel("dBZ", color='b')
    ax1.tick_params(axis='x', labelcolor='b')
    ax1.set_ylabel("height [m.a.g.l]")
    h1 = ax1.axhline(hmaxZe[id_sat_time], ls='--', alpha=0.7, label='max(Ze)')
    h2 = ax1.axhline(hmaxgradZe[id_sat_time], ls='--', color='g', alpha=0.7, label=r'max($\partial_{z}{Ze}$)')
    h3 = ax1.axhline(1100, ls='--', color='orange', alpha=0.7, label=r'min($\partial_{z}{Ze}$)')
    plt.legend(handles=p1+p2+[h1]+[h2]+[h3])
    ax1.grid()
    fig.tight_layout()
    plt.show()


"""                             VELOCITY                                  """

def remove_invalid_time_steps(W_filt):
    for i in range(len(W_filt)):
        if ma.count_masked(W_filt[i, :]) / len(W_filt[i, :]) > 0.5:
            W_filt[i, :] = ma.masked_where(True, W_filt[i, :])
    return W_filt

def filter_gradient(gradW):
    maxgradW = ma.MaskedArray.max(gradW, axis=1)
    maxgradW = ma.masked_where((maxgradW < 0.25) | (maxgradW > 4), maxgradW)
    imaxgrad = ma.masked_array(ma.MaskedArray.argmax(gradW, axis=1), fill_value=0, mask=maxgradW.mask)
    return maxgradW, imaxgrad

def plot_max_gradient(time, hmaxgrad, date_event):
    plt.plot_date(pd.to_datetime(time, unit='s'), hmaxgrad, '+', markersize=2, markeredgewidth=1)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title(f'MRR  {date_event} \n Maximum velocity gradient position')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def remove_outliers(hmaxgrad, std_gradW=1):
    z_score = np.abs(stats.zscore(ma.filled(hmaxgrad, fill_value=np.nan), nan_policy='omit'))
    return ma.masked_where(z_score > std_gradW, hmaxgrad)

def maxcompute_second_gradient(gradW, height, h_smooth):
    grad2W = -np.gradient(gradW, axis=1)
    hML_lim = h_smooth.mean() + 500
    maxgrad2W = ma.MaskedArray.max(ma.masked_where(height > hML_lim, grad2W), axis=1)
    imaxgrad2W = ma.masked_array(ma.MaskedArray.argmax(ma.masked_where(height > hML_lim, grad2W), axis=1), fill_value=0, mask=maxgrad2W.mask)
    hmaxgrad2W = ma.masked_array(height[:, imaxgrad2W.data][0], mask=imaxgrad2W.mask)
    return hmaxgrad2W

def mincompute_second_gradient(gradW, height, h_smooth):
    grad2W = -np.gradient(gradW, axis=1)
    hML_lim = h_smooth.mean() + 500
    mingrad2W= ma.MaskedArray.min(ma.masked_where(height>hML_lim,grad2W),axis=1)
    mingrad2W= ma.masked_where((np.abs(mingrad2W)<0.2) | (np.abs(mingrad2W)>3),mingrad2W) #only keep the strong gradients
    imingrad2W=ma.masked_array(ma.MaskedArray.argmin(ma.masked_where(height>hML_lim,grad2W),axis=1),
    fill_value=0, mask=mingrad2W.mask) #array of index of the maxima
    hmingrad2W= ma.masked_array(height[:,imingrad2W.data][0],mask=imingrad2W.mask)
    hmingrad2W= ma.masked_where(hmingrad2W<=500, hmingrad2W)
    return hmingrad2W

def plot_max_second_gradient(time, hmaxgrad2W, date_event):
    plt.plot_date(pd.to_datetime(time, unit='s'), hmaxgrad2W, '+', markersize=2, markeredgewidth=1)
    plt.ylim(0, 3000)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title(f'MRR  {date_event} \n Maximum second velocity gradient position')
    plt.grid()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def remove_second_gradient_outliers(hmaxgrad2W, std_grad2W=1.5):
    z_score2 = np.abs(stats.zscore(ma.filled(hmaxgrad2W, fill_value=np.nan), nan_policy='omit'))
    return ma.masked_where(z_score2 > std_grad2W, hmaxgrad2W)

def maxsmooth_second_gradient(hmaxgrad2W,h_smooth, sigma=20):
    hmaxgrad2W_interpol = pd.DataFrame(hmaxgrad2W).interpolate(method='linear', limit=90, limit_area='inside')
    hmaxgrad2W_smooth = gaussian_filter1d(hmaxgrad2W_interpol.values[:, 0], sigma=sigma)
    hmaxgrad2W_smooth = ma.masked_values(hmaxgrad2W_smooth, np.nan)
    hmaxgrad2W_smooth = ma.masked_array(data=hmaxgrad2W_smooth.data, mask=np.isnan(hmaxgrad2W_smooth.data))
    return ma.masked_where(hmaxgrad2W_smooth < h_smooth, hmaxgrad2W_smooth)

def minsmooth_second_gradient(hmingrad2W,h_smooth, sigma=20):
    hmingrad2W_interpol=pd.DataFrame(hmingrad2W).interpolate(method='linear', limit=90, limit_area='inside')
    hmingrad2W_smooth = gaussian_filter1d(hmingrad2W_interpol.values[:,0],sigma=sigma)
    hmingrad2W_smooth = ma.masked_values(hmingrad2W_smooth, np.nan)
    hmingrad2W_smooth = ma.masked_array(data=hmingrad2W_smooth.data, mask=np.isnan(hmingrad2W_smooth.data))
    return ma.masked_where(hmingrad2W_smooth>h_smooth, hmingrad2W_smooth)

def plot_combined_gradients(time, h_smooth, hmaxgrad2W_smooth, hmingrad2W_smooth, date_event):
    plt.plot_date(pd.to_datetime(time, unit='s'), h_smooth, '-', label=r'max($\nabla W$)', markersize=2, markeredgewidth=1)
    plt.plot_date(pd.to_datetime(time, unit='s'), hmaxgrad2W_smooth, '-', label=r'max($\nabla^{2} W$)', markersize=1, markeredgewidth=1)
    plt.plot_date(pd.to_datetime(time, unit='s'), hmingrad2W_smooth, '-', label=r'min($\nabla^{2} W$)', markersize=1, markeredgewidth=1)
    plt.ylim(0, 3000)
    plt.xlim(pd.to_datetime(time[0], unit='s'), pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Height (m.a.i.l)')
    plt.xlabel('Time [UTC]')
    plt.title(r'MRR ' + date_event + r'\n $\partial_{z}{W}$ and $\partial_{z}^{2}{W}$ extrema positions')
    plt.legend(fontsize=10)
    plt.grid(which='major')
    plt.grid(which='minor', alpha=0.5)
    plt.minorticks_on()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%H:00')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.show()

def plot_VELspe_time(time, W_filt, gradW, height, hmaxgrad, hmaxgrad2W, hmingrad2W, specific_time):

    id_sat_time = np.argwhere(pd.to_datetime(time[:].data, unit='s') == specific_time)[0][0]
    grad2W= -np.gradient(gradW, axis=1)
    W_Avg = W_filt[id_sat_time-5:id_sat_time+5,:].mean(axis=0)

    fig = plt.figure(facecolor='white', figsize=(5,5))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax3 = ax1.twiny()
    ax3.spines['top'].set_position(('outward', 30))
    ax3.set_xlim(-1, 2)
    ax2.set_xlim(-0.3, 2.7)
    ax1.set_xlim(0, 6)
    ax1.set_ylim(0, 2000)
    p1=ax1.plot(W_Avg, height[id_sat_time,:], '-+',c='b', label='W')
    p2=ax2.plot(gradW[id_sat_time,:], height[id_sat_time,:], '-+',c='r', label=r'$\partial_{z}{W}$')
    p3=ax3.plot(grad2W[id_sat_time,:], height[id_sat_time,:], '-+',c='purple', label=r'$\partial_{z}^2{W}$')
    ax2.set_xlabel(r"$10^{-2}.s^{-1}$", color='r')
    ax2.tick_params(axis='x', labelcolor='r')
    ax1.set_xlabel(r"$m.s^{-1}$", color='b')
    ax1.tick_params(axis='x', labelcolor='b')
    ax3.set_xlabel(r"$10^{-4}.m^{-1}.s^{-1}$", color='purple')
    ax3.tick_params(axis='x', labelcolor='purple')
    ax1.set_ylabel("height [m.a.g.l]")
    h1 = ax1.axhline(hmaxgrad[id_sat_time], color='b', linestyle='-', label=r'max($\partial_{z}W$)')
    h2 = ax1.axhline(hmaxgrad2W[id_sat_time], color='r', linestyle='-', label=r'max($\partial_{z}^2 W$)')
    h3 = ax1.axhline(hmingrad2W[id_sat_time], color='purple', linestyle='-', label=r'min($\partial_{z}^2 W$)')
    ax1.grid()
    plt.legend(handles=p1+p2+p3+[h1]+[h2]+[h3], ncol=2)
    plt.show()

def velocity_reflectivity_boundaries(time, h_smooth, hmaxgrad2W_smooth, hmingrad2W_smooth, hmaxZe_smooth, hmingradZe_smooth, hmaxgradZe_smooth, date_event):
    alt_station = 0 # meter above sea level ( a.s.l ) 
    plt.figure(facecolor='white', figsize=(8,5))
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), h_smooth + alt_station, '-', label=r'$max(\partial_z{W})$', markersize=1, markeredgewidth=1, alpha=0.9)
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmaxgrad2W_smooth + alt_station, '-', label=r'$max(\partial_z^2{W})$', markersize=1, markeredgewidth=1, alpha=0.9)
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmingrad2W_smooth + alt_station, '-', label=r'$min(\partial_z{W})$', markersize=1, markeredgewidth=1, alpha=0.9)
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmaxZe_smooth + alt_station, '--', label=r'$max(Ze)$', color='tab:blue', markersize=1, markeredgewidth=1, alpha=0.9)
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmingradZe_smooth + alt_station, '--', label=r'$min(\partial_z{Ze})$', color='tab:orange', markersize=1, markeredgewidth=1, alpha=0.9)
    plt.plot_date(pd.to_datetime(time[:].data, unit='s'), hmaxgradZe_smooth + alt_station, '--', label=r'$max(\partial_z{Ze})$', color='tab:green', markersize=1, markeredgewidth=1, alpha=0.9)
    plt.ylim(0+alt_station, 3100+alt_station)
    plt.xlim(pd.to_datetime(time[0], unit='s'),pd.to_datetime(time[-1], unit='s'))
    plt.ylabel('Elevation (m.a.s.l)')
    plt.xlabel('Time [UTC]')
    plt.title('MRR '+ date_event +'\n ML boundaries elevation')
    plt.legend(fontsize=10, ncol=2)
    plt.grid(which='minor', alpha=0.1)
    plt.grid(which='major', alpha=0.3)
    plt.minorticks_on()
    plt.gcf().autofmt_xdate()
    date_format = mdate.DateFormatter('%d%b \n %H:%M')
    plt.gca().xaxis.set_major_formatter(date_format)
    plt.tight_layout()
    plt.rc('xtick', labelsize=10) 
    plt.show()



"""                               MAIN                                  """

def main():
    file_path = input("Entrez le chemin du fichier à lire: ")
    event_start = input("Entrez la date et l'heure de début de l'événement (YYYY-MM-DD HH:MM:SS): ")
    event_stop = input("Entrez la date et l'heure de fin de l'événement (YYYY-MM-DD HH:MM:SS): ")

    # Read and process reflectivity data

    time, Ze, W, height, date_event = read_data(file_path, event_start, event_stop)
    Ze_filtered = filter_reflectivity(Ze) 
    maxZe = ma.MaskedArray.max(Ze_filtered, axis=1)
    imaxZe = ma.masked_array(ma.MaskedArray.argmax(Ze_filtered, axis=1), fill_value=0, mask=maxZe.mask)
    hmaxZe = ma.masked_array(height[:, imaxZe.data][0], mask=imaxZe.mask)

    gradZe = np.gradient(Ze_filtered, axis=1)
    mingradZe = ma.MaskedArray.min(ma.masked_where(height > hmaxZe.mean() + 600, gradZe[:,:]), axis=1)
    imingradZe = ma.masked_array(ma.MaskedArray.argmin(gradZe[:,:], axis=1), fill_value=0, mask=mingradZe.mask)
    hmingradZe = ma.masked_array(height[:, imingradZe.data][0], mask=imingradZe.mask)

    hmaxZe_interpol = pd.DataFrame(hmaxZe).interpolate(method='linear', limit=90, limit_area='inside')
    hmaxZe_smooth = gaussian_filter1d(hmaxZe_interpol.values[:, 0], sigma=20)
    hmaxZe_smooth = ma.masked_values(hmaxZe_smooth, np.nan)
    hmaxZe_smooth = ma.masked_array(data=hmaxZe_smooth.data, mask=np.isnan(hmaxZe_smooth.data))

    hmingradZe_interpol = pd.DataFrame(hmingradZe).interpolate(method='linear', limit=90, limit_area='inside')
    hmingradZe_smooth = gaussian_filter1d(hmingradZe_interpol.values[:, 0], sigma=20)
    hmingradZe_smooth = ma.masked_values(hmingradZe_smooth, np.nan)
    hmingradZe_smooth = ma.masked_array(data=hmingradZe_smooth.data, mask=np.isnan(hmingradZe_smooth.data))

    gradZe_formax = gradZe
    for i in range(0, len(time)):
        gradZe_formax[i,:] = ma.masked_where(np.arange(0,31) >= imaxZe[i], gradZe[i,:])
    maxgradZe = ma.MaskedArray.max(gradZe_formax[:,:], axis=1)
    imaxgradZe = ma.masked_array(ma.MaskedArray.argmax(gradZe_formax[:,:], axis=1), fill_value=0, mask=maxgradZe.mask)
    hmaxgradZe = ma.masked_array(height[:, imaxgradZe.data][0], mask=imaxgradZe.mask)

    hmaxgradZe_interpol = pd.DataFrame(hmaxgradZe).interpolate(method='linear', limit=90, limit_area='inside')
    hmaxgradZe_interpol[hmaxgradZe_interpol >= hmaxZe_interpol] = np.nan
    hmaxgradZe_smooth = gaussian_filter1d(hmaxgradZe_interpol.values[:,0],sigma=20)
    hmaxgradZe_smooth = ma.masked_values(hmaxgradZe_smooth, np.nan)
    hmaxgradZe_smooth = ma.masked_array(data=hmaxgradZe_smooth.data, mask=np.isnan(hmaxgradZe_smooth.data))

    # Read and process velocity data

    W_filt = ma.masked_outside(W, -6, 25)
    W_filt = remove_invalid_time_steps(W_filt)

    k = Gaussian2DKernel(0.5)
    z_filtered = astropy.convolution.convolve(W_filt, k, boundary='fill')
    z_filtered = ma.array(data=z_filtered, mask=W_filt.mask, fill_value=-9999)

    gradW = -np.gradient(z_filtered, axis=1)
    maxgradW, imaxgrad = filter_gradient(gradW)
    hmaxgrad = ma.masked_array(height[:, imaxgrad.data][0], mask=imaxgrad.mask)
    hmaxgrad = remove_outliers(hmaxgrad)
    
    h_interpol = pd.DataFrame(hmaxgrad).interpolate(method='linear', limit=90, limit_area='inside')
    h_smooth = gaussian_filter1d(h_interpol.values[:, 0], sigma=20)
    h_smooth = ma.masked_values(h_smooth, np.nan)
    h_smooth = ma.masked_array(data=h_smooth.data, mask=np.isnan(h_smooth.data))

    hmaxgrad2W = maxcompute_second_gradient(gradW, height, h_smooth)
    hmaxgrad2W = remove_second_gradient_outliers(hmaxgrad2W)
    hmaxgrad2W_smooth = maxsmooth_second_gradient(hmaxgrad2W,h_smooth, sigma=20)

    hmingrad2W = mincompute_second_gradient(gradW, height, h_smooth)
    hmingrad2W = remove_second_gradient_outliers(hmingrad2W)
    hmingrad2W_smooth = minsmooth_second_gradient(hmingrad2W,h_smooth, sigma=20)

    while True:
        print("1: Reflectivité")
        print("2: Velocity")
        print("3: All")
        print("4: Essential")
        print("5: Quitter")

        menu_choice = input("Choisissez un menu: ")

        if menu_choice == '1':
            reflectivity_options = {
                "1": ("Maximum de Ze", plot_maxZe),
                "2": ("Maximum de maxZe_elevation", plot_maxZe_elevation),
                "3": ("Filtrage et tracé de maxZe_elevation sans valeurs aberrantes", plot_filtered_maxZe_elevation),
                "4": ("Interpolation et lissage de maxZe_elevation", plot_smoothed_maxZe_elevation),
                "5": ("Afficher la température", plot_min_gradZe),
                "6": ("Afficher les précipitations et la température", plot_min_gradZe_elevation),
                "7": ("Interpolation et lissage de hmingradZe", plot_hmingradZe),
                "8": ("Gradient maximum de gradZe", plot_max_gradZe),
                "9": ("Filtrage et tracé de hmaxgradZe sans valeurs aberrantes", plot_filtered_hmaxgradZe),
                "10": ("Interpolation et lissage de hmaxgradZe", plot_smoothed_hmaxgradZe),
                "11": ("Profile de Reflectivity à un temps particulier", None)
            }

            print("Options de Reflectivité disponibles:")
            for key, (description, _) in reflectivity_options.items():
                print(f"{key}: {description}")

            choices = input("Choisissez une ou plusieurs options (séparées par des virgules): ").split(',')

            for choice in choices:
                if choice in reflectivity_options:
                    description, func = reflectivity_options[choice]
                    if choice == '11':
                        specific_time = input("Entrez l'horaire spécifique (HH:MM:SS): ")
                        specific_datetime = pd.to_datetime(f"{event_start.split(' ')[0]} {specific_time}")
                        if pd.to_datetime(event_start) <= specific_datetime <= pd.to_datetime(event_stop):
                            plot_REFspe_time(time, Ze_filtered, gradZe, height, hmaxZe, hmaxgradZe, specific_datetime)
                        else:
                            print("L'horaire spécifique est en dehors de l'intervalle de l'événement.")
                    else:
                        func(time, maxZe, date_event)
                else:
                    print(f"Option {choice} non reconnue.")

        elif menu_choice == '2':
            velocity_options = {
                "1": ("Maximum du gradient de velocity en position", plot_max_gradient),
                "2": ("Maximum de la derivée de velocity en position", plot_max_second_gradient),
                "3": ("Combinaison des gradient de velocity", plot_combined_gradients),
                "4": ("Profile de Velocity à un temps particulier", None)
            }

            print("Options de Velocity disponibles:")
            for key, (description, _) in velocity_options.items():
                print(f"{key}: {description}")

            choices = input("Choisissez une ou plusieurs options (séparées par des virgules): ").split(',')

            for choice in choices:
                if choice in velocity_options:
                    description, func = velocity_options[choice]
                    if choice == '4':
                        specific_time = input("Entrez l'horaire spécifique (HH:MM:SS): ")
                        specific_datetime = pd.to_datetime(f"{event_start.split(' ')[0]} {specific_time}")
                        if pd.to_datetime(event_start) <= specific_datetime <= pd.to_datetime(event_stop):
                            plot_VELspe_time(time, W_filt, gradW, height, hmaxgrad, hmaxgrad2W, hmingrad2W, specific_datetime)
                        else:
                            print("L'horaire spécifique est en dehors de l'intervalle de l'événement.")
                    else:
                        func(time, maxgradW, date_event)
                else:
                    print(f"Option {choice} non reconnue.")

        elif menu_choice == '3':
            # Plot all reflectivity figures
            plot_maxZe(time, maxZe, date_event)
            plot_maxZe_elevation(time, hmaxZe, date_event)
            plot_filtered_maxZe_elevation(time, hmaxZe, date_event)
            plot_smoothed_maxZe_elevation(time, hmaxZe_smooth, date_event)
            plot_min_gradZe(time, mingradZe, date_event)
            plot_min_gradZe_elevation(time, hmingradZe, date_event)
            plot_filtered_min_gradZe_elevation(time, hmingradZe, hmaxZe_smooth, date_event)
            plot_max_gradZe(time, gradZe, imaxZe, height, date_event)
            plot_filtered_hmaxgradZe(time, hmaxgradZe, hmaxZe_smooth, date_event)
            plot_smoothed_hmaxgradZe(time, hmaxgradZe, hmaxZe_interpol, hmaxZe_smooth, hmingradZe_smooth, date_event)

            # Plot all velocity figures
            plot_max_gradient(time, hmaxgrad, date_event)
            plot_max_second_gradient(time, hmaxgrad2W, date_event)
            plot_combined_gradients(time, h_smooth, hmaxgrad2W_smooth, hmingrad2W_smooth, date_event)

            specific_time = input("Entrez l'horaire spécifique pour les profils de velocity (HH:MM:SS): ")
            specific_datetime = pd.to_datetime(f"{event_start.split(' ')[0]} {specific_time}")
            if pd.to_datetime(event_start) <= specific_datetime <= pd.to_datetime(event_stop):
                plot_REFspe_time(time, Ze_filtered, gradZe, height, hmaxZe, hmaxgradZe, specific_datetime)
                plot_VELspe_time(time, W_filt, gradW, height, hmaxgrad, hmaxgrad2W, hmingrad2W, specific_datetime)  
            else:
                print("L'horaire spécifique est en dehors de l'intervalle de l'événement.")

            velocity_reflectivity_boundaries(time, h_smooth, hmaxgrad2W_smooth, hmingrad2W_smooth, hmaxZe_smooth, hmingradZe_smooth, hmaxgradZe_smooth, date_event)

        elif menu_choice == '4':

            velocity_reflectivity_boundaries(time, h_smooth, hmaxgrad2W_smooth, hmingrad2W_smooth, hmaxZe_smooth, hmingradZe_smooth, hmaxgradZe_smooth, date_event)
            specific_time = input("Entrez l'horaire spécifique pour les profils de velocity (HH:MM:SS): ")
            specific_datetime = pd.to_datetime(f"{event_start.split(' ')[0]} {specific_time}")
            if pd.to_datetime(event_start) <= specific_datetime <= pd.to_datetime(event_stop):
                plot_REFspe_time(time, Ze_filtered, gradZe, height, hmaxZe, hmaxgradZe, specific_datetime)
                plot_VELspe_time(time, W_filt, gradW, height, hmaxgrad, hmaxgrad2W, hmingrad2W, specific_datetime)  
            else:
                print("L'horaire spécifique est en dehors de l'intervalle de l'événement.")

        elif menu_choice == '5':
            print("Quitter le programme.")
            break
        else:
            print("Choix de menu non valide. Veuillez réessayer.")

if __name__ == "__main__":
    main()

