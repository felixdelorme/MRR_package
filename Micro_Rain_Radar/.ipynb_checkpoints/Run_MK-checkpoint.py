
import numpy as np
from  netCDF4 import Dataset
import calendar, datetime, time, glob, os, sys, shutil
from pathlib import Path

sys.path.append("lib/") # adding lib path
sys.path.append("lib/IMProToo/") # adding lib path
sys.path.append("lib/ReportSender/") # adding lib path

from MRR_functions import raw2snow # Process raw data into Doppler moments using MK2012
from MRR_functions import MRRQlooks # Make the Quicklooks  
from MRR_functions import merge_files
import core3 as IMProToo
from emailsender_function import sendreport

np.warnings.filterwarnings('ignore')#to avoid the error messages

df = open('default_parameters.txt','r')

df_lines = df.readlines()

df.close()

#Descriptions
Root_desc = df_lines[0].replace('\t','').replace('\n','').split('=')[0]
OUTpath_desc = df_lines[1].replace('\t','').replace('\n','').split('=')[0]

MODE_desc = df_lines[2].replace('\t','').replace('\n','').split('=')[0]

TRES_desc = df_lines[3].replace('\t','').replace('\n','').split('=')[0]
Short_name_station_desc = df_lines[4].replace('\t','').replace('\n','').split('=')[0]

KtoX_desc = df_lines[5].replace('\t','').replace('\n','').split('=')[0]
KtoX_a_desc = df_lines[6].replace('\t','').replace('\n','').split('=')[0]
KtoX_b_desc = df_lines[7].replace('\t','').replace('\n','').split('=')[0]

ZeToS_desc = df_lines[8].replace('\t','').replace('\n','').split('=')[0]
ZeToS_A_desc = df_lines[9].replace('\t','').replace('\n','').split('=')[0]
ZeToS_B_desc = df_lines[10].replace('\t','').replace('\n','').split('=')[0]

QLooks_desc = df_lines[11].replace('\t','').replace('\n','').split('=')[0]

Overwrite_desc = df_lines[12].replace('\t','').replace('\n','').split('=')[0]

Prefix_desc = df_lines[13].replace('\t','').replace('\n','').split('=')[0]
Sufix_desc = df_lines[14].replace('\t','').replace('\n','').split('=')[0]

removeOverlap_desc = df_lines[15].replace('\t','').replace('\n','').split('=')[0]

selectAll_desc = df_lines[16].replace('\t','').replace('\n','').split('=')[0]


#User values
Root = df_lines[0].replace('\t','').replace('\n','').split('=')[1]
OUTpath = df_lines[1].replace('\t','').replace('\n','').split('=')[1]

MODE = df_lines[2].replace('\t','').replace('\n','').replace(' ','').split('=')[1]

TRES = int(df_lines[3].replace('\t','').replace('\n','').replace(' ','').split('=')[1])
Short_name_station = df_lines[4].replace('\t','').replace('\n','').split('=')[1]

KtoX = df_lines[5].replace('\t','').replace('\n','').replace(' ','').split('=')[1]
KtoX_a = float(df_lines[6].replace('\t','').replace('\n','').replace(' ','').split('=')[1])
KtoX_b = float(df_lines[7].replace('\t','').replace('\n','').replace(' ','').split('=')[1])

ZeToS = df_lines[8].replace('\t','').replace('\n','').replace(' ','').split('=')[1]
ZeToS_A = float(df_lines[9].replace('\t','').replace('\n','').replace(' ','').split('=')[1])
ZeToS_B = float(df_lines[10].replace('\t','').replace('\n','').replace(' ','').split('=')[1])

QLooks = df_lines[11].replace('\t','').replace('\n','').replace(' ','').split('=')[1]

Overwrite = df_lines[12].replace('\t','').replace('\n','').replace(' ','').split('=')[1]

Prefix = df_lines[13].replace('\t','').replace('\n','').replace(' ','').split('=')[1]
Sufix = df_lines[14].replace('\t','').replace('\n','').replace(' ','').split('=')[1]

removeOverlap = df_lines[15].replace('\t','').replace('\n','').replace(' ','').split('=')[1]

selectAll = df_lines[16].replace('\t','').replace('\n','').replace(' ','').split('=')[1]

#Parameters
print("Define parameters (data path, temporal resolution, etc.). YES (Y,y,yes): update parameters, or NOT (N,n,not,Enter): Use default parameters (already defined by the user).")
answer =  input()  #input from the user
if (answer == "Y") or (answer == "y") or (answer == "YES") or (answer == "yes"):

	print('Insert input Data path (Press Enter for default = ' + Root + '):')
	answer =  input()
	if answer != '': Root=answer  #input from the user 

	print('Insert output Data path (Press Enter for default = ' + OUTpath + '):')
	answer =  input()
	if answer != '': OUTpath=answer  #input from the user 

	print('Select Processing MODE (SingleFolder: 1, SubFolders: 2, MergingFiles: 3, netCDF2ASCII: 4) (Press Enter for default = ' + MODE + '):')
	answer =  input()
	if answer != '': MODE=answer  #input from the user 

	print("Insert output temporal resolution in seconds (Press Enter for default = " + str(TRES) + "s):")
	answer =  input()  #input from the user
	if answer != '': TRES =  int(answer)

	print('Insert short name of the station (Press Enter for default = ' + Short_name_station + '):')
	answer =  input()
	if answer != '': Short_name_station = answer

	print('Perform Linear correction of radome attenuation (K to X band conversion)? (True or False) (Press Enter for default = ' + KtoX + '):')
	answer =  input()
	if answer != '': KtoX = answer

	if (KtoX == 'True') or (KtoX == 'T') or (KtoX == 'TRUE') or (KtoX == 'true') or (KtoX == 't'): 
		print('Insert the slope parameter for the radome attenuation correction (Press Enter for default = ' + str(KtoX_a) + '):')
		answer =  input()
		if answer != '': KtoX_a = float(answer)
		print('Insert the intercept parameter for the radome attenuation correction (Press Enter for default = ' + str(KtoX_b) + '):')
		answer =  input()
		if answer != '': KtoX_b = float(answer)

	print('Convert Ze to Precipitation rate? (True or False) (Press Enter for default = ' + ZeToS + '):')
	answer =  input()
	if answer != '': ZeToS = answer

	if (ZeToS == 'True') or (ZeToS == 'T') or (ZeToS == 'TRUE') or (ZeToS == 'true') or (ZeToS == 't'): 
		print('Insert the "A" parameter (constant in Ze-S relationship) (Press Enter for default = ' + str(ZeToS_A) + '):')
		answer =  input()
		if answer != '': ZeToS_A = float(answer)
		print('Insert the "B" parameter (Exponent in Ze-S relationship) (Press Enter for default = ' + str(ZeToS_B) + '):')
		answer =  input()
		if answer != '': ZeToS_B = float(answer)

	print('Do you want to make a quicklook of each output file? (True or False) (Press Enter for default = ' + QLooks + '):')
	answer =  input()
	if answer != '': QLooks = answer

	print('Do you want to Overwrite existing output files? (True or False) (Press Enter for default = ' + Overwrite + '):')
	answer =  input()
	if answer != '': Overwrite = answer

	if MODE == '3':
		print('PREFIX for output files (Only for merged files) (Press Enter for default = ' + Prefix + '):')
		answer =  input()
		if answer != '': Sufix = answer

	print('SUFIX of output files (Press Enter for default = ' + Sufix + '):')
	answer =  input()
	if answer != '': Sufix = answer

	print('Remove next day profiles (after 24h) (True or False) (Press Enter for default = ' + removeOverlap + '):')
	answer =  input()
	if answer != '': removeOverlap = answer

	print('Select All variables. (If False, only "Ze", "W", "specWidth", "SnowfalRate" and "eta" are selected)  (True or False) (Press Enter for default = ' + selectAll + '):')
	answer =  input()
	if answer != '': selectAll = answer

	df = open('default_parameters.txt','w')
	df.write(Root_desc + "\t" + "=" + '\t' + Root + "\n")
	df.write(OUTpath_desc + "\t" + "=" + '\t' + OUTpath + "\n")

	df.write(MODE_desc + "\t" + "=" + '\t' + MODE + "\n")

	df.write(TRES_desc+"\t"+"="+'\t'+str(TRES)+"\n")
	df.write(Short_name_station_desc+"\t"+"="+'\t'+Short_name_station+"\n")

	df.write(KtoX_desc+"\t"+"="+'\t'+KtoX+"\n")
	df.write(KtoX_a_desc+"\t"+"="+'\t'+str(KtoX_a)+"\n")
	df.write(KtoX_b_desc+"\t"+"="+'\t'+str(KtoX_b)+"\n")

	df.write(ZeToS_desc+"\t"+"="+'\t'+ZeToS+"\n")
	df.write(ZeToS_A_desc+"\t"+"="+'\t'+str(ZeToS_A)+"\n")
	df.write(ZeToS_B_desc+"\t"+"="+'\t'+str(ZeToS_B)+"\n")

	df.write(QLooks_desc+"\t"+"="+'\t'+str(QLooks)+"\n")

	df.write(Overwrite_desc+"\t"+"="+'\t'+str(Overwrite)+"\n")

	df.write(Prefix_desc+"\t"+"="+'\t'+str(Prefix)+"\n")
	df.write(Sufix_desc+"\t"+"="+'\t'+str(Sufix)+"\n")

	df.write(removeOverlap_desc+"\t"+"="+'\t'+str(removeOverlap)+"\n")

	df.write(selectAll_desc+"\t"+"="+'\t'+str(selectAll)+"\n")

	df.close()

os.chdir(Root)

if (KtoX == 'True') or (KtoX == 'T') or (KtoX == 'TRUE') or (KtoX == 'true') or (KtoX == 't'):
	KtoX_bool = True
else:
	KtoX_bool = False

if (ZeToS == 'True') or (ZeToS == 'T') or (ZeToS == 'TRUE') or (ZeToS == 'true') or (ZeToS == 't'):
	ZeToS_bool = True
else:
	ZeToS_bool = False

if (removeOverlap == 'True') or (removeOverlap == 'T') or (removeOverlap == 'TRUE') or (removeOverlap == 'true') or (removeOverlap == 't'):
	removeOverlap_bool = True
else:
	removeOverlap_bool = False

if (selectAll == 'True') or (selectAll == 'T') or (selectAll == 'TRUE') or (selectAll == 'true') or (selectAll == 't'):
	selectAll_bool = True
else:
	selectAll_bool = False	

name_station = '_'.join(Short_name_station.split(' '))

Descr = "MRR data at " + name_station + ", first MRR processed with MK12 method v.0.103."

if Sufix != '': Sufix = "_" + Sufix

folder=Root

if (MODE=='1'): 
	dircf=glob.glob(Root+'*.raw')
	input_format = "raw"	
elif (MODE=='2'): 
	dircf=glob.glob(Root+'**/*.raw')
	input_format = "raw"	
elif (MODE in ['3','4']):
	dircf=glob.glob(Root+'*.nc')
	if np.size(dircf)==0:
		dircf=glob.glob(Root+'**/*.nc')

	input_format = "nc"	

else:
	print('MODE '+MODE+' is not available. Check "default_parameter.txt".') 
	sys.exit()

dircf=np.sort(dircf)

if len(dircf) == 1:
	print('In this folder there is '+str(len(dircf)) + ' ' + input_format + ' file')
else:
	print('In this folder there are '+str(len(dircf)) + ' ' + input_format + ' files')

if len(dircf) > 0:
	if not os.path.exists(OUTpath):
		os.mkdir(OUTpath)

if MODE in ["1","2"]:

	for name in dircf:

		NameFile=name 

		count=0

		if MODE == "1":
			NameFile_out = OUTpath+os.path.basename(NameFile)[:-4]+'_'+str(TRES)+'s'+Sufix+'.nc' #create a new file
			NameFile_fig = NameFile_out[:-3]

		if MODE == "2":
			path = Path(NameFile)
			OutFolder = os.path.basename(path.parent)
			OUTpath2 = OUTpath+OutFolder +'/'
			
			if not os.path.exists(OUTpath2):
				os.mkdir(OUTpath2)

			NameFile_out = OUTpath2+os.path.basename(NameFile)[:-4]+'_'+str(TRES)+'s'+Sufix+'.nc' #create a new file
			NameFile_fig = NameFile_out[:-3]

		if selectAll_bool == True:
			varsToSave = 'all'
		else:
			varsToSave = ['Ze','W','specWidth','SNR','quality']	

		if Overwrite == 'True':
			raw2snow(NameFile,NameFile_out, TRES = TRES, Descr = Descr, removeOverlap = removeOverlap_bool,varsToSave=varsToSave) # Convert Raw into Doppler Moments using MK2012
		else:
			if (os.path.isfile(NameFile_out) == False): 
				raw2snow(NameFile,NameFile_out, TRES = TRES, Descr = Descr, removeOverlap = removeOverlap_bool,varsToSave=varsToSave) # Convert Raw into Doppler Moments using MK2012
			else:
				print("The rawfile was not processed.",NameFile_out, "already exists. Change the parameters if you want to Overwrite output files.")

		##Including S estimates using Grazioli et al, 2017.

		if (KtoX == 'True') or (KtoX == 'T') or (KtoX == 'TRUE') or (KtoX == 'true') or (KtoX == 't'): 
			#print("Converting K to X band using Grazioli et al, 2017, TC.")
			ds = Dataset(NameFile_out,'a')
			Ze = ds.variables["Ze"][:]
			try:
				ZeX = ds.createVariable('ZeX', 'f', ('time', 'range',), fill_value=-9999.)
				#print(type(KtoX_a),type(KtoX_b))
				ZeX[:] = KtoX_a*Ze+KtoX_b
				ZeX.description = "Ze converted into X-band to take into accound the radome attenuation (see Grazioli et al. 2017, TC)"
				ZeX.units = "dBZe"
			except:
				print("> Ze radome attenuationed already exists. Change the parameters if you want to Overwrite output files.")

			if (ZeToS == 'True') or (ZeToS == 'T') or (ZeToS == 'TRUE') or (ZeToS == 'true') or (ZeToS == 't'):	
				##print("Including S estimates using Grazioli et al, 2017.")
				try:		
					S = ds.createVariable('SnowfallRate', 'f', ('time', 'range',),fill_value=-9999.)
					S[:] = ((10**(ZeX[:]/10.))/(1.*ZeToS_A))**(1./ZeToS_B) 	
					S.description = "Snowfall rate derived from S-Ze relationship in Grazioli et al. (2017, TC)"
					S.units = "mm h-1"
				except:
					print("> The Snowfall rate already exists. Change the parameters if you want to Overwrite output files.")
				ds.close()    

				if Overwrite == 'True':
					MRRQlooks(NameFile_out, NameFile_fig, name_station=Short_name_station)
				else:
					if (os.path.isfile(NameFile_fig+'.png') == False): 
						MRRQlooks(NameFile_out, NameFile_fig, name_station=Short_name_station)
					else:
						print("> The quicklook already exists. Change the parameters if you want to Overwrite output files.")
				
			else:
				ds.close()    
				if Overwrite == 'True':
					MRRQlooks(NameFile_out, NameFile_fig, name_station=Short_name_station, IncludeS=ZeToS_bool)
				else:
					if (os.path.isfile(NameFile_fig+'.png') == False): 
						MRRQlooks(NameFile_out, NameFile_fig, name_station=Short_name_station, IncludeS=ZeToS_bool)
					else:
						print("> The quicklook was not produced.",NameFile_out, "already exists. Change the parameters if you want to Overwrite output files.")

		else:
			if (ZeToS == 'True') or (ZeToS == 'T') or (ZeToS == 'TRUE') or (ZeToS == 'true') or (ZeToS == 't'):	#In case there is not a Radome only perform the Z-S conversion
				print("Including S estimates using Grazioli et al, 2017, TC.")
				ds = Dataset(NameFile_out,'a')
				Ze = ds.variables["Ze"][:]
				
				varkeys = list(ds.variables.keys())

				if "SnowfallRate" not in varkeys: 

					S = ds.createVariable('SnowfallRate', 'f', ('time', 'range',),fill_value=-9999.)
					S[:] = ((10**(Ze[:]/10.))/(1.*ZeToS_A))**(1./ZeToS_B) 	
					S.description = "Snowfall rate derived from S-Ze relationship in Grazioli et al. (2017, TC)"
					S.units = "mm h-1"
				else:
					print('Variable "SnowfallRate" already exists in the netCDF file. Change the parameters if you want to Overwrite output files.')	
				ds.close()


				if Overwrite == 'True':
					MRRQlooks(NameFile_out, NameFile_fig, name_station=Short_name_station, Zecorrected=KtoX_bool)
				else:
					if (os.path.isfile(NameFile_fig+'.png') == False): 
						MRRQlooks(NameFile_out, NameFile_fig, name_station=Short_name_station, Zecorrected=KtoX_bool)
					else:
						print("> The quicklook was not produced.",NameFile_out, "already exists. Change the parameters if you want to Overwrite output files.")

			else:
				try:
					ds.close()
				except NameError:
					print('')


				if Overwrite == 'True':
					MRRQlooks(NameFile_out, NameFile_fig, name_station=Short_name_station, Zecorrected=KtoX_bool, IncludeS=ZeToS_bool)
				else:
					if (os.path.isfile(NameFile_fig+'.png') == False): 
						MRRQlooks(NameFile_out, NameFile_fig, name_station=Short_name_station, Zecorrected=KtoX_bool, IncludeS=ZeToS_bool)
					else:
						print("> The quicklook was not produced.",NameFile_out, "already exists. Change the parameters if you want to Overwrite output files.")

if MODE in ["3"]: # Mergin All NC files within the directory
	fnames = sorted(dircf) #Files sorted by name
	fname_out = OUTpath + Prefix + Sufix + ".nc"

	print(OUTpath)

	merge_files(fnames,fname_out,name_station=Short_name_station)

	print("All netCDF files were merged in : ", fname_out)


if MODE in ["4"]: #netCDF2ASCII conversion

	for name in dircf:
		fname_out = name[:-3] + ".txt"

		os.system("ncdump " + name + " > "+fname_out)
		print(name, " was converted to ASCII format")

print("done!")

##default Parameters

#Input data path	=	/home/claudio/Projects/git/MK_packaged/Data/RawSpectra/201912/
#Output data path	=	/home/claudio/Projects/git/MK_packaged/Data/MK_processed/
#MODE (SingleFolder: 1, SubFolders: 2, MergingFiles: 3, netCDF2ASCII: 4)	=	1
#output temporal resolution in seconds	=	60
#Short name of the station	=	DDU
#Radome attenuation correction (k to x band) (True or False)	=	True
##Radome attenuation (k to x band), a slope (dBZ) 	=	0.99
##Radome attenuation (k to x band), b intercept (dBZ) 	=	6.14
#Snowfall rate conversion (Z-S)  (True or False) 	=	True
##Z-S relationship, A parameter (constant) 	=	76.0
##Z-S relationship, B parameter (exponent) 	=	0.91
#Make quicklook (True or False)	=	True
#Overwrite (True or False)	=	False
#PREFIX for output files (Only for merged files)	=	201912
#SUFIX for output files	=	MK_Ze-corrected_S-included
#Remove previous/next day profiles (before/after OO/24h) (True or False) 	=	True
#Select All variables. (If False, only "Ze", "W", "specWidth" and "SnowfalRate" are selected)  (True or False) 	=	False