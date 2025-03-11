#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 15:14:29 2021

@author: magnaldom

Algorithm that
- detects clear skies in obs and model
- calculates the bias and standard deviation for each case from ground measurements
- works using the hour-by-hour method
- detection with SAT data for observations, and cloud cover in AROME
- on zones and not point by point comparison
"""
### ============================================================================
### LIBRARY AND TOOLS IMPORTATION

from scipy import stats
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfea
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib import cm
import netCDF4 as nc
from netCDF4 import Dataset
import time
import array as arr
from array import *
from sklearn.metrics import mean_squared_error


import sys; sys.path.insert(0,'../tools')
from TOA import TOA_moyen_heure
from plot_AROME import zone_M_N_lon_lat_opt
from zenith_angle import zenith_angle_ephem
from get_solar_constant import get_solar_constant
from plot_BDClim_v2 import plot_BDClim_station_h
from num_jour_between import num_jour_between
import sys; sys.path.insert(0,'../cloud_detection')
from detection_cc_AROME_SAT_opt import dcc_AROME_SAT_cumul_zone
from detection_cc_CEMS import dcc_CEMS_cumul_moypoints, fournearestpoints






start = time.time()

### ============================================================================
### DONNEES ENTREES

date_temp1 = sys.argv[1]
date_temp2 = sys.argv[2]
y1 = int(date_temp1[:4])
m1 = int(date_temp1[4:6])
d1 = int(date_temp1[6:8])
y2 = int(date_temp2[:4])
m2 = int(date_temp2[4:6])
d2 = int(date_temp2[6:8])
date1 = datetime(y1,m1,d1,0,0,0)
date2 = datetime(y2,m2,d2,0,0,0)

K = int(sys.argv[3]) #Defaults value : 2
fn_lim = float( sys.argv[4]) #Defaults value : 0.02


### ============================================================================
### LECTURE DONNES RELATIVES AUX STATIONS

fichier_station = open("/cnrm/phynh/data1/magnaldom/BDClim/station_glo_v2" ,"r")
datas = np.loadtxt(fichier_station, dtype = 'str', delimiter=' ', usecols=(0,1,2, 3, 4, 5, 6, 7, 8), unpack=True)
fichier_station.close()
num_station = datas[0]
lon = [float(i) for i in datas[1]]
lat = [float(i) for i in datas[2]]
alti = [float(i) for i in datas[3]]
M = [float(i) for i in datas[5]]
N = [float(i) for i in datas[6]]
M_sat = [float(i) for i in datas[7]]
N_sat = [float(i) for i in datas[8]]
classe = [float(i) for i in datas[4]]
nb_station = len(num_station)

### ============================================================================
### DEFINITIONS MATRICES DE SCORES

score = np.zeros((nb_station,4))
score_pos = np.zeros((nb_station,4))
score_neg = np.zeros((nb_station,4))
ecart_type = np.zeros((nb_station,4))
biais = np.zeros((nb_station,4))
biais_pos = np.zeros((nb_station,4))
biais_neg = np.zeros((nb_station,4))
biais_norm = np.zeros((nb_station,4))
flu_moy = np.zeros((nb_station,4))
flu_moy_norm = np.zeros((nb_station,4))
MAE = np.zeros((nb_station,4))
MAPE = np.zeros((nb_station,4))
flu_moy_mod = np.zeros((nb_station,4))
flu_moy_mod_norm = np.zeros((nb_station,4))
std_deviation = np.zeros((nb_station,4))
std_deviation_norm = np.zeros((nb_station,4))
nb_heure = np.zeros((nb_station,4))
array_cas1 = np.empty((0,4), float) #arr.array('f')
array_cas2 = np.empty((0,4), float)
array_cas3 = np.empty((0,4), float)
array_cas4 = np.empty((0,4), float)




distribution_1 = []
distribution_2 = []
distribution_3 = []
distribution_4 = []

## -------------------------------
### DEBUT DE LA BOUCLE TEMPORELLE/ BEGINNING OF TEMPORAL LOOP

num_jour = num_jour_between(date1,date2)
for d in range(0,num_jour):

	#----------------------------------------
	#--- LECTURE JOUR ET FICHIER AROME/AROME FILES LECTURE-------
	#----------------------------------------
	date_v = date1 + timedelta(days = d)
	annee = date_v.year
	mois = date_v.month
	jour = date_v.day
	if len(expe)>1:
		fichier = "/cnrm/phynh/data1/magnaldom/STOCKAGE_AROME_UP/%s_light/%s%s/AROME_%s%s%s_light.nc" %(expe,str(date_v.year).zfill(4),str(date_v.month).zfill(2), str(date_v.year).zfill(4),str(date_v.month).zfill(2),str(date_v.day).zfill(2))
	else :
		fichier = "/cnrm/phynh/data1/magnaldom/STOCKAGE_AROME_UP/DATA_light/%s%s/AROME_%s%s%s_light.nc" %(str(date_v.year).zfill(4),str(date_v.month).zfill(2), str(date_v.year).zfill(4),str(date_v.month).zfill(2),str(date_v.day).zfill(2))
	arome = Dataset(fichier, "r", format="NETCDF4")
	h_min = 6
#####----------------- Loop on hours
	for h in range(h_min, 23):
		nebul_cems_h = np.zeros((nb_station)) #Save data for each station each hour (initialisation chaque heure)
		csd_arome_h = [np.inf for st in range(nb_station )] #Nan by default
		csd_cems_h = [np.inf for st in range(nb_station )] #Nan by default
		for m in ['00', '15', '30', '45']: #Satellite product frequency
			#------------------------------
			#------ LECTURE CEMS-----------
			#------------------------------
			if m != '00':
				h_cems = h-1
			else :
				h_cems = h
			try :
				if date_v.month == 1 :
					try :
						file = 'S_NWC_%s_MSG4_globeM-VISIR_%s%s%sT%s%s00Z.nc.v2018'%("CMA", str(date_v.year), str(date_v.month).zfill(2), str(date_v.day).zfill(2),str(h_cems).zfill(2), m)
						fichier = "/cnrm/phynh/data1/magnaldom/CEMS/DATA/%s_%s/%s/%s"  %(str(date_v.year), str(date_v.month).zfill(2), 'CMA', file )
						CEMS = Dataset(fichier, "r", format="NETCDF4")
						#print("Trouvé : " + fichier)
					except :
						file = 'S_NWC_%s_MSG2_globeM-VISIR_%s%s%sT%s%s00Z.nc.v2018'%("CMA", str(date_v.year), str(date_v.month).zfill(2), str(date_v.day).zfill(2),str(h_cems).zfill(2), m)
						fichier = "/cnrm/phynh/data1/magnaldom/CEMS/DATA/%s_%s/%s/%s"  %(str(date_v.year), str(date_v.month).zfill(2), 'CMA', file )
						CEMS = Dataset(fichier, "r", format="NETCDF4")
						#print("Trouvé : " + fichier)
				if date_v.month == 3 :
					try :
						file = 'S_NWC_%s_MSG4_globeM-VISIR_%s%s%sT%s%s00Z.nc'%("CMA", str(date_v.year), str(date_v.month).zfill(2), str(date_v.day).zfill(2),str(h_cems).zfill(2), m)
						fichier = "/cnrm/phynh/data1/magnaldom/CEMS/DATA/%s_%s/%s/%s"  %(str(date_v.year), str(date_v.month).zfill(2), 'CMA', file )
						CEMS = Dataset(fichier, "r", format="NETCDF4")
						#print("Trouvé : " + fichier)
					except :
						file = 'S_NWC_%s_MSG4_globeM-VISIR_%s%s%sT%s%s00Z.nc.v2018'%("CMA", str(date_v.year), str(date_v.month).zfill(2), str(date_v.day).zfill(2),str(h_cems).zfill(2), m)
						fichier = "/cnrm/phynh/data1/magnaldom/CEMS/DATA/%s_%s/%s/%s"  %(str(date_v.year), str(date_v.month).zfill(2), 'CMA', file )
						CEMS = Dataset(fichier, "r", format="NETCDF4")
						#print("Trouvé : " + fichier)

				if  date_v.month == 4 or date_v.month == 5 or date_v.month == 6 or date_v.month == 7 or date_v.month == 8 or date_v.month == 9 or date_v.month == 10 or date_v.month == 11 or date_v.month == 12 : #pas le même nom selon les mois
					file = 'S_NWC_%s_MSG4_globeM-VISIR_%s%s%sT%s%s00Z.nc'%("CMA", str(date_v.year), str(date_v.month).zfill(2), str(date_v.day).zfill(2),str(h_cems).zfill(2), m)
					fichier = "/cnrm/phynh/data1/magnaldom/CEMS/DATA/%s_%s/%s/%s"  %(str(date_v.year), str(date_v.month).zfill(2), 'CMA', file )
					CEMS = Dataset(fichier, "r", format="NETCDF4")
				if date_v.month == 2:
					file = 'S_NWC_%s_MSG4_globeM-VISIR_%s%s%sT%s%s00Z.nc.v2018'%("CMA", str(date_v.year), str(date_v.month).zfill(2), str(date_v.day).zfill(2),str(h_cems).zfill(2), m)
					fichier = "/cnrm/phynh/data1/magnaldom/CEMS/DATA/%s_%s/%s/%s"  %(str(date_v.year), str(date_v.month).zfill(2), 'CMA', file )
					CEMS = Dataset(fichier, "r", format="NETCDF4")


				for i in range(0,nb_station):
					cos_zen =  zenith_angle_ephem(lat[i], lon[i], datetime(annee, mois, jour,h)+timedelta(minutes = -30)) #Only when zenith angle is lower than a threshold
					if classe[i]<4 and num_station[i]!='6088001' and alti[i]<1000 and cos_zen > 0.1 : #Selection of the stations

						Z =  CEMS.variables['cma'][int(M_sat[i]),int(N_sat[i])]
						nebul_cems_h[i] = nebul_cems_h[i] + Z

						if m=='00': #on identifie pour AROME une seule fois par heure

							#---------------------------
							#-----------AROME-----------
							#---------------------------
							csd_a, nebul = dcc_AROME_SAT_cumul_zone(arome, datetime(annee, mois, jour, h), lon[i], lat[i], int(M[i]), int(N[i]), int(K), fn_lim)
							csd_arome_h[i] = csd_a

			except :
				print("Fichier non trouvé : S_NWC_%s_MSG4_globeM-VISIR_%s%s%sT%s%s00Z.nc" %("CMA", str(date_v.year), str(date_v.month).zfill(2), str(date_v.day).zfill(2),str(h).zfill(2), m))

		#--------------------------------------------
		#- DEBUT POSTPROCESSING (chaque heure)-------
		#--------------------------------------------
		for i in range(0,nb_station):
			if csd_arome_h[i]==0 or csd_arome_h[i]==1: #To take only valid station/cosSZA
				if nebul_cems_h[i]>0:
					csd_cems_h[i]=1
				else :
					csd_cems_h[i]=0
				
				#--------------------------------------------------------
				#---------Statistiques calculation--------------------------------
				#--------------------------------------------------------
				cos_zen =  zenith_angle_ephem(lat[i], lon[i], datetime(annee, mois, jour,h)+timedelta(minutes = -30))
				S0 = get_solar_constant(date_v)
				S = 1365*S0
				TOA_estimation = S*cos_zen
				flux_obs_norm = ghi_b/TOA_estimation
				#Recover information each hour and each station
				ghi_b = plot_BDClim_station_h(num_station[i],date_v, h) #inf si pas present
				csd_a = csd_arome_h[i] 
				csd_b = csd_cems_h[i]

				#Calculate GHI in AROME with the spatial tolerance method
				Z, lon_zone, lat_zone = zone_M_N_lon_lat_opt(arome, M[i],N[i], date_v, h,K)
				Z_tot = Z.flatten()
				ecart = [np.abs(Z_tot[bou]-(ghi_b) )for bou in range(len(Z_tot))]
				val_min = stats.scoreatpercentile(ecart, 10, interpolation_method="lower")
				position = ecart.index(val_min)
				ghi_zone = Z_tot[position]
				flux_mod_norm = ghi_zone/TOA_estimation


				if csd_a == 1 and csd_b== 1 and ghi_b<20000 and ghi_b !=0: #If clouds present in the mod and obs
					score[i][0] = score[i][0] + 1
					biais_norm[i][0] = biais_norm[i][0] + flux_mod_norm - flux_obs_norm
					biais[i][0] = biais[i][0] + ghi_zone - ghi_b
					MAPE[i][0] = MAPE[i][0] + (ghi_zone - ghi_b)/ghi_b
					MAE[i][0] = MAE[i][0] + np.abs(ghi_zone - ghi_b)
					if ghi_zone - ghi_b > 0:
						biais_pos[i][0] = biais_pos[i][0] + ghi_zone - ghi_b
						score_pos[i][0] = score_pos[i][0] + 1
					elif ghi_zone - ghi_b <0:
						biais_neg[i][0] = biais_neg[i][0] + ghi_zone - ghi_b
						score_neg[i][0] = score_neg[i][0] + 1
					flu_moy_norm[i][0] = flu_moy_norm[i][0] + flux_obs_norm
					flu_moy[i][0] = flu_moy[i][0] + ghi_b
					ecart_type[i][0] = ecart_type[i][0] + (ghi_zone - ghi_b)**2
					array_cas1 = np.vstack((array_cas1, [ghi_zone, flux_mod_norm, ghi_b,flux_obs_norm ]))
					distribution_1.append(ghi_zone - ghi_b)
					
				if csd_a == 1 and csd_b == 0 and ghi_b<20000 and ghi_b !=0: #If clouds present in the mod but not in the obs
					score[i][1] = score[i][1] + 1
					biais_norm[i][1] = biais_norm[i][1] + flux_mod_norm - flux_obs_norm
					MAE[i][1] = MAE[i][1] + np.abs(ghi_zone - ghi_b)
					MAPE[i][1] = MAPE[i][1] + (ghi_zone - ghi_b)/ghi_b
					biais[i][1] = biais[i][1] + ghi_zone - ghi_b
					if ghi_zone - ghi_b > 0:
						biais_pos[i][1] = biais_pos[i][1] + ghi_zone - ghi_b
						score_pos[i][1] = score_pos[i][1] + 1
					elif ghi_zone - ghi_b <0:
						biais_neg[i][1] = biais_neg[i][1] + ghi_zone - ghi_b
						score_neg[i][1] = score_neg[i][1] + 1
					flu_moy_norm[i][1] = flu_moy_norm[i][1] + flux_obs_norm
					flu_moy[i][1] = flu_moy[i][1] + ghi_b
					ecart_type[i][1] = ecart_type[i][1] + (ghi_zone - ghi_b)**2
					array_cas2 = np.vstack((array_cas2, [ghi_zone, flux_mod_norm, ghi_b,flux_obs_norm ]))
					distribution_2.append(ghi_zone - ghi_b)
					
				if csd_a == 0 and csd_b == 1 and ghi_b<20000 and ghi_b !=0: #If clouds present in the obs but not in the mod
					score[i][2] = score[i][2] + 1
					biais_norm[i][2] = biais_norm[i][2] + flux_mod_norm - flux_obs_norm
					MAE[i][2] = MAE[i][2] + np.abs(ghi_zone - ghi_b)
					MAPE[i][2] = MAPE[i][2] + (ghi_zone - ghi_b)/ghi_b
					biais[i][2] = biais[i][2] + ghi_zone - ghi_b
					if ghi_zone - ghi_b > 0:
						biais_pos[i][2] = biais_pos[i][2] + ghi_zone - ghi_b
						score_pos[i][2] = score_pos[i][2] + 1
					elif ghi_zone - ghi_b <0:
						biais_neg[i][2] = biais_neg[i][2] + ghi_zone - ghi_b
						score_neg[i][2] = score_neg[i][2] + 1
					flu_moy_norm[i][2] = flu_moy_norm[i][2] + flux_obs_norm
					flu_moy[i][2] = flu_moy[i][2] + ghi_b
					ecart_type[i][2] = ecart_type[i][2] + (ghi_zone - ghi_b)**2
					array_cas3 = np.vstack((array_cas3, [ghi_zone, flux_mod_norm, ghi_b,flux_obs_norm ]))
					distribution_3.append(ghi_zone - ghi_b)
					
				if csd_a == 0 and csd_b == 0 and ghi_b<20000 and ghi_b !=0: #If clouds not present in the mod and obs
					score[i][3] = score[i][3] + 1
					biais_norm[i][3] = biais_norm[i][3] + flux_mod_norm - flux_obs_norm
					MAE[i][3] = MAE[i][3] + np.abs(ghi_zone - ghi_b)
					MAPE[i][3] = MAPE[i][3] + (ghi_zone - ghi_b)/ghi_b
					biais[i][3] = biais[i][3] + ghi_zone - ghi_b
					if ghi_zone - ghi_b > 0:
						biais_pos[i][3] = biais_pos[i][3] + ghi_zone - ghi_b
						score_pos[i][3] = score_pos[i][3] + 1
					elif ghi_zone - ghi_b <0:
						biais_neg[i][3] = biais_neg[i][3] + ghi_zone - ghi_b
						score_neg[i][3] = score_neg[i][3] + 1
					flu_moy_norm[i][3] = flu_moy_norm[i][3] + flux_obs_norm
					flu_moy[i][3] = flu_moy[i][3] + ghi_b
					ecart_type[i][3] = ecart_type[i][3] + (ghi_zone - ghi_b)**2
					array_cas4 = np.vstack((array_cas4, [ghi_zone, flux_mod_norm, ghi_b,flux_obs_norm ]))
					distribution_4.append(ghi_zone - ghi_b)




	arome.close()


### ============================================================================
### SCORES CALCULATION


tc1 = 0
tc2 = 0
tc3 = 0
tc4 = 0

b1 = 0
b2 = 0
b3 = 0
b4 = 0
b1_pos = 0
b2_pos = 0
b3_pos = 0
b4_pos = 0
b1_neg = 0
b2_neg = 0
b3_neg = 0
b4_neg = 0
b1_norm = 0
b2_norm = 0
b3_norm = 0
b4_norm = 0

bt = 0
fm1 = 0
fm2 = 0
fm3 = 0
fm4 = 0

ft = 0
MAEt = 0
MAPEt = 0
ft_norm = 0
ett = 0

fm1_norm = 0
fm2_norm = 0
fm3_norm = 0
fm4_norm= 0


et1 = 0
et2 = 0
et3 = 0
et4 = 0



for i in range(nb_station):
	tc1 = tc1 + score[i][0]
	tc2 = tc2 + score[i][1]
	tc3 = tc3 + score[i][2]
	tc4 = tc4 + score[i][3]

	b1 = b1 + biais[i][0]
	b2 = b2 + biais[i][1]
	b3 = b3 + biais[i][2]
	b4 = b4 + biais[i][3]
	b1_pos = b1_pos + biais_pos[i][0]
	b2_pos = b2_pos + biais_pos[i][1]
	b3_pos = b3_pos + biais_pos[i][2]
	b4_pos = b4_pos + biais_pos[i][3]
	b1_neg = b1_neg + biais_neg[i][0]
	b2_neg = b2_neg + biais_neg[i][1]
	b3_neg = b3_neg + biais_neg[i][2]
	b4_neg = b4_neg + biais_neg[i][3]
	b1_norm = b1_norm + biais_norm[i][0]
	b2_norm = b2_norm + biais_norm[i][1]
	b3_norm = b3_norm + biais_norm[i][2]
	b4_norm = b4_norm + biais_norm[i][3]

	bt = bt + biais[i][0]+ biais[i][1] + biais[i][2] + biais[i][3]
	MAEt = MAEt + MAE[i][0]+ MAE[i][1] + MAE[i][2] + MAE[i][3]
	MAPEt = MAPEt + MAPE[i][0]+ MAPE[i][1] + MAPE[i][2] + MAPE[i][3]
	ft = ft + flu_moy[i][0]+ flu_moy[i][1] + flu_moy[i][2] + flu_moy[i][3]
	ett = ett + ecart_type[i][0]+ ecart_type[i][1] + ecart_type[i][2] + ecart_type[i][3]
	ft_norm = ft_norm + flu_moy_norm[i][0]+ flu_moy_norm[i][1] + flu_moy_norm[i][2] + flu_moy_norm[i][3]
	fm1 = fm1 + flu_moy[i][0]
	fm2 = fm2 + flu_moy[i][1]
	fm3 = fm3 + flu_moy[i][2]
	fm4 = fm4 + flu_moy[i][3]

	fm1_norm = fm1_norm + flu_moy_norm[i][0]
	fm2_norm = fm2_norm + flu_moy_norm[i][1]
	fm3_norm = fm3_norm + flu_moy_norm[i][2]
	fm4_norm = fm4_norm + flu_moy_norm[i][3]


	et1 = et1 + ecart_type[i][0]
	et2 = et2 + ecart_type[i][1]
	et3 = et3 + ecart_type[i][2]
	et4 = et4 + ecart_type[i][3]


Nh = tc1 + tc2 + tc3 + tc4

freq = [tc1/Nh, tc2/Nh, tc3/Nh, tc4/Nh]

biais = [b1/tc1, b2/tc2, b3/tc3, b4/tc4]


biais_pos = [b1_pos/tc1, b2_pos/tc2, b3_pos/tc3, b4_pos/tc4]


biais_neg = [b1_neg/tc1, b2_neg/tc2, b3_neg/tc3, b4_neg/tc4]


biais_norm = [b1_norm/tc1, b2_norm/tc2, b3_norm/tc3, b4_norm/tc4]

cont1 = b1/Nh
cont2 = b2/Nh
cont3 = b3/Nh
cont4 = b4/Nh

cont = [cont1, cont2, cont3, cont4]


ecart_type = [np.sqrt(et1/tc1), np.sqrt(et2/tc2), np.sqrt(et3/tc3), np.sqrt(et4/tc4)]


biais_rel = [(b1/tc1)/(fm1/tc1),(b2/tc2)/(fm2/tc2),(b3/tc3)/(fm3/tc3),(b4/tc4)/(fm4/tc4)]


biais_rel_norm = [(b1_norm/tc1)/(fm1_norm/tc1), (b2_norm/tc2)/(fm2_norm/tc2), (b3_norm/tc3)/(fm3_norm/tc3), (b4_norm/tc4)/(fm4_norm/tc4)]


ecart_type_rel = [np.sqrt(et1/tc1)/(fm1/tc1), np.sqrt(et2/tc2)/(fm2/tc2),np.sqrt(et3/tc3)/(fm3/tc3), np.sqrt(et4/tc4)/(fm4/tc4)]


flux_moyen = [fm1/tc1, fm2/tc2, fm3/tc3, fm4/tc4]


flux_moyen_norm = [fm1_norm/tc1, fm2_norm/tc2, fm3_norm/tc3, fm4_norm/tc4]


standard_deviation_1 = np.std(array_cas1, 0)
standard_deviation_2 = np.std(array_cas2, 0)
standard_deviation_3 = np.std(array_cas3, 0)
standard_deviation_4 = np.std(array_cas4, 0)

mat1 = np.corrcoef(array_cas1[:,0],array_cas1[:,2], rowvar = False)
mat2 = np.corrcoef(array_cas2[:,0],array_cas2[:,2], rowvar = False)
mat3 = np.corrcoef(array_cas3[:,0],array_cas3[:,2], rowvar = False)
mat4 = np.corrcoef(array_cas4[:,0],array_cas4[:,2], rowvar = False)

corrcoeff = [mat1[0][1],mat2[0][1],mat3[0][1],mat4[0][1]]

RMSE_norm = [np.sqrt(mean_squared_error(array_cas1[:,1],array_cas1[:,3])), np.sqrt(mean_squared_error(array_cas2[:,1],array_cas2[:,3])),np.sqrt(mean_squared_error(array_cas3[:,1],array_cas3[:,3])),np.sqrt(mean_squared_error(array_cas4[:,1],array_cas4[:,3]))]


all_obs =  np.empty((0,1), float)
all_obs = np.append(all_obs, array_cas1[:,2])
all_obs = np.append(all_obs, array_cas2[:,2])
all_obs = np.append(all_obs, array_cas3[:,2])
all_obs = np.append(all_obs, array_cas4[:,2])

STD_all_obs = np.std(all_obs, 0)

all_obs_norm =  np.empty((0,1), float)
all_obs_norm = np.append(all_obs_norm, array_cas1[:,3])

all_obs_norm = np.append(all_obs_norm, array_cas2[:,3])
all_obs_norm = np.append(all_obs_norm, array_cas3[:,3])
all_obs_norm = np.append(all_obs_norm, array_cas4[:,3])
STD_all_obs_norm = np.std(all_obs_norm, 0)


### ============================================================================
# ### WRITING IN A .NC

file_nc = nc.Dataset("RESULT/"+str(date_temp1)+"_"+str(date_temp2)+"_"+str(K)+"_"+str(fn_lim)+".nc","w",format="NETCDF4")
file_nc.createDimension('cas', 4)
file_nc.createDimension('aux', 1)
file_nc.createDimension('const', 1)
K_nc   = file_nc.createVariable('K', 'f4',dimensions=('aux'))
fn_lim_nc   = file_nc.createVariable('fn_lim', 'f4', dimensions=('aux'))
Nh_nc = file_nc.createVariable('Nh', 'f4', dimensions=('aux'))
biais_total = file_nc.createVariable('biais_total', 'f4', dimensions=('const'))
flux_total = file_nc.createVariable('flux_total', 'f4', dimensions=('const'))
MAE_total = file_nc.createVariable('MAE_total', 'f4', dimensions=('const'))
MAPE_total = file_nc.createVariable('MAPE_total', 'f4', dimensions=('const'))
RMSE_total = file_nc.createVariable('RMSE_total', 'f4', dimensions=('const'))
flux_norm_total = file_nc.createVariable('flux_total_norm', 'f4', dimensions=('const'))
STD = file_nc.createVariable('STD_all_obs', 'f4', dimensions=('const'))
STD[:] = STD_all_obs
STD_norm = file_nc.createVariable('STD_all_obs_norm', 'f4', dimensions=('const'))
STD_norm[:] = STD_all_obs_norm
K_nc[:] = K
fn_lim_nc[:] = fn_lim
Nh_nc[:] = Nh
biais_total[:] = bt/Nh
MAE_total[:] = MAEt/Nh
MAPE_total[:] = MAPEt/Nh
flux_total[:] = ft/Nh
flux_norm_total[:] = ft_norm/Nh
RMSE_total[:] = np.sqrt(ett/Nh)


frequence_nc = file_nc.createVariable('frequence', 'f4', dimensions=('cas'), zlib=True, complevel=1)
frequence_nc[:] = freq
biais_nc = file_nc.createVariable('biais', 'f4', dimensions=('cas'), zlib=True, complevel=1)
biais_nc[:] = biais
biais_pos_nc = file_nc.createVariable('biais_positif', 'f4', dimensions=('cas'), zlib=True, complevel=1)
biais_pos_nc[:] = biais_pos
biais_neg_nc = file_nc.createVariable('biais_negatif', 'f4', dimensions=('cas'), zlib=True, complevel=1)
biais_neg_nc[:] = biais_neg
biais_norm_nc = file_nc.createVariable('biais_norm', 'f4', dimensions=('cas'), zlib=True, complevel=1)
biais_norm_nc[:] = biais_norm
cont_nc = file_nc.createVariable('contribution_au_biais', 'f4', dimensions=('cas'), zlib=True, complevel=1)
cont_nc[:] = cont
ecart_type_nc = file_nc.createVariable('ecart_type', 'f4', dimensions=('cas'), zlib=True, complevel=1)
ecart_type_nc[:] = ecart_type
biais_rel_nc = file_nc.createVariable('biais_relatif', 'f4', dimensions=('cas'), zlib=True, complevel=1)
biais_rel_nc[:] = biais_rel
biais_rel_norm_nc = file_nc.createVariable('biais_relatif_norm', 'f4', dimensions=('cas'), zlib=True, complevel=1)
biais_rel_norm_nc[:] = biais_rel_norm
ecart_type_rel_nc = file_nc.createVariable('ecart_type_relatif', 'f4', dimensions=('cas'), zlib=True, complevel=1)
ecart_type_rel_nc[:] = ecart_type_rel

RMSE_norm_nc = file_nc.createVariable('RMSE_norm', 'f4', dimensions=('cas'), zlib=True, complevel=1)
RMSE_norm_nc[:] = RMSE_norm
flux_moyen_nc = file_nc.createVariable('flux_moyen', 'f4', dimensions=('cas'), zlib=True, complevel=1)
flux_moyen_nc[:] = flux_moyen
flux_moyen_norm_nc = file_nc.createVariable('flux_moyen_norm', 'f4', dimensions=('cas'), zlib=True, complevel=1)
flux_moyen_norm_nc[:] = flux_moyen_norm

#flux_moyen_mod_nc = file_nc.createVariable('flux_moyen_mod', 'f4', dimensions=('cas'), zlib=True, complevel=1)
#flux_moyen_mod_nc[:] =
#
#flux_moyen_mod_norm_nc = file_nc.createVariable('flux_moyen_mod_norm', 'f4', dimensions=('cas'), zlib=True, complevel=1)
#flux_moyen_mod_norm_nc[:] = standard_deviation_1

standard_deviation_mod_nc = file_nc.createVariable('standard_deviation_mod', 'f4', dimensions=('cas'), zlib=True, complevel=1)
standard_deviation_mod_nc[:] = [standard_deviation_1[0],standard_deviation_2[0],standard_deviation_3[0],standard_deviation_4[0]]

standard_deviation_mod_norm_nc = file_nc.createVariable('standard_deviation_mod_norm', 'f4', dimensions=('cas'), zlib=True, complevel=1)
standard_deviation_mod_norm_nc[:] = [standard_deviation_1[1],standard_deviation_2[1],standard_deviation_3[1],standard_deviation_4[1]]

standard_deviation_obs_nc = file_nc.createVariable('standard_deviation_obs', 'f4', dimensions=('cas'), zlib=True, complevel=1)
standard_deviation_obs_nc[:] = [standard_deviation_1[2],standard_deviation_2[2],standard_deviation_3[2],standard_deviation_4[0]]

standard_deviation_obs_norm_nc = file_nc.createVariable('standard_deviation_obs_norm', 'f4', dimensions=('cas'), zlib=True, complevel=1)
standard_deviation_obs_norm_nc[:] = [standard_deviation_1[3],standard_deviation_2[3],standard_deviation_3[3],standard_deviation_4[3]]

corrcoeff_nc = file_nc.createVariable('corrcoeff', 'f4', dimensions=('cas'), zlib=True, complevel=1)
corrcoeff_nc[:] = corrcoeff



file_nc.close()

