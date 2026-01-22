#dirList = ['C:\\Users\\jared\\Downloads\\rhoLBT500_6Percent','C:\\Users\\jared\\Downloads\\rhoLBT500_7Percent','C:\\Users\\jared\\Downloads\\rhoLBT500_8Percent',
#          'C:\\Users\\jared\\Downloads\\rhoLBT500_9Percent']
#DensityList = ['6Percent','7percent','8Percent','9Percent']
#DensityList = ['15Percent']

#dirList = ['C:\\Users\\jared\\Downloads\\15Percent']

import os
import numpy as np
import pandas as pd
import scipy
from scipy.integrate import cumulative_trapezoid


dirList  = ['MC_Simulations/Simulation_Results/DataFiles/']
OutDir = 'Free_Energies/'


FreeEnergyList = list()
AnglesList = list()
for i in range(len(dirList)):
    dir1=dirList[i]
    files = os.listdir(dir1)

    thetas = list()
    meanVeff = list()
    atMin = list()
    forces = list()
    fileName = list()
    free_Energy=list()
    mean_Sep = list()
    distanceR = list()
    min_Sep = list()
    for file in files:
        filePath = os.path.join(dir1,file)
        df = pd.read_csv(filePath, skiprows=39, delimiter='\\s+')
        substring = 'theta0_'
        index = file.find(substring)
        str1=str(pd.read_csv(filePath).loc[3])
        #biaStr=float(str1[38:40])
        biaStr=20
        my_string=pd.read_csv(filePath).iloc[4,0]
        theta=float(my_string.split("   Bias location:   ",1)[1])
        thetas.append(theta)
        df_run = df[df['#Sweep(1)']>2E5]
        fBias = -2*biaStr*(np.mean(df_run['ThetaValTop(16)']-float(theta)))
        meanSep = np.mean(df_run['ThetaValTop(16)'])
        mean_Sep.append(meanSep)
        forces.append(fBias)
        meanVeff.append(np.mean(df_run['Veff(8)']))
    thetas2=thetas
    res=pd.DataFrame({'Angle':thetas2, 'atMin':forces,'MeanVeff':meanVeff,'MeanThetaTop':mean_Sep})
    res['Angle'] = res['Angle'].astype(int)
    res.sort_values(by='Angle')
    res2 = res[res['Angle']>15]
    
    newDF=res2.sort_values(by = ['Angle'])
    newDF.reset_index(drop=True, inplace=True)
    new_order = ['Angle', 'atMin', 'MeanVeff','MeanThetaTop']
    df_reordered = newDF[new_order]
    Mean_FE=df_reordered
    
    dfSlope=df_reordered[:len(df_reordered)-1]
    dfSlope=dfSlope[dfSlope['MeanVeff'] == 0]
    if(len(dfSlope > 0)):
        r_bg = dfSlope['MeanThetaTop'].values
        f_bg = dfSlope['atMin'].values
        a,b=np.polyfit(r_bg, f_bg, deg=1)
    else:
        a = 1
        b = 0
    newDF['Force_corr'] = newDF['atMin'] - (a * newDF['MeanThetaTop'] + b)
    newDF =newDF[:len(newDF)-1]
    ForceatPlateau = np.mean(newDF[newDF['MeanVeff'] == 0]['Force_corr'])
    newDF['Force_corr'] = newDF['Force_corr'] - ForceatPlateau
    print(newDF)
    FreeEnergy6=cumulative_trapezoid(newDF['Force_corr'],x=newDF['MeanThetaTop'], initial=0.0)
    F6 = FreeEnergy6 - FreeEnergy6[-1]

    
    
    FreeEnergyList.append(F6)
    AnglesList.append(newDF['MeanThetaTop'])
    meanFE=pd.DataFrame({'Distance':newDF['MeanThetaTop'],'Force':newDF['atMin'],'FE':F6})
    #print(meanFE['FE'])
    filename = "Free_Energy.dat"
    meanFE.to_csv(OutDir+ filename,sep = ' ', index=False)
