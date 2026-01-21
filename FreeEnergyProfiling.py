#dirList = ['C:\\Users\\jared\\Downloads\\rhoLBT500_6Percent','C:\\Users\\jared\\Downloads\\rhoLBT500_7Percent','C:\\Users\\jared\\Downloads\\rhoLBT500_8Percent',
#          'C:\\Users\\jared\\Downloads\\rhoLBT500_9Percent']
#DensityList = ['6Percent','7percent','8Percent','9Percent']
#DensityList = ['15Percent']

#dirList = ['C:\\Users\\jared\\Downloads\\15Percent']

dirList  = ['rhoLBT500/']
OutDir = 'Free_Energies/'


FreeEnergyList = list()
AnglesList = list()
densityList = list()
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
        biaStr=float(str1[38:40])
        #theta = file[index + len(substring):index + len(substring)+2]
        my_string=pd.read_csv(filePath).iloc[4,0]
        theta=float(my_string.split("   Bias location:   ",1)[1])
        thetas.append(theta)
        df_run = df[df['#Sweep(1)']>2E5]
        #BiasStrength = 50
        #avgSep = np.mean(df['ZTop(10)']-df['ZBot(9)'])
        #fBias = -2*BiasStrength*(np.mean(df_run['ThetaValTop(16)']))
        fBias = -2*biaStr*(np.mean(df_run['ThetaValTop(16)']-float(theta)))
        meanSep = np.mean(df_run['ThetaValTop(16)'])
        mean_Sep.append(meanSep)
        forces.append(fBias)
        meanVeff.append(np.mean(df_run['Veff(8)']))
        print(len(df))
    #thetas2 = [s.replace(".", "") for s in thetas]
    thetas2=thetas
    res=pd.DataFrame({'Angle':thetas2, 'atMin':forces,'MeanVeff':meanVeff,'MeanThetaTop':mean_Sep})
    res['Angle'] = res['Angle'].astype(int)
    res.sort_values(by='Angle')
    res2 = res[res['Angle']>15]
    
    newDF=res2.sort_values(by = ['Angle'])
    newDF.reset_index(drop=True, inplace=True)
    new_order = ['Angle', 'atMin', 'MeanVeff','MeanThetaTop','MinThetaTop']
    df_reordered = newDF[new_order]
    Mean_FE=df_reordered
    
    dfSlope=df_reordered[:len(df_reordered)-1]
    dfSlope=dfSlope[dfSlope['MeanVeff'] == 0]
    #r_bg = dfSlope['Angle'].values
    r_bg = dfSlope['MeanThetaTop'].values
    f_bg = dfSlope['atMin'].values
    a,b=np.polyfit(r_bg, f_bg, deg=1)
    #newDF['Force_corr'] = newDF['atMin'] - (a * newDF['Angle'] + b)
    newDF['Force_corr'] = newDF['atMin'] - (a * newDF['MeanThetaTop'] + b)
    newDF =newDF[:len(newDF)-1]
    ForceatPlateau = np.mean(newDF[newDF['MeanVeff'] == 0]['Force_corr'])
    newDF['Force_corr'] = newDF['Force_corr'] - ForceatPlateau

    FreeEnergy6=cumulative_trapezoid(newDF['Force_corr'],x=newDF['MeanThetaTop'], initial=0.0)
    F6 = FreeEnergy6 - FreeEnergy6[-1]

    
    
    FreeEnergyList.append(F6)
    AnglesList.append(newDF['MeanThetaTop'])
    density = [DensityList[i]] * len(F6)
    densityList.append(density)
    meanFE=pd.DataFrame({'Distance':newDF['MeanThetaTop'],'Force':newDF['atMin'],'FE':F6})
    filename = f"{DensityList[i]}.dat"
    meanFE.to_csv(OutDir+ filename,sep = ' ', index=False)
