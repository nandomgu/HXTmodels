##HXT model selection. last update 20170315
from numpy import linspace, array
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from math import exp
from numpy import matrix,column_stack, where, ones, zeros, array, size, concatenate, nan
import numpy as np

####getting the experiment Data. experimentData is a dictionary of the form
#####experimentData[EXPT_DATE]		['time']	the timepoints of the experiment in hrs
#					['input']	the normalised input 
#					['strains']	the name of each of the strains in the experiment.
#					['means']	the mean expression data for each of the 'strains', by column

from AllExperimentData import experimentData, getStrains, getStrain
###ACCESSORY FUNCTIONS BLOCK.
##if in list return index
def inList(a, lst):
	if a in lst:
		return(array(where([a==j for j in fitStrains])).flatten())
	else:
		return nan

def removeInfinites(a):
	return(a[isfinite(a)])
  
## this generates input interpolation functions for any experiment in experimentData,
## such that the input can be retrieved at any timempoint
def inputFunction(date):
	return interp1d(experimentData[date]['time'],experimentData[date]['input']) 



#Function to create simulated step Inputs.
def step(len1, len2, len3):
	 return concatenate([zeros(len1), ones(len2), zeros(len3)])
#example: generating a step input with 180 min lag time, 480 min step length, and 500 min post step
stepInput=step(180, 480, 500)

tstep= array([j for j in range(0, size(stepInput))])
f= interp1d( tstep, stepInput)

###paramInf is a dictionary with details about each *possible* parameter in the model. different models might need access
###to the same or different parameter Information, and some parameters could potentially be repurposed.

paramInf={};
paramInf['VStd1']			={'defaultValue': 1, 'priorType': 'uniformExp', 'priorRange': [-6, 6]}; #%maximum production rate of Std1
paramInf['KSelfStd1']		= {'defaultValue': .4, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#threshold at which Std1 represses itself.
paramInf['KDegStd1']		= {'defaultValue': 1, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};# maximum glucose induced degradation rate of Std1
paramInf['threshStd1'] 	    = {'defaultValue': .8, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#threshold at which Std1 is degraded by glucose
paramInf['basalDegStd1'] 	= {'defaultValue': .05, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#basal degradation rate for Std1
paramInf['VMth1'] 			= {'defaultValue': 1, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#Maximum production rate for Mth1
paramInf['basalDegMth1']	= {'defaultValue': .05, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#basal degradation rate for Mth1
paramInf['KDegMth1']		= {'defaultValue': 1, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#maximum glucose induced deg. rate for Mth1
paramInf['threshMth1']		= {'defaultValue': .8, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};# threshold for Mth1 to be degraded by glucose
paramInf['KRepMth1']		= {'defaultValue': .01, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#Threshold for Mth1 to repress Std1
paramInf['thresholdStd1Mig1']= {'defaultValue': .1, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#threshold for std to inhibit mig1/ raise the activation threshold of mig1
paramInf['hillactiveMig1']	= {'defaultValue': 4, 'priorType': 'discrete', 'priorValues': [1, 4]};#hill coefficient for mig1 to become active upon glucose
paramInf['KSnf1']			= {'defaultValue': .01, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#threshold at which glucose inhibits snf1
paramInf['KRepSnf1Mig1']	= {'defaultValue': .9, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#threshold of snf1 needed to inhibit Mig1
paramInf['KEnhance']		= {'defaultValue': .1, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#threshold at which Std1 enhances Snf1

paramInf['KStd1HXT4']		= {'defaultValue': .2, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};		#	.2,[-6,6]};
paramInf['KMth1HXT1']		= {'defaultValue': .005, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#threshold at which Mth1 Represses Mth1
paramInf['KMth1HXT2']		={'defaultValue': 10, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};
paramInf['KMth1HXT3']		={'defaultValue': 10, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};
paramInf['KMth1HXT4']		= {'defaultValue': .2, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};		#	.2,[-6,6]};

paramInf['KMig1HXT2']		= {'defaultValue': 2, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#.005,[-6,6]};
paramInf['KMig1HXT4']		= {'defaultValue': 2, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#.005,[-6,6]};
paramInf['KMig1HXT7']		= {'defaultValue': 2, 'priorType': 'uniformExp', 'priorRange': [-6, 6]};#.005,[-6,6]};

paramInf['KDegSnf1']		={'defaultValue': 0.1, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['basalDegSnf1']	={'defaultValue': 0.1, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHxt1']		={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHxt2']		={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHxt3']		={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHxt4']		={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHxt7']		={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};

paramInf['KDegHxt1']		={'defaultValue': .05, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KDegHxt2']		={'defaultValue': .05, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KDegHxt3']		={'defaultValue': .05, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KDegHxt4']		={'defaultValue': .05, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KDegHxt7']		={'defaultValue': .05, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};

paramInf['VSnf1']		={'defaultValue': .05, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KFBSnf1']		={'defaultValue': .01, 'priorType': 'uniformExp', 'priorRange': [-6, -1]};
paramInf['VHXT1']		={'defaultValue': .6, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHXT2']		={'defaultValue': .6, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHXT3']		={'defaultValue': .6, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHXT4']		={'defaultValue': .6, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['VHXT7']		={'defaultValue': .6, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};

paramInf['KDegHXT1'] ={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KDegHXT2'] ={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KDegHXT3'] ={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KDegHXT4'] ={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};
paramInf['KDegHXT7'] ={'defaultValue': .3, 'priorType': 'uniformExp', 'priorRange': [-6, 0]};

## From the information in paramInf, Here we assign starting value for each of the parameters
## either the default from above or sampling from the distribution/range. right now we go for the default value.


params={};
params['VStd1']					=paramInf['VStd1']['defaultValue']
params['KSelfStd1']				=paramInf['KSelfStd1']['defaultValue']
params['KDegStd1']				=paramInf['KDegStd1']['defaultValue']
params['threshStd1']				=paramInf['threshStd1']['defaultValue'] 	 
params['basalDegStd1']				=paramInf['basalDegStd1']['defaultValue'] 
params['VMth1']					=paramInf['VMth1']['defaultValue'] 	
params['basalDegMth1']				=paramInf['basalDegMth1']['defaultValue'] 
params['KDegMth1']				=paramInf['KDegMth1']['defaultValue'] 	
params['threshMth1']				=paramInf['basalDegStd1']['defaultValue'] 	
params['KRepMth1']				=paramInf['KRepMth1']['defaultValue'] 
params['thresholdStd1Mig1']			=paramInf['thresholdStd1Mig1']['defaultValue'] 
params['hillactiveMig1']			=paramInf['hillactiveMig1']['defaultValue'] 
params['KSnf1']					=paramInf['KSnf1']['defaultValue'] 		
params['KRepSnf1Mig1']				=paramInf['KRepSnf1Mig1']['defaultValue'] 
params['KEnhance']				=paramInf['KEnhance']['defaultValue'] 

params['KStd1HXT4']				=paramInf['KStd1HXT4']['defaultValue'] 	
params['KMth1HXT1']				=paramInf['KMth1HXT1']['defaultValue'] 	
params['KMth1HXT2']				=paramInf['KMth1HXT2']['defaultValue'] 	
params['KMth1HXT3']				=paramInf['KMth1HXT3']['defaultValue'] 	
params['KMth1HXT4']				=paramInf['KMth1HXT4']['defaultValue'] 	

params['KMig1HXT2']				=paramInf['KMig1HXT2']['defaultValue'] 
params['KMig1HXT4']				=paramInf['KMig1HXT4']['defaultValue'] 	
params['KMig1HXT7']				=paramInf['KMig1HXT7']['defaultValue'] 	

params['KDegSnf1']				=paramInf['KDegSnf1']['defaultValue'] 	
params['basalDegSnf1']				=paramInf['basalDegSnf1']['defaultValue']

params['VHxt1']					=paramInf['VHxt1']['defaultValue']
params['VHxt2']					=paramInf['VHxt2']['defaultValue']
params['VHxt3']					=paramInf['VHxt3']['defaultValue']
params['VHxt4']					=paramInf['VHxt4']['defaultValue']
params['VHxt7']					=paramInf['VHxt7']['defaultValue']

params['KDegHxt1']				=paramInf['KDegHxt1']['defaultValue']
params['KDegHxt2']				=paramInf['KDegHxt2']['defaultValue']
params['KDegHxt3']				=paramInf['KDegHxt3']['defaultValue']
params['KDegHxt4']				=paramInf['KDegHxt4']['defaultValue']
params['KDegHxt7']				=paramInf['KDegHxt7']['defaultValue']

params['KDegHXT1']				=paramInf['KDegHXT1']['defaultValue']
params['KDegHXT2']				=paramInf['KDegHXT2']['defaultValue']
params['KDegHXT3']				=paramInf['KDegHXT3']['defaultValue']
params['KDegHXT4']				=paramInf['KDegHXT4']['defaultValue']
params['KDegHXT7']				=paramInf['KDegHXT7']['defaultValue']

params['VSnf1']					=paramInf['VSnf1']['defaultValue']
params['KFBSnf1']				=paramInf['KFBSnf1']['defaultValue']
params['VHXT1']					=paramInf['VHXT1']['defaultValue']
params['VHXT2']					=paramInf['VHXT2']['defaultValue']
params['VHXT3']					=paramInf['VHXT3']['defaultValue']
params['VHXT4']					=paramInf['VHXT4']['defaultValue']
params['VHXT7']					=paramInf['VHXT7']['defaultValue']

#Declaring the order of the variables each equation represents.

varnames=['Std1',
'Mth1',
'Snf1',
'Mig1',
'HXT1',
'HXT2',
'HXT3',
'HXT4',
'HXT7',
'Hxt1',
'Hxt2',
'Hxt3',
'Hxt4',
'Hxt7']

### declaring a default initial condition for every variable.
initialConditionsDict={
'Std1':0,
'Mth1':0,
'Snf1':0,
'Mig1':0,
'HXT1':0,
'HXT2':0,
'HXT3':0,
'HXT4':0,
'HXT7':0,
'Hxt1':0,
'Hxt2':0,
'Hxt3':0,
'Hxt4':0,
'Hxt7':0,
}

### we create an array from the dictionary above (the model probably requires an array and dicts are disordered)
initialConditions= array([initialConditionsDict[j] for j in varnames])

###Here the parameters have individual variable names. Not currently in use, and not all parameters are represented.
VStd1=params['VStd1'];
KSelfStd1=params['KSelfStd1'];
KDegStd1=params['KDegStd1'];
threshStd1=params['threshStd1'];
basalDegStd1=params['basalDegStd1'];
VMth1= params['VMth1'];
basalDegMth1=params['basalDegMth1'];
KDegMth1=params['KDegMth1'];
threshMth1=params['threshMth1'];
KRepMth1=params['KRepMth1'];
KSnf1= params['KSnf1'];
KRepSnf1Mig1=params['KRepSnf1Mig1'];
KEnhance=params['KEnhance'];
hillactiveMig1=params['hillactiveMig1'];
KStd1HXT4=params['KStd1HXT4'];

KMth1HXT1=params['KMth1HXT1'];
KMth1HXT2=params['KMth1HXT2'];
KMth1HXT3=params['KMth1HXT3'];
KMth1HXT4=params['KMth1HXT4'];

KMig1HXT2=params['KMig1HXT2'];
KMig1HXT4=params['KMig1HXT4'];
KMig1HXT7=params['KMig1HXT7'];
KMig1HXT4=params['KMig1HXT4'];

KDegSnf1=params['KDegSnf1'];
basalDegSnf1=params['basalDegSnf1'];

#### assigning the starting values to a parameter array that will be given to the model. we explicitly write the numbers
#### for reference

prm=zeros(47);
prm[0]	=params['VStd1']		
prm[1]	=params['KSelfStd1']	
prm[2]	=params['KDegStd1']	
prm[3]	=params['threshStd1']
prm[4]	=params['basalDegStd1']
prm[5]	=params['VMth1']		
prm[6]	=params['basalDegMth1']
prm[7]	=params['KDegMth1']	
prm[8]	=params['threshMth1']
prm[9]	=params['KRepMth1']	
prm[10]	=params['thresholdStd1Mig1']
prm[11]	=params['hillactiveMig1']
prm[12]	=params['KSnf1']		
prm[13]	=params['KRepSnf1Mig1']
prm[14]	=params['KEnhance']	
prm[15]	=params['KStd1HXT4']	
prm[16]	=params['KMth1HXT1']	
prm[17]	=params['KMth1HXT2']	
prm[18]	=params['KMth1HXT3']	
prm[19]	=params['KMth1HXT4']	
prm[20]	=params['KMig1HXT2']	
prm[21]	=params['KMig1HXT4']	
prm[22]	=params['KMig1HXT7']	
prm[23]	=params['KDegSnf1']	
prm[24]	=params['basalDegSnf1']
prm[25]	=params['VHxt1']		
prm[26]	=params['VHxt2']		
prm[27]	=params['VHxt3']		
prm[28]	=params['VHxt4']		
prm[29]	=params['VHxt7']		
prm[30]	=params['KDegHxt1']	
prm[31]	=params['KDegHxt2']	
prm[32]	=params['KDegHxt3']	
prm[33]	=params['KDegHxt4']	
prm[34]	=params['KDegHxt7']	
prm[35]	=params['KDegHXT1']	
prm[36]	=params['KDegHXT2']	
prm[37]	=params['KDegHXT3']	
prm[38]	=params['KDegHXT4']	
prm[39]	=params['KDegHXT7']	
prm[40]	=params['VSnf1']		
prm[41]	=params['KFBSnf1']	
prm[42]	=params['VHXT1']		
prm[43]	=params['VHXT2']		
prm[44]	=params['VHXT3']		
prm[45]	=params['VHXT4']		
prm[46]	=params['VHXT7']		




##generate model function for HXT model 1 for any input function I. 
##Here we assume that the function environment has access to the params variable.
#the output is a model function where the glucose can be different.

def inputModel1(I):
	def func(incons, t=0):
		return array(
		[params['VStd1']/(1+(incons[0]/params['KSelfStd1']) +(incons[1]/params['KRepMth1'])) - params['basalDegStd1']*incons[0] -params['KDegStd1']*I(t)*incons[0], #incons[0]=Std1 #DStd1
		params['VMth1']- params['basalDegMth1']*incons[1]- ((params['KDegMth1']*I(t))**4/(params['threshMth1']**4+I(t)**4))*incons[1],  #incons[1]= Mth1 #DMth1
		params['VSnf1']+params['KFBSnf1']*incons[2]+I(t)/(params['KDegSnf1']+I(t))  - params['basalDegSnf1']*incons[2], ##Snf1 equation incons[2]= Snf1
		I(t)**4/(0.1**4+I(t)**4)  -.6*incons[3], ##Mig1 equation. incons[3] = Mig1
		params['VHXT1']/(1+ incons[1]/.0001)- params['KDegHXT1']*incons[4], ##HXT1 equation. incons[4]= HXT1. regulated by Mth1 
		params['VHXT2']/(1+((incons[1]/params['KMth1HXT2'])**4)+(incons[3]/params['KMig1HXT2']))-params['KDegHXT2']*incons[5], ##HXT2 equation.##incons[5]= HXT2
		params['VHXT3']/(1+ incons[1]/params['KMth1HXT3']**4)- params['KDegHXT3']*incons[6],##HXT3 equation.#incons[6]= HXT3
		params['VHXT4']/(1+(incons[1]/params['KMth1HXT4'])**4+(incons[3]/params['KMig1HXT4'])**4)- params['KDegHXT4']*incons[7], ##HXT4 equation. incons[7]=HXT4
		params['VHXT7']/(1+ incons[3]/params['KMig1HXT7']**4)- params['KDegHXT7']*incons[8],##HXT7 equation. incons[8]=HXT7
		params['VHxt1']*incons[4]- params['KDegHxt1']*incons[4],##Hxt1 equation.
		params['VHxt2']*incons[5]- params['KDegHxt2']*incons[5],##Hxt2 equation.
		params['VHxt3']*incons[6]- params['KDegHxt3']*incons[6],##Hxt3 equation.
		params['VHxt4']*incons[7]- params['KDegHxt4']*incons[7],##Hxt4 equation.
		params['VHxt7']*incons[8]- params['KDegHxt7']*incons[8]	#Hxt7 equation
		])
	return func

#prm is the parameter array. expd is the experimentData structure. we pass the values from prm to internal structure params.
#fitStrains specifies the strains we are interested in fitting. 
#outputs:
#	totalerrs= sum of all errors from all gene specific fittings. not normalised in any way
#	lsqerrs= dictionary where the fitting error is stored by experiment and strain.
#	ints= full results of the simulation of every experiment, for plotting purposes

def model0(prm,experimentData, fitStrains, plotSim=0): 
	params['VStd1']			=prm[0]	
	params['KSelfStd1']		=prm[1]	
	params['KDegStd1']		=prm[2]	
	params['threshStd1']		=prm[3]	
	params['basalDegStd1']		=prm[4]	
	params['VMth1']			=prm[5]	
	params['basalDegMth1']		=prm[6]	
	params['KDegMth1']		=prm[7]	
	params['threshMth1']		=prm[8]	
	params['KRepMth1']		=prm[9]	
	params['thresholdStd1Mig']	=prm[10]	
	params['hillactive']		=prm[11]	
	params['KSnf1']			=prm[12]	
	params['KRepSnf1Mig1']		=prm[13]	
	params['KEnhance']		=prm[14]	
	params['KStd1HXT4']		=prm[15]	
	params['KMth1HXT1']		=prm[16]	
	params['KMth1HXT2']		=prm[17]	
	params['KMth1HXT3']		=prm[18]	
	params['KMth1HXT4']		=prm[19]	
	params['KMig1HXT2']		=prm[20]	
	params['KMig1HXT4']		=prm[21]	
	params['KMig1HXT7']		=prm[22]	
	params['KDegSnf1']		=prm[23]	
	params['basalDegSnf1']		=prm[24]	
	params['VHxt1']			=prm[25]	
	params['VHxt2']			=prm[26]	
	params['VHxt3']			=prm[27]	
	params['VHxt4']			=prm[28]	
	params['VHxt7']			=prm[29]	
	params['KDegHxt1']		=prm[30]	
	params['KDegHxt2']		=prm[31]	
	params['KDegHxt3']		=prm[32]	
	params['KDegHxt4']		=prm[33]	
	params['KDegHxt7']		=prm[34]	
	params['KDegHXT1']		=prm[35]	
	params['KDegHXT2']		=prm[36]	
	params['KDegHXT3']		=prm[37]	
	params['KDegHXT4']		=prm[38]	
	params['KDegHXT7']		=prm[39]	
	params['VSnf1']			=prm[40]	
	params['KFBSnf1']		=prm[41]	
	params['VHXT1']			=prm[42]	
	params['VHXT2']			=prm[43]	
	params['VHXT3']			=prm[44]	
	params['VHXT4']			=prm[45]	
	params['VHXT7']			=prm[46]	

	#Define initial conditions, default is all zeros
	totalerr=0
	ints={} ###ints is a list that contains the numerical integration results
	lsqerrs={} ### simulation error dictionary
	####for each of the experiments, simulate model1
	##it seems you are not allowed to interpolate at the very edge of the experiment time so we go one point before.
	##we find what variables in the model correspond to the strains we are interested in fitting
	fitStrainIndices=array([where([j==strain for j in varnames]) for strain in fitStrains]).flatten()
	##we build a dictionary to access these indices by name.
	simStrainIndices=dict(zip( fitStrains, fitStrainIndices))
	j=0;
	for date in experimentData.keys():
		#for every experiment we simulate all variables
		ints[date]=integrate.odeint(inputModel1(inputFunction(date)),initialConditions,experimentData[date]['time'][0:-1]) 
		##indices of strains to be fitted
		##
		dataDict=getStrains(date, fitStrains)
		simLSQ={}
		for strn in dataDict.keys():
			if size(dataDict[strn])==0:
				simLSQ[strn]=nan
			else:
				simLSQ[strn]= sum((ints[date][:, simStrainIndices[strn]]-dataDict[strn][0:-1])**2)
				totalerr+=sum((ints[date][:, simStrainIndices[strn]]-dataDict[strn][0:-1])**2)
		if plotSim==1:
			plt.plot(experimentData[date]['time'][0:-1], ints[date][:, [simStrainIndices[j] for j in fitStrains]])
		lsqerrs[date]=simLSQ	

	return totalerr, lsqerrs, ints
