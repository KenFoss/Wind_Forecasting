import types
from numpy.lib import nanfunctions
import scipy.stats as s
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from scipy.special import gamma
from tabulate import tabulate
import scipy.integrate as integrate

#This Model takes a csv files of wind speeds in the first column, in m/s, and calculates a Weibull distribution with R^2 and RMSE Error. 
#Future additions include other methods of calculations for the k and c parameters, and other models such as the Rayleigh Distribution


#change wroking directory to find csv
#When debugging eliminate path and os.chdir(path), comment them out, it will find the file. It does not when run normally.
#The above applies to running this code in the command prompt as well
#This code will now run without the following code in any aspect... but it could be useful and will be left
import os
#path = ('Desktop\Python\WeibullOffline')  #set path
#os.chdir(path)
print("Current working directory: {0}".format(os.getcwd()))

#List for speeds to be read into from CSV
Tv = []

#Values needed for Maximum Likelihood method
vzero = 0

#For Probabilities of wind speeds of zero 
prob = 0

#For Error Calculations, generic between methods
yi = [] #frequencies of observed wind speeds at each step


n = 0

#Adding P(v>0) as P, probability wind speed exceeds zero, k is the shape factor and c is scale factor. All calculated on excel.
#potentially will add in incremental wind speed
def weib(p, v, c, k):
    return p * (k/c) * (v/c) ** (k-1) * np.exp(-(v/c)**k)

sendesign = 'Collins_Wind_Data.csv'
research = 'Nebraska_Data.csv'

with open(sendesign, 'r') as read_obj:
    reader = csv.reader(read_obj)
    for row in reader:
        n += 1;
        speed = float(row[0]); #We read speeds as dictionary with header, row[0] gets value as string, convert to float 
        Tv.append(speed);
        if speed > 0:
            vzero +=1
        prob = float(vzero/n) #probability wind speed exceeds zero

            
print("iteration complete for data import")
        


#Weibull distribution over given range for error calculation
kmlm = 0
cmlm = 0
#for testing only
#kmlm = 2.486088774563999 #Test
#cmlm = 3.6310269503718406 #Test

#Calculate the shape parameter, Iterative
def kmlmcalc():
    krange = 2 #Guess for K value
    i = 0
    lnv = []
    vlnv = []
    vtok = []
    v = []
    while i < len(Tv):
        speed = Tv[i]
        if speed == 0: #zero wind speed makes log indefinite
            i += 1
        else:
            # I would like to see if using a krange calculated for all the data fixes my error issues
            if krange < 0: 
                print("i = {}, lnv = {}, vlnv = {}, vtok = {}".format(i, lnv[i-1], vlnv[i-1], vtok[i-1]))
            lnv.append(np.log(speed))
            vtok.append(speed ** krange)  #current wind speed to power of k, for scale factor calc
            vlnv.append((speed ** krange) * (np.log(speed)))      
            v.append(speed)
            if i >= 1:
                krange = (((sum(vlnv) / sum(vtok)) - (sum(lnv) / n)) ** (-1)) #This code will find the k factor, must be solved iteratively with an original k value guess
            #print("speed = {}, krange = {}. sumvlv = {}, sumvtok = {}, sumlnv = {}".format(speed, krange, sum(vlnv), sum(vtok), sum(lnv)))  #TESTING FUNCTION
            if math.isnan(krange): #why is krange coming out as nan, TESTING FUNCTION, problem solved n was not implemented
                print("krange is nan at iteration {}. sumvlv = {}, sumvtok = {}, sumlnv = {}".format(krange, sum(vlnv), sum(vtok), sum(lnv)))
            
            i += 1

    cv = []
    j = 0
    while j < len(Tv):
        speed = Tv[j]
        if speed == 0: #zero wind speed makes log indefinite
            j += 1
        else:
            cv.append(speed**krange)
        j += 1
    crange = ((1/n) * sum(cv)) ** (1/krange)  #This code will find c scale factor, can be solved explicitly
    #Assign global values to be used elsewhere
    global kmlm
    kmlm = krange
    global cmlm
    cmlm = crange

#Calculate k and c for maximum likelihood method
kmlmcalc()

def weibull(v): 
    return prob * (kmlm/cmlm) * ((v/cmlm) ** (kmlm - 1)) * np.exp(-(v/cmlm)**kmlm)


def weibmlm(x, y, use): #Weibull Within Time Range for error calc
    v = []

    points = 100
    errorv = np.linspace(x, y, points)

    #if ner == 1:
        #With one value in the bin it will scew data, better to return zero
    #    return 0

    if use == 'error':
        #From probability textbook, probability that X takes a value in interval [a,b] is the integral from a to b of the function
        result = integrate.quad(weibull, x, y)
        return result[0]
       
        
    elif use == 'plot':
        print("iteration complete for plotting weibull, maximum likelihood method")
        xv = np.linspace(min(Tv), max(Tv), 1000)
        return prob * (kmlm/cmlm) * ((xv/cmlm) ** (kmlm-1)) * np.exp(-(xv/cmlm)**kmlm)

#The following function models a weibull distribution by moment method estimation of the c and k parameters
kmm = int #k value for moment method
cmm = int #c value for moment method
def weibmm(x, y, use):
    v = []
    i = 0
    ner = 0  #n starts at zero, but every entry to v adds 1, n= true n
    while i < len(Tv):
        speed = Tv[i]
        i += 1
        if x <= speed <= y:
            ner += 1 #Whend doing error calculations, there will be less then the total entries being evaluated, this will account for that in mean calculation
            v.append(speed)
    if len(v) <= 1:
        return 0
    vmean = sum(v)/(n)
    stdadd = 0
    j = 0  
    while j < len(v): #sum vi-vmean for each wind speed, to be used in k value calculation
        stdadd += ((v[j] - vmean) ** 2) 
        j += 1
    if n == 0 or n == 1:
        print("ALL STOP")
    stdev = ((1 / n) * stdadd) ** (0.5) 
    k = (stdev/vmean) ** (-1.086)
    c = vmean/(gamma(1 + (1/k)))
    if use == 'error':
        #print("Iteration complete weibull moment method, error calculation for range {} - {}".format(x, y))
        if len(v) != 0 and c != 0:
            weiber = 0
            for i in v:
                weiber += prob * (k/c) * (i/c) ** (k-1) * np.exp(-(i/c)**k) #weibull distribution, must add probability for all values of v in this range
            return weiber/ner
        else:
            return 0
    elif use == 'plot':
        print("iteration complete for plotting weibull, moment method")
        global kmm #Set global kmlm, to change the value of kmlm outside of this function...same for cmlm
        kmm = k #Set global variables for later display
        global cmm 
        cmm = c #''
        xv = np.linspace(min(Tv), max(Tv), 1000)
        return prob * (k/c) * (xv/c) ** (k-1) * np.exp(-(xv/c)**k)
   

#Gets frequencies of wind speeds in a range based on fint, up to ceiling of max wind speed. If inc is 2, we get frequencies of speeds from 0-2 m/s, 2-4 m/s, and so on...
def freqrange(fint, max):
    lst = []
    inc = 0
    while inc <= max:
        current = 0
        for i in range(len(Tv)):
            speed = Tv[i]
            if inc <= speed <= inc + fint:
                current += 1
        inc += fint
        lst.append(current/n)
    return lst

BinNum = 1000 #How many Bins do you want
incery = max(Tv)/BinNum  #Change increments of frequency range. If 2, we will break speeds into incs of 2 m/s. Frequencies include wind speeds of 0-2 m/s, 2-4 m/s, and on...
Ery = freqrange(incery, max(Tv))
#List of the Ranges for our bins, for table to be printed. NOT BEING USED
def EryBins(inc):
    lst = []
    k = 0
    while k <= max(Tv):
        lst.append("{} - {}".format(k, k+inc))
        k+=inc


yav = sum(Ery) / len(Ery)
#Numerator for R^2 value, also used in RMSE calc, for Maximum Likelihood Method
def numWeib(use):
    num = 0
    w = 0
    if use == 'mlm':
        for j in Ery:
            we = weibmlm(w, w + incery, 'error')
            num += (we - j)**2
            w += incery
        return num
    elif use == 'mm':
        for j in Ery:
            num += (j - weibmm(w, w + incery, 'error'))**2
            w += incery
        return num 

#Denominator for R^2 value, for all methods of weibull distribution, only reliant on frequency of wind speeds which are the same for all methods
def demWeib():
    f = 0
    dem = 0
    for j in Ery:
        dem += (j - yav) ** 2
        f += 1
    return dem 

#Number of "bins" our speeds are broken up into for frequency calculations. It was putinto the list so this is the length of the list
freqbins = len(Ery)

#Error Calculations for maximum likelihood Method
num1 = numWeib('mlm')
demW= demWeib()
R2mlm = 1 -  (num1/demW)  
RMSEmlm = ((1/len(Ery)) * num1) ** 0.5 

#Error Calculations for Moment Method
num2 = numWeib('mm')
R2mm = 1 - (num2/demW)
RMSEmm = ((1/freqbins) * num2) ** 0.5



#Enter Turbine Parameters for Calculations of AEP, power curve, capacity factor, region probability
Ci = float(input("What is the cut-in speed of the turbine? "))
Co = float(input("What is the cut-out speed of the turbine? "))
Pr = float(input("What is the rated power of the turbine? (in kW)"))
Vr = float(input("What is the rated wind speed fo the turbine? "))

#for power curve, normalized power
def normPow(v):
    return ((v**(kmlm)) - (Ci ** (kmlm))) / ((Vr ** (kmlm)) - (Ci ** (kmlm)))

#for capacity factor calculation, normalized power times weibull
def CapCalc(v):
    nP = ((v**(kmlm)) - (Ci ** (kmlm))) / ((Vr ** (kmlm)) - (Ci ** (kmlm)))
    Wei = prob * (kmlm/cmlm) * ((v/cmlm) ** (kmlm - 1)) * np.exp(-(v/cmlm)**kmlm)
    return nP*Wei

#calculate capactiy factor
def capFact():
    Unrated = integrate.quad(CapCalc, Ci, Vr)
    Rated = integrate.quad(weibull, Vr, Co)
    return Unrated[0] + Rated[0]

CF = capFact()
Peavg = Pr * CF
AEP = Peavg * 8760 #Annual energy production, Peavg in kW * hours in year

#Calculate probability density of being in the different regions
def regionWeib(min, max, reg):
    weib = 0
    space = np.linspace(min, max, 1000)  
    if reg == 'unrated':
        #goes to max-0.001 as this region excends from Vci <= V < Vr
        result = integrate.quad(weibull, min, max-0.001)
        return result[0] 
    if reg == 'rated':
        #Goes to max as this region excends from Vr <= V <= Vco
        result = integrate.quad(weibull, min, max)
        return result[0] 
    if reg == 'cutin':
        #goes to max-0.001 as this region excends from 0 <= V < Vci
        result = integrate.quad(weibull, min, max-0.001)
        return result[0]
    if reg == 'cutout':
        #Starts at min + 0.001 as this region excends from Vco < V <= Vmax
        result = integrate.quad(weibull, min+0.001, max)
        return result[0] 
    


#perhaps change this later for general power curve
def genpow():
    vpow = np.linspace(0, 30, 1000)
    lst = []
    for j in vpow:
        speed = j
        if speed < Ci:
            lst.append(0)
        elif Ci <= speed < Vr: #MOdel says the speed should cut off at Vr, which it does but not mathematically
            lst.append(normPow(speed) * Pr)
        elif Vr <= speed <= Co: #Also changed from Vr
            lst.append(Pr)
        elif speed > Co:
            lst.append(0)
    array = np.array(lst)
    return array

freqheader = [EryBins(incery)]


#evenly spaces wind speed over 1000 points based on its minimum and maximum
xv = np.linspace(min(Tv), max(Tv), 1000)

plt.plot(xv, weibmlm(0, 24, 'plot'), label = " Maximum Likelihood Method: k = {}, c = {}, R^2 = {}, RMSE = {}".format(kmlm, cmlm, R2mlm, RMSEmlm))
#plt.plot(xv, weibmm(0, 24, 'plot'), label = " Moment Method: k = {}, c = {}, R^2 = {}, RMSE = {}".format(kmm, cmm, R2mm, RMSEmm)) 
inchist = [i for i in range(math.ceil(max(Tv)))] #x values for histogram, if max speed is 24 m/s, will give frequencies from 1 to 24 in increments of 1, so the first will be speeds from 0-1, then 1-2 and so on.
histweight = freqrange(1, math.ceil(max(Tv) - 1)) #x and weights must be lists of same length, use ceiling of wind speed max to ensure this as we will get frequencies in increments of 1.
plt.hist(inchist, weights = histweight, bins = range(len(inchist))) #plot observed wind speed frequencies
plt.xlabel("Wind Speed at 100m, evenly spaced over 1000 points (m/s)")
plt.ylabel("Weibull PDF")
plt.title("Weibull Distribution on Historical Wind data, 2007-2012, 30 Minute Increments, Parameters Estimated by the Maximum Liklihood Method")
plt.legend()
plt.show()

i = 0
for i in range(math.ceil(max(Tv))):
    print(inchist[i], histweight[i])
    i += 1

print('Unrated region average probability:')
print(regionWeib(Ci, Vr, 'unrated'))

print('Rated region average probability:')
print(regionWeib(Vr, Co, 'rated'))

print('Cut-in average probability:')
print(regionWeib(0, Ci, 'cutin'))

print('Cut-out region average probability:')
print(regionWeib(Co, max(Tv), 'rated'))

print('Capacity Factor')
print(CF)

print('AEP')
print(AEP)

print('Average Expected Output Power')
print(Peavg)

xp = np.linspace(0, 30, 1000)
#Power Curve, Experimental not correct, Not being used
#power = np.piecewise(xp, [xp < Ci, ((Ci <= xp) & (xp < Vr)), ((Vr <= xp) & (xp <= Co)), xp > Co], [0, genpow(), Pr, 0])
pw = genpow()
plt.plot(xp, pw)
#plt.vlines(Ci, 0, Pr, colors = '#008000', linestyles = 'dashed', label = 'cut-in speed')
#plt.vlines(Co, 0, Pr, colors = '#008000', linestyles = 'dashed', label = 'cut-out speed')
plt.xlabel("Wind speed at 100m, evenly spaced over 1000 points (m/s)")
plt.ylabel("Power (kW)")
plt.title("Power Curve")
plt.show()

print("DONE")