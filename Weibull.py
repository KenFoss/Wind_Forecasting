import types
from numpy.lib import nanfunctions
import scipy.stats as s
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from scipy.special import gamma
from tabulate import tabulate

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


Tv = []

#Values needed for Maximum Likelihood method
vzero = 0
speeds = [0,0,0,0,0,0] #This list gives amount of speeds within increments of 4 m/s, so 0-4m/s, 4-8 m/s, 8-12, 12-16, 16-20, 20-24
prob = 0

#For Error Calculations, generic between methods
yi = [] #frequencies of observed wind speeds at each step


n = 0

#Adding P(v>0) as P, probability wind speed exceeds zero, k is the shape factor and c is scale factor. All calculated on excel.
#potentially will add in incremental wind speed
def weib(p, v, c, k):
    return p * (k/c) * (v/c) ** (k-1) * np.exp(-(v/c)**k)



with open('windspeedtest.csv', 'r') as read_obj:
    reader = csv.reader(read_obj)
    for row in reader:
        n += 1;
        speed = float(row[0]); #We read speeds as dictionary with header, row[0] gets value as string, convert to float 
        Tv.append(speed);
        if speed >= 0:
            vzero +=1
        prob = float(vzero/n) #probability wind speed exceeds zero

        #The following is for the frequencies of observed wind speed data. Will gather frequency of speeds in the 1 m/s ranges for error calc
        if 0 <= speed <= 4:
            speeds[0] += 1
            yi.append(speeds[0])

        elif speed <= 8:
            speeds[1] += 1
            yi.append(speeds[1])
            
        elif speed <= 12:
            speeds[2] += 1
            yi.append(speeds[2])
            
        elif speed <= 16:
            speeds[3] += 1
            yi.append(speeds[3])    
        elif speed <= 20:
            speeds[4] += 1
            yi.append(speeds[4])    
        elif speed <= 24:
            speeds[5] += 1
            yi.append(speeds[5]/n)
        #print("On iteration {}".format(n)) FOR TESTING
            
print("iteration complete for data import")
        


#Weibull distribution over given range for error calculation
kmlm = 0
cmlm = 0

def weibmlm(x, y, use): #Weibull Within Time Range for error calc
    v = []
    lnv = []
    vtok = []
    vlnv = []
    vzero = 0
    krange = 2
    ner = 0 #How many speed points are in the range we are looking for error in
    i=0

    for i in range(len(Tv)): #n value will not be equal to global n if we look in a range
        speed = Tv[i]
        if x <= speed <= y:
            ner += 1

    i = 0

    if ner == 1:
        #With one value in the bin it will scew data, better to return zero
        return 0

    while i < len(Tv):
        speed = Tv[i]
        if x <= speed <= y:  
            if krange < 0: 
                print("i = {}, lnv = {}, vlnv = {}, vtok = {}".format(i, lnv[i-1], vlnv[i-1], vtok[i-1]))
            lnv.append(np.log(speed))
            vtok.append(speed ** krange)  #current wind speed to power of k, for scale factor calc
            vlnv.append((speed ** krange) * (np.log(speed)))      
            v.append(speed)
            krange = (((sum(vlnv) / sum(vtok)) - (sum(lnv) / ner)) ** (-1)) #This code will find the k factor, must be solved iteratively with an original k value guess
            #print("speed = {}, krange = {}. sumvlv = {}, sumvtok = {}, sumlnv = {}".format(speed, krange, sum(vlnv), sum(vtok), sum(lnv)))  #TESTING FUNCTION
            if math.isnan(krange): #why is krange coming out as nan, TESTING FUNCTION, problem solved n was not implemented
                print("krange is nan at iteration {}. sumvlv = {}, sumvtok = {}, sumlnv = {}".format(krange, sum(vlnv), sum(vtok), sum(lnv)))
        i += 1
        #print("Iteration in Weibmlm, on iteration {}, k value is {}".format(i, krange))
    if use == 'error':
        print("Iteration complete maximum likelihood method, error calculation for range {} - {}".format(x, y))
        if len(v) != 0 or len(v):
            crange = ((1/ner) * sum(vtok)) ** (1/krange)  #This code will find c factor, can be solved explicitly
            weiberror = 0
            for i in v:
                weiberror +=  prob * (krange/crange) * ((i/crange) ** (krange-1)) * np.exp(-(i/crange)**krange) #weibull distribution, prob is the probability the wind speed is not zero, then we have c and k factors
            return weiberror/ner
        else:
            return 0
    elif use == 'plot':
        print("iteration complete for plotting weibull, maximum likelihood method")
        crange = ((1/ner) * sum(vtok)) ** (1/krange)  #This code will find c factor, can be solved explicitly
        global kmlm #Set global kmlm, to change the value of kmlm outside of this function...same for cmlm
        kmlm = krange #Set global variables for later display
        global cmlm 
        cmlm = crange #''
        xv = np.linspace(min(Tv), max(Tv), 1000)
        return prob * (krange/crange) * ((xv/crange) ** (krange-1)) * np.exp(-(xv/crange)**krange)

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
        print("Iteration complete weibull moment method, error calculation for range {} - {}".format(x, y))
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
    #NOT TESTED

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


incery = 4  #Change increments of frequency range. If 2, we will break speeds into incs of 2 m/s. Frequencies include wind speeds of 0-2 m/s, 2-4 m/s, and on...
Ery = freqrange(incery, max(Tv))
#List of the Ranges for our bins, for table to be printed.
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
#RMSEmlm = ((1/6) * num1) ** 0.5 

#Error Calculations for Moment Method
num2 = numWeib('mm')
R2mm = 1 - (num2/demW)
RMSEmm = ((1/freqbins) * num2) ** 0.5

Cap = input("What is the turbines power coefficient? Entery nothing if this data is not avalible. ")
Cp = 0
if Cap == '':
    print("Since you do not have a power coefficient, the approximate cubic method will be used. The power coefficient has been set to betz limit of 0.593.")
    Cp = 0.593 #Betz limit, physacist betz calculated no turbine could convert more then 59.3% of kinetic energy from wind into mechanical energy from turbine
else:
    Cp = float(Cap)

A = float(input("What is the rotor area of the turbine (in m^2)? "))
Ci = float(input("What is the cut-in speed of the turbine? "))
Co = float(input("What is the cut-out speed of the turbine? "))
Pr = float(input("What is the rated power of the turbine? "))
Vr = float(input("What is the rated wind speed fo the turbine? "))

#This may be deleted, ideally I would like to have multiple points of air density data to average
p = float(input("What is the air density? "))

def cubpow(p, A, Cp):
    vpow = np.linspace(0, 23, 1000)
    lst = []
    for j in vpow:
        speed = j
        if speed < Ci:
            lst.append(0)
        elif Ci <= speed < Vr: #MOdel says the speed should cut off at Vr, which it does but not mathematically
            lst.append(((0.5) * p * A * Cp * (speed ** 3)) * (0.000001))
        elif Vr <= speed <= Co: #Also changed from Vr
            lst.append(Pr)
        elif speed > Co:
            lst.append(0)
    array = np.array(lst)
    return array

freqheader = [EryBins(incery)]


#evenly spaces wind speed over 1000 points based on its minimum and maximum
xv = np.linspace(min(Tv), max(Tv), 1000)

#plt.plot(xv, weibmlm(0, 24, 'plot'), label = " Maximum Likelihood Method: k = {}, c = {}, R^2 = {}, RMSE = {}".format(kmlm, cmlm, R2mlm, RMSEmlm))
plt.plot(xv, weibmm(0, 24, 'plot'), label = " Moment Method: k = {}, c = {}, R^2 = {}, RMSE = {}".format(kmm, cmm, R2mm, RMSEmm)) 
inchist = [i for i in range(math.ceil(max(Tv)))] #x values for histogram, if max speed is 24 m/s, will give frequencies from 1 to 24 in increments of 1, so the first will be speeds from 0-1, then 1-2 and so on.
histweight = freqrange(1, math.ceil(max(Tv) - 1)) #x and weights must be lists of same length, use ceiling of wind speed max to ensure this as we will get frequencies in increments of 1.
plt.hist(inchist, weights = histweight, bins = range(len(inchist))) #plot observed wind speed frequencies
plt.xlabel("Wind Speed at 100m, evenly spaced over 1000 points (m/s)")
plt.ylabel("Weibull PDF")
plt.title("Weibull Distribution on Historical Wind data, 2007-2012, 30 Minute Increments, Parameters Estimated by the Maximum Liklihood Method")
plt.legend()
plt.show()

xp = np.linspace(0, 30, 1000)
#power = np.piecewise(xp, [xp < Ci, ((Ci <= xp) & (xp < Vr)), ((Vr <= xp) & (xp <= Co)), xp > Co], [0, cubpow(p,A,Cp), Pr, 0])
pw = cubpow(p,A,Cp)
plt.plot(xp, pw)
plt.vlines(Ci, 0, Pr, colors = '#008000', linestyles = 'dashed', label = 'cut-in speed')
plt.vlines(Co, 0, Pr, colors = '#ff0000', linestyles = 'dashed', label = 'cut-in speed')
plt.xlabel("Wind speed at 100m, evenly spaced over 1000 points (m/s)")
plt.ylabel("Power (MW)")
plt.show()

print("DONE")