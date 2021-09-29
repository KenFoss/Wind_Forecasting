import types
import scipy.stats as s
import numpy as np
import matplotlib.pyplot as plt
import math
import csv

#This Model takes a csv files of wind speeds in the first column, in m/s, and calculates a Weibull distribution with R^2 and RMSE Error. 
#Future additions include other methods of calculations for the k and c parameters, and other models such as the Rayleigh Distribution


#change wroking directory to find csv
#When debugging eliminate path and os.chdir(path), comment them out, it will find the file. It does not when run normally.
#The above applies to running this code in the command prompt as well
import os
path = ('Desktop\Python\WeibullOffline')  #set path
os.chdir(path)
print("Current working directory: {0}".format(os.getcwd()))


kmlm = 2  #This guess is for the maximum likelihood Method, from papers k is a suitable initial guess.
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
        print("On iteration {}".format(n))
            
print("iteration complete for data import")

#ERROR Calculations for mlmv (R^2) (RMSE)
yav = sum(yi)/n #Average frequency of wind speed 

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
    n = 0
    i=0
    while i < len(Tv):
        speed = Tv[i]
        n += 1 #Must increment n again, as we go through the list generated from our csv wind speed data
        if x <= speed <= y:  
            lnv.append(np.log(speed))
            vtok.append(speed ** krange)  #current wind speed to power of k, for scale factor calc
            vlnv.append((speed ** krange) * (np.log(speed)))      
            v.append(speed)
            krange = ((sum(vlnv) / sum(vtok)) - (sum(lnv) / n) ** (-1)) #This code will find the k factor, must be solved iteratively with an original k value guess
            #print("speed = {}, krange = {}. sumvlv = {}, sumvtok = {}, sumlnv = {}".format(speed, krange, sum(vlnv), sum(vtok), sum(lnv)))  TESTING FUNCTION
            #if math.isnan(krange): #why is krange coming out as nan, TESTING FUNCTION, problem solved n was not implemented
                #print("krange is nan at iteration {}. sumvlv = {}, sumvtok = {}, sumlnv = {}".format(krange, sum(vlnv), sum(vtok), sum(lnv)))
        i += 1
        #print("Iteration in Weibmlm, on iteration {}, k value is {}".format(i, krange))
    if use == 'error':
        print("Iteration complete, error calculation for range {} - {}".format(x, y))
        if len(v) != 0:
            crange = ((1/n) * sum(vtok)) ** (1/krange)  #This code will find c factor, can be solved explicitly
            return prob * (krange/crange) * (y/crange) ** (krange-1) * np.exp(-(y/crange)**krange) #weibull distribution, prob is the probability the wind speed is not zero, then we have c and k factors
        else:
            return 0
    elif use == 'plot':
        print("iteration complete for plotting")
        crange = ((1/n) * sum(vtok)) ** (1/krange)  #This code will find c factor, can be solved explicitly
        global kmlm #Set global kmlm, to change the value of kmlm outside of this function...same for cmlm
        kmlm = krange #Set global variables for later display
        global cmlm 
        cmlm = crange #''
        xv = np.linspace(min(Tv), max(Tv), 1000)
        return prob * (krange/crange) * (xv/crange) ** (krange-1) * np.exp(-(xv/crange)**krange)


#Error Calculations for Maximum Likelihood Method
def freq(speed):
    return speeds[speed]/n

#Numerator for R^2 value, also used in RMSE calc, for Maximum Likelihood Method
def nummlm():
    f = 0
    wx = 0
    wy = 4
    num = 0
    while f <= 5:
        num += (freq(f) - weibmlm(wx, wy, 'error'))
        f += 1
        wx += 4
        wy += 4
    return num ** 2

#Denominator for R^2 value, for Maximum Likelihood Method
def demmlm():
    f = 0
    dem = 0
    while f <= 5:
        dem += (freq(f) - yav)
        f += 1
    return dem ** 2

#Error Calculations for maximum likelihood Method
num1 = nummlm()
R2mlm = 1 -  (num1/demmlm())  
RMSEmlm = ((1/6) * num1) ** 0.5 






#evenly spaces wind speed over 1000 points based on its minimum and maximum
xv = np.linspace(min(Tv), max(Tv), 1000)

plt.plot(xv, weibmlm(0, 24, 'plot'), label = "k = {}, c = {}, R^2 = {}, RMSE = {}".format(kmlm, cmlm, R2mlm, RMSEmlm))
plt.xlabel("Wind Speed at 100m, evenly spaced over 1000 points (m/s)")
plt.ylabel("Weibull PDF")
plt.title("Weibull Distribution on Historical Wind data, 2007-2012, 30 Minute Increments, Parameters Estimated by the Maximum Liklihood Method")
plt.legend()
plt.show()
print("DONE")