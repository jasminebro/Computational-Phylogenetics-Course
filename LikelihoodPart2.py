"""
Created on Mon Feb 02 16:11:07 2015

@author: Jasmine Brown
"""
from __future__ import division
import matplotlib.pyplot as plt
from scipy.stats import binom

#n = 5
#p = 0.5 # Change this and repeat
#data = binom.rvs(n,p)

data =4 # number of success
numTrials = 5 #number of trials

def multiple (max,min):
    total = 1 #set total to 1, not zero
    if min <= 0 or max <=0: #make sure your range doesnt equal zero
        print "your range contains the number 0 and will therefore return zero"
    else:
        for num in range(max,min-1,-1): #this works because the second argument is non-inclusive, so even if it's 0 then it will work
            total *= num #multiply all numbers together
        return total 
 

def BiCoefficient(n,k): #def for my bionomial coefficient equation
    y=multiple(k,1)
    max=n
    min=(n-k+1)
    j=multiple(max,min)
    binomial=j/y
    #note that n should be your max value and k the min. you want to see how many sets are in the larger
    return binomial


def PMF(k,n,p): 
    ProbMassFunction= BiCoefficient(n, k)*pow(p,k)*pow((1-p),(n-k))
    
    return ProbMassFunction
    
#Now we need to calculate likelihoods for a series of different values for p to compare likelihoods. There
#are an infinite number of possible values for p, so let's confine ourselves to steps of 0.05 between 0 and 1.

# Set up a list with all relevant values of p
Pvalues=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]
# Calculate the likelihood scores for these values of p
likelihood=[]

for x in Pvalues:
    p=PMF(12,20,x)
    likelihood.append(p)
print (likelihood)
#I used this equation for both parts the 20 and the 5 trials and observed the differences in having a larger set size

# Find the maximum likelihood value of p (at least, the max in this set)
likelimax=likelihood.index(max(likelihood))

print("the value of p with the max likelihood is", Pvalues[likelimax],"its likelihood being", max(likelihood))
plt.scatter(Pvalues,likelihood)
plt.xlabel('p values')
plt.ylabel('Likelihood')
"""
Sometimes it will not be feasible or efficient to calculate the likelihoods for every
value of a parameter in which we're interested. Also, that approach can lead to large
gaps between relevant values of the parameter. Instead, we'd like to have a 'hill
climbing' function that starts with some arbitrary value of the parameter and finds
values with progressively better likelihood scores. This is an ML optimization
function. There has been a lot of work on the best way to do this. We're going to try
a fairly simple approach that should still work pretty well, as long as our likelihood
surface is unimodal (has just one peak). Our algorithm will be:
(1) Calculate the likelihood for our starting parameter value (we'll call this pCurr)
(2) Calculate likelihoods for the two parameter values above (pUp) and below (pDown)
our current value by some amount (diff). So, pUp=pCurr+diff and pDown=pCurr-diff. To
start, set diff=0.1, although it would be nice to allow this initial value to be set
as an argument of our optimization function.
(3) If either pUp or pDown has a better likelihood than pCurr, change pCurr to this
value. Then repeat (1)-(3) until pCurr has a higher likelihood than both pUp and
pDown.
(4) Once L(pCurr) > L(pUp) and L(pCurr) > L(pDown), reduce diff by 1/2. Then repeat
(1)-(3).
(5) Repeat (1)-(4) until diff is less than some threshold (say, 0.001).
(6) Return the final optimized parameter value.
Write a function that takes some starting p value and observed data (k,n) for a
binomial as its arguments and returns the ML value for p.
To write this function, you will probably want to use while loops. The structure of
these loops is
while (someCondition):
code line 1 inside loop
code line 2 inside loop
As long as the condition remains True, the loop will continue executing. If the
condition isn't met (someCondition=False) when the loop is first encountered, the
code inside will never execute.
If you understand recursion, you can use it to save some lines in this code, but it's
not necessary to create a working function.
"""
# Write a function that finds the ML value of p for a binomial, given k ann. 
#Something is wrong with my equations below. It runs now but nothing is being returned
def HillClimber (n,k,pCurr, difference): 
    pUp= pCurr + difference
    pDown= pCurr - difference 
    PMFCurr=PMF(n,k,pCurr)
    PMFUp=PMF(n,k,pUp)
    PMFDown=PMF(n,k,pDown)
    
    while difference > 0.001:  #Try to make while difference > ).001  
       
        if  PMFUp> PMFCurr: #checking is a value above or below is better 
        #make it so you can update the p value
            pCurr=pUp #changed so it calculates the up value I had a random PMF function here
            PMFCurr=PMF(n,k,pCurr)
            pUp=pCurr+difference 
            PMFUp=PMF(n,k,pUp)
            
        elif PMFDown > PMFCurr:
            pCurr=pDown
            PMFCurr= PMF(n,k,pCurr)
            pDown=pCurr-difference
            PMFDown=PMF(n,k,pDown)
        else:
            difference *=0.5 #this will split the pCurr value in half which reduces the values and lets me look at small values around p
    return pCurr
        
print HillClimber(5,4,.74,.1)      
        
"""
In the exercise above, you tried to find an intuitive cutoff for likelihood ratio
scores that would give you a reasonable interval in which to find the true value of
p. Now, we will empirically determine one way to construct such an interval. To do
so, we will ask how far away from the true value of a parameter the ML estimate
might stray. Use this procedure: (1) start with a known value for p, (2) simulate
a bunch of datasets, (3) find ML parameter estimates for each simulation, and then
(4) calculate the likelihood ratios comparing the true parameter values and the ML
estimates. When you do this, you will be constructing a null distribution of
likelihood ratios that might be expected if the value of p you picked in (1)
was true. Note that the ML values for these replicates are very often greater than
L(true value of P), because the ML value can only ever be >= L(true value). Once
you have this distribution, find the likelihood ratio cutoff you need to ensure
that the probability of seeing an LR score that big or greater is <= 5%.
"""
# Set a starting, true value for p
trueP = .74

# Simulate 1,000 datasets of 200 trials from a binomial with this p
# If you haven't already done so, you'll want to import the binom class from scipy:
#from scipy.stats import binom #Importing this so I can use the randon variable function
data=[] #creating a list to store my randome variables generated

for BinomRVS in range (0,1000): #ranging from 0 to 1000 to get 1000 data sets 
    BinomRVS=binom.rvs(200,trueP) # I will do 200 trials for the true p value I assigned 
    data.append(BinomRVS) 
print data

     
# Now find ML parameter estimates for each of these trials
MLdata=[]  #list created to hold the Maximum likelihood values 
for x in data: 
    ML=PMF(x,200,trueP) #for every value in the random variable data set will be the success k values 200 is n my trials and true p is my p value using my PMF def from above
    MLdata.append(ML)
    print MLdata #append the new list of ML data and print it 
# Calculate likelihood ratios comparing L(trueP) in the numerator to the maximum
# likelihood (ML) in the denominator. Sort the results and find the value
# corresponding to the 95th percentile.
 
Ratios=[] # a list to hold te ratios of the values I got above 

for r in range(0,200): #going from a range of 0 to number of trials (my n value)
    LikeTrueP=PMF(data[r],200,trueP) # grabbing the likelihood of the true P value
    Ratio=LikeTrueP/MLdata[r] #the ratio will be the likelihood of the true p divided by the MLdata from the table(each value will be r for the loop)
    Ratios.append(Ratio) #appending the ratios to the list 
    
Ratios.sort() #this sorts my list in numerical order
        
     
# Now, convert the likelihood ratios (LRs) to -2ln(LRs) values.
#OMG Now I'm getting all negative zeros here....something must be wrong with a preceeding equation 
import math
#this should be equivalent to the chi squared data 
LogRatios=[]
for j in range(0,200) :
    LikeliLog=-2*(math.log(Ratios[j]))
    LogRatios.append(LikeliLog)
print LogRatios
# Find the 95th percentile of these values. Compare these values to this table:
# https://people.richland.edu/james/lecture/m170/tbl-chi.html. In particular, look
# at the 0.05 column. Do any of these values seem similar to the one you calculated?
# Any idea why that particular cell would be meaningful?
     
 #now I'm going to do the percentile using numpy since that's the easiest method I found
          
import numpy 
#this part works fine...used numpy to take the 95 percentile of the ratios 
percentage=numpy.percentile(Ratios,95)

print "this is my percentile:",percentage 
     
     
     
     
# Based on your results (and the values in the table), what LR statistic value
# [-2ln(LR)] indicates that a null value of p is far enough away from the ML value
# that an LR of that size is <=5% probable if that value of p was true?
# Using this cutoff, what interval might you report for the 5- and 20-trial data
# sets above?
# We've talked in previous classes about two ways to interpret probabilities. Which
# interpretation are we using here to define these intervals?
