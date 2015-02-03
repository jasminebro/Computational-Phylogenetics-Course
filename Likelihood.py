from __future__ import division
import matplotlib.pyplot as plt
from scipy.stats import binom

n = 5
p = 0.5 # Change this and repeat
data = binom.rvs(n,p)

data =4 # Supply observed number of successes here.
numTrials = 5

def multiple (max,min):
    total = 1 #set total to 1, not zero
    if min <= 0 or max <=0: #make sure your range doesnt equal zero
        print "your range contains the number 0 and will therefor return zero"
    else:
        for num in range(max,min-1,-1): #this works because the second argument is non-inclusive, so even if it's 0 then it will work
            total *= num #multiply all numbers together
        return total 
 

def BiCoefficient(n,k): 
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
    p=PMF(4,5,x)
    likelihood.append(p)
print (likelihood)


# Find the maximum likelihood value of p (at least, the max in this set)
likelimax=likelihood.index(max(likelihood))

print("the value of p with the max likelihood is", Pvalues[likelimax],"its likelihood being", max(likelihood))
plt.scatter(Pvalues,likelihood)
plt.xlabel('p values')
plt.ylabel('Likelihood')
