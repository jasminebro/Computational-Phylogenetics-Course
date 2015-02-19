"""
Created on Tue Feb 10 11:22:01 2015

@author: jasminebrown
"""
"""
Recall from your reading that any irreducible and aperiodic Markov chain has a
stationary distribution. To convince ourselves that things will converge for
such a chain with arbitrary transition probabilities, let's give it a try.
Work in pairs for this. It's more fun to be social.
"""
# Paste your Markov chain simulation function below, where the starting state
# is drawn with uniform probability from all possible states. Remember to also
# copy any import statements or other functions on which your simulator is
# dependent.
def DiscreteDistr(events,probabilities):
    """This is my discreet probability function""" 
    #you need to import random to get arbritrary random numbers generated
    import random
    x= random.choice(events)
    #this pulls a random number from a list 
    #the lists must be the same length and I should use the lists index to match
    #events with its proper probability 
    index=events.index(x)
    #this indexes the random events 
    probability=probabilities[index]
    #this is the probabilities of the events given 
    
    return "this is the event", x, "and this is the probability", probability

def Markov(states,MatrixData,steps):
    """Now I am going to create a function for a Markov Chain. states can use a list for the space. MatrixData
    can be used to represent a list of tuples that represent the transitions, steps is how many times this will run""" 
    import random
    list=[]
    #I am making a list that will hold my current states as they transition     
    Currently=random.choice(states)
    #choice chooses between whatever is in that space, my choices are A and B from the Chain list above
  
    for i in range (steps):
        #for some variable in a range of steps that will be given      
        if Currently ==states[0]: #if it equals A
            Currently= DiscreteDistr(states,MatrixData[0])[1] #I need the 1 index to get the x value from my discrete function above
        elif Currently ==states[1]: #if it equals B
            Currently =DiscreteDistr(states,MatrixData[1])[1]
        list.append(Currently)
    return list 


# Define a 2x2 transition matrix. For fun, don't make all the probabilities
# equal. Also, don't use any 0s or 1s (to make sure the chain is irreducible
# and aperiodic).
Chain=("A","B")
Matrix=[[0.6,0.4], [0.3,0.7]]
# Simulate a single chain for three time steps and print the states

print Markov(Chain,Matrix,3)
#[A B B]
# Analytically calculate the progression of states for this chain.
# Calculate the probability of observing the state in step 3, given the initial
# state in step 1 (i.e., as if you didn't know the state in step 2).
"""

P(x1=state1,[x2=state2,x3=state3])
P(A)*P(B|A)
P(state1)*P(state2,state3|state1)
P(state1)* P(state3|state1, state2)* P(state2|state1) ---state 1 disappears because state3 depends on 2 not 1 according to the memeoryless fucntiong of markov
P(state1)* P(state2|state1) *P(state3|state2)
#WHAT IF WE DIDN'T KNOW THE MIDDLE B existed??? 
#can only be AAB or ABB 
"""
"""
#find the probabilities of both of these combinations #n step transition probability in the reading 
P(A,B,B)=0.5*0.7*0.7 
#the first probability if 0.5 because of having a uniform distribution

P(A,?,B)=(0.5*0.3*0.7)+(0.5*0.7*0.7)=0.35
P(A,A,B)+P(B,B,B)
#this is the probability of both sets

# Now think of the chain progressing in the opposite direction. What is the
# probability of the progression through all 3 states in this direction? How
# does this compare to the original direction?
#This means the probability of
P(B,A,A)+P(B,B,B)=(0.5*0.3*0.3)+(0.5*0.7*0.7)=0.245

The probabilities are different 
"""

# Try the same "forward" and "reverse" calculations as above, but with this
# transition matrix:
revMat = [[0.77,0.23], [0.39,0.61]]
# and these starting frequencies for "a" and "b"
freqa = 0.63 
freqb = 0.37
import numpy as np
import random
mt=np.matrix(revMat)

ran=random.random()

if ran <= freqa:
    i=0
elif ran >= freqb:
    i=1
print i

def transitionmatrix(matrix,i,k,n):
    mat=np.matrix(mt)
    transition=mat**n
    probability=(transition[i,k])
    print "my matrix is" +str(transition)
    return "The probability is"+ str(probability)
    
transitionmatrix(revMat,i,1,1000)
# What is (roughly) true about these probabilities?
# Simulate 1,000 replicates (or 10K if your computer is fast enough) of 25
# steps. What are the frequencies of the 2 states across replicates through time?
# NOTE: Here is a function that reports the frequencies of a state through time
# for replicate simulations. You'll need to do this several times during this exercise.
def mcStateFreqSum(sims,state=revMat):
    """
    Pass this function a list of lists. Each individual list should be the
    states of a discrete-state Markov chain through time (and all the same
    length). It will return a list containing the frequency of one state
    ("a" by default) across all simulations through time.
    """
    freqs = []
    for i in range(len(sims[0])): # Iterate across time steps
        stateCount = 0
    for j in range(len(sims)): # Iterate across simulations
        if sims[j][i] == state:
            stateCount += 1
            freqs.extend([float(stateCount)/float(len(sims))])
            return freqs
# Run replicate simulations
# Summarize the frequency of one state through time
# What do you notice about the state frequencies through time? Try another round
# of simulations with a different transition matrix. How do the state freq.
# values change?
# Now, calculate a vector of probabilities for the focal state (e.g., 'a')
# based on the transition matrix directly (not by simulation). How do these
# values compare to the simulated frequencies?
