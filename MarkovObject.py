# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 09:40:36 2015

@author: Jasminebrown
"""
"""
try creating a Markov chain object based on the code you wrote in Exercise 4.
 What types of information might you want to associate with a Markov chain 
 (e.g., transition matrix, state space, etc.)? What types of methods might you want a Markov chain object 
 to have (e.g., run a simulation, give me the state at some particular time step)?
 We’ll do some brainstorming in class. Creating this object shouldn’t actually involve writing much new code. 
 It’s just a matter of reorganizing code into a cohesive, logical unit.
 """
 
class ContinMarkovChain(object):
    """this is me defining a new class"""
    def __init__(self,QMatrix=[[-1.59,0.54,0.27,0.78],[0.67,-1.87,0.45,0.75],[0.78,0.51,-1.61,0.32],[0.32,0.75,0.45,-1.52]],v=12.0,StateSpace=("a","t","c","g"),NumSimulations=1):
        self.QMatrix=QMatrix #the transition matrix
        self.v=v #branch length
        self.StateSpace=StateSpace #my state space of possible states
        self.NumSimulations=NumSimulations #number of simulations I will do
        self.ChainStates=[] #empty list to hold states as they change
        self.WaitingTimes=[] #empty list to hold waiting times between state changes 
        #note:make sure the length of the total chain is the length of the # of wait times 
                      
    def Simulate(self):
        """This is my Markov simulation method"""
        import numpy as np
        import scipy
        
        Q= np.array(self.QMatrix)
        if np.array(self.QMatrix) is Q:
            print "You have turned your QMatrix into an Array"
        else:
            print "You DO NOT have an Array! Go fix it" 
      
        StationaryProbs=scipy.linalg.expm(self.QMatrix*self.v)
        #you need stationary probabilities for continuous markov chains.
        
        def ProbFunction(self,events=StateSpace,StationaryProbs):
            """This is my probability function. This should pull random values for my Q matrix.""" 
            #you need to import random to get arbritrary random numbers generated
            import numpy
            x= numpy.random.choice(self.events)
            #this pulls a random number from a list
            #the lists must be the same length and I should use the lists index to match
            #events with its proper probability 
            index=self.events.index(x)
            #this indexes the random events 
            probability=self.StationaryProbs[index]
            #this is the probabilities of the events given 
    
            return "this is the event", x, "and this is the probability", probability
        
        
def Markov(states,MatrixData,steps):
    import random
    """Now I am going to create a function for a Markov Chain. states can use a list for the space. MatrixData
    can be used to represent a list of tuples that represent the transitions, steps is how many times this will run""" 
    if steps < 0: #dummy check for invalid value for steps
        return 'Invalid number of steps. Must be positive integer.'

    initialState = random.choice(states) #draws the first state with equal probability
    markov_chain = [] #creates the list to store the Markov chain as it transitions
    markov_chain.append(initialState) #appends the first state to the Markov chain
    n = 0 #sets the starting index of the Markov list
    while n <= steps-1: #sets the number of "steps" to run the chain
        if markov_chain[n] == states[0]: #checks if the current element is state[0]
            nextState = DiscreteDistr(states, MatrixData[0])[1] #if it's true, draws the next state based on the transition probabilities of state[0]
            markov_chain.append(nextState) #appends it to the Markov chain
    else:
        nextState = DiscreteDistr(states, MatrixData[1])[1]
        markov_chain.append(nextState)
    n += 1 #increases the value of n to look at the next index
    return markov_chain
                
            
    def StateTime(self):
        
        
    def EndFreq(self): # this is the method for calculating the ending frequencies
        
    def ProbHistory(self): # this is the method for probability histories
        
        
    def MarginalProb(self): # this is the method for the marginal probabilities
