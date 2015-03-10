"""
Created on Thu Mar 05 13:52:33 2015
@author: jasminebrown

I am going to re-do my object oriented assignment following the directions posted on 3/5 
and by adding the normalization step at the beginning 

My matrix should look like this 
   A       C        G      T
A

C   Place your values here! 

G

T
"""
class ContinMarkov(object):
    def __init__(self,v=12.0,StateSpace=("a","t","c","g"),NumSimulations=1,EquilFreq=None,R=None,QMatrix=None):
                     
            self.v=v #branch length
            self.StateSpace=StateSpace #my state space of possible states
            self.NumSimulations=NumSimulations #number of simulations I will do
            self.ChainStates=[] #empty list to hold states as they change
            self.WaitingTimes=[]#empty list to hold waiting times between state changes
            #note:make sure the length of the total chain is the length of the # of wait times
            
            if EquilFreq==None:    
                FreqA=float(input("Please enter the frequency for A: "))#0.27 
                FreqC=float(input("Please enter the frequency for C: "))#0.34
                FreqG=float(input("Please enter the frequency for G: "))#0.20
                FreqT=round(1-(FreqA+FreqC+FreqG),2)
                self.EquilFreq=[FreqA,FreqC,FreqG,FreqT]
            if R==None: #now I am going to define individual r values for transition
                #user can set these values 
                rAC=0.45
                rAG=0.67
                rAT=0.76
                rCG=0.39
                rCT=0.59
                rGT=0.43
                self.R=[rAC,rAG,rAT,rCG,rCT,rGT] #making the values a list so I can index them
            if QMatrix==None:
                self.QMatrix=normalize()
            
    def normalize(self,NormQMatrix): 
        """I am creating a function to normalize my matrix before simulating transtions""" 
        if sum(self.EquilFreq) ==1:
            print "Your equilibrium frequency rates equal 1"
            print "this is your frquecy list for A,C,G and T", EquilFreq
        else:
            print "Your frequencies DO NOT EQUAL 1...better start over!"
        #this is my lists of lists that contains my q matrix values
        Qmatrix=[[-1*(R[0]*freqlist[1]+R[1]*freqlist[2]+R[2]*freqlist[3]),R[0]*freqlist[1],R[1]*freqlist[2],R[2]*freqlist[3]],
                 [R[0]*freqlist[0],-1*(R[0]*freqlist[0]+R[3]*freqlist[2]+R[4]*freqlist[3]),R[3]*freqlist[2],R[4]*freqlist[3]],
                 [R[1]*freqlist[0],R[3]*freqlist[1],-1*(R[1]*freqlist[0]+R[3]*freqlist[1]+R[5]*freqlist[3]),R[5]*freqlist[3]],
                 [R[2]*freqlist[0],R[4]*freqlist[1],R[5]*freqlist[2],-1*(R[2]*freqlist[0]+R[4]*freqlist[1]+R[5]*freqlist[2])]] 
                 
        for l in Qmatrix:
            #map(a function,a sequence)---applies your function over the entire sequence 
            print ', '.join(map(str, l))
            #my Qmatrix will be displayed in a "matrix" like format(each row on it's own line)
            #now I need to calculate the weighted average of rates equation
            waverage=(freqlist[0]*-Qmatrix[0][0])+(freqlist[1]*-Qmatrix[1][1])+(freqlist[2]*-Qmatrix[2][2])+(freqlist[3]*-Qmatrix[3][3])
            print "My weighted average of rates is:" ,waverage 
            import numpy as np
            #now I will use the array method of numpy so I can divide all the values of my matrix by the weighted average
            QArray=np.array(Qmatrix)
            NormQMatrix=QArray/waverage    
        return NormQMatrix
        
        def simulate():
            import scipy
            import random
            #create a list of the marginal probabilties associated with the possible transitions(6)
            MargProbs=[(self.R[0]*self.EquilFreq[1]/-NormQMatrix[1][1]),(self.R[1]*self.EquilFreq[2]/-NormQMatrix[2][2]),
                       (self.R[2]*self.EquilFreq/[-NormQMatrix[3]][3]),(self.R[3]*self.EquilFreq[2]/-NormQMatrix[2[2]]),
                       (self.R[4]*self.EquilFreq[3]/-NormQMatrix[3][3]),(self.R[5]*self.EquilFreq[3]/-NormQMatrix[3][3])]
            states=[] #this is a list to hold that states of the chain 
            waitTimes=[]
            # Draw a STARTING STATE from the equilibrium frequencies
            randomNum = scipy.random.random()
           
            if randomNum <= self.EquilFreq[0]: 
                states.extend(self.statespace[0]) #should return state a
            elif self.EquilFreq[0]< randomNum <=self.EquilFreq[1]: 
                states.extend(self.statespace[1]) #should return state c
            elif self.EquilFreq[3]> randomNum >self.EquilFreq[1]:
                states.extend( self.statespace[2]) #should return state g
            else:
                states.extend(self.statespace[3]) #should return t
            # Draw a waiting time from the appropriate exponential distribution.
            while sum(waitTimes)<self.v:
                for i in states[-1]: #this will pull the last element of the states list
                  if states[-1] == "a": 
                      waitTime=random.expovariate(-NormQMatrix[0][0])
                      waitTimes.extend(waitTime)
                      # Draw a new state from the marginal probabilities associated with the current state.
                      newstate=max(MargProbs[0],MargProbs,[1] and MargProbs[2])
                      states.extend(newstate)
                  elif states[-1] =="c":
                      waitTime=random.expovariate(-NormQMatrix[1][1])
                      waitTimes.extend(waitTime)
                      # Draw a new state from the marginal probabilities associated with the current state.
                      newstate=max(MargProbs[0],MargProbs[3] and MargProbs[4])
                      states.extend(newstate)
                  elif states[-1] == "g":
                      waitTime=random.expovariate(-NormQMatrix[2][2])
                      waitTimes.extend(waitTime)
                      # Draw a new state from the marginal probabilities associated with the current state.
                      newstate=max(MargProbs[1], MargProbs[3] and MargProbs[5])
                      states.extend(newstate)
                  else: #this will give us the wait time for "t"
                      waitTime=random.expovariate(-NormQMatrix[3][3])
                      waitTimes.extend(waitTime) 
                      # Draw a new state from the marginal probabilities associated with the current state.
                      newstate=max(MargProbs[2], MargProbs[4] and MargProbs[5])
                      states.extend(newstate)
                  
            # Keep drawing waiting times and new states until the total time exceeds
            # the branch length.
        
