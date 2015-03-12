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
    def __init__(self,v=0.6,statespace=("a","t","c","g"),EquilFreq=None,R=None,QMatrix=None,NumSimulations=None):
                     
            self.v=v #branch length
            self.statespace=statespace #my state space of possible states
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
                
            if NumSimulations==None:
                self.NumSimulations=input("Please enter the number of simulations:")
                  
            
    def normalize(self,NormQMatrix=None): 
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
               
            if NormQMatrix==None:
            #now I will use the array method of numpy so I can divide all the values of my matrix by the weighted average
                QArray=np.array(Qmatrix)
                self.NormQMatrix=QArray/waverage 
                NormQMatrix=self.NormQMatrix
        return NormQMatrix
    def discSamp(self,probs):
        import scipy
        """
        Author:Jeremy Brown
        This function samples from a list of discrete events provided in the events argument, using the event
        probabilities provided in the probs argument. These lists must:
        - Be the same length
        - Be in corresponding orders
        Also, the probabilities in probs must sum to 1.
        """
        ranNum = scipy.random.random()
        cumulProbs = []
        cumulProbs.extend([probs[0]])
        for i in range(1,len(probs)):
            cumulProbs.extend([probs[i]+cumulProbs[-1]])
        for i in range(0,len(probs)):
            if ranNum < cumulProbs[i]:
                return self.statespace[i]
        return None
        
    def NumSimulations():
        for i in range(self.NumSimulations): 
            simulate() 
    def simulate(self,MargProbs):
        import scipy
        import random
        #create a list of the marginal probabilties associated with each nucleotide
        MargProbs=[NormQMatrix[0][1]/-NormQMatrix[0][0],NormQMatrix[0][2]/-NormQMatrix[0][0],NormQMatrix[0][3]/-NormQMatrix[0][0],
                   NormQMatrix[1][0]/-NormQMatrix[1][1],NormQMatrix[1][2]/-NormQMatrix[1][1],NormQMatrix[1][3]/-NormQMatrix[1][1],
                   NormQMatrix[2][0]/-NormQMatrix[2][2],NormQMatrix[2][1]/-NormQMatrix[2][2],NormQMatrix[2][3]/-NormQMatrix[2][2],
                   NormQMatrix[3][0]/-NormQMatrix[3][3],NormQMatrix[3][1]/-NormQMatrix[3][3],NormQMatrix[3][2]/-NormQMatrix[3][3]]
        states=[] #this is a list to hold that states of the chain 
        waitTimes=[]
        # Draw a STARTING STATE from the equilibrium frequencies              
      
        randomNum = scipy.random.random()
       
        if randomNum <= self.EquilFreq[0]: 
            states.extend(self.statespace[0]) #should return state a
        elif self.EquilFreq[0]< randomNum <=sum(self.EquilFreq[1]+self.EquilFreq[0]): 
            states.extend(self.statespace[1]) #should return state c
        elif self.EquilFreq[3]> randomNum >self.EquilFreq[1]:
            states.extend( self.statespace[2]) #should return state g
        else:
            states.extend(self.statespace[3]) #should return t
        # Draw a waiting time from the appropriate exponential distribution.
        while sum(waitTimes)<self.v:
            
            for i in states[-1]: #this will pull the last element of the states list
            
              if states[-1] == self.statespace[0]: #we are in statespace a
                  waitTime=random.expovariate(-NormQMatrix[0][0])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[0]) #draw from row 1 of list
                  if newstate==NormQMatrix[0][1]/-NormQMatrix[0][0]:
                      newstate=self.statespace[1]
                      states.extend(newstate)
                  elif newstate==NormQMatrix[0][2]/-NormQMatrix[0][0]:
                      newstate=self.statespace[2]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[3]
                      states.extend(newstate)
                      
              elif states[-1] ==self.statespace[1]: #we are in statespace c
                  waitTime=random.expovariate(-NormQMatrix[1][1])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[1]) #draw from row 2 of list 
                  if newstate==NormQMatrix[1][0]/-NormQMatrix[1][1]:
                      newstate=self.statespace[0]
                      states.extend(newstate)
                  elif newstate==NormQMatrix[1][2]/-NormQMatrix[1][1]:
                      newstate=self.statespace[2]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[3]
                      states.extend(newstates)
                      
              elif states[-1] == self.statespace[2]: #we are in statespace g
                  waitTime=random.expovariate(-NormQMatrix[2][2])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[2])
                  if newstate==NormQMatrix[2][0]/-NormQMatrix[2][2]:
                      newstate=self.statespace[0]
                      states.extend(newstate)
                  elif newstate==NormQMatrix[2][1]/-NormQMatrix[2][2]:
                      newstate=self.statespace[1]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[3]
                      states.extend(newstates)
                      
              else: #this will give us the wait time for "t"
                  waitTime=random.expovariate(-NormQMatrix[3][3])
                  waitTimes.extend(waitTime) 
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[3])
                  if newstate ==NormQMatrix[3][0]/-NormQMatrix[3][3]:
                      newstate=self.statespace[0]
                      states.extend(newstate)
                  elif newstate ==NormQMatrix[3][1]/-NormQMatrix[3][3]:
                      newstate=self.statespace[1]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[2]
                      states.extend(newstate)
              
            # Keep drawing waiting times and new states until the total time exceeds
            # the branch length.
            #Now I need estimate current likelihood of branchlength 
              
        def ProbMatrix():
            import scipy
            ProbM=scipy.linalg.expm(NormQMatrix*self.v)
            
            return ProbM
         def IndivStates(self,Chains):
            """from itertools import *

            for i in izip([1, 2, 3], ['a', 'b', 'c']):
                print i
            
            (1, 'a')
            (2, 'b')
            (3, 'c')
            
            """
         try:
                from itertools import izip as zip
            except ImportError: # will be 3.x series
                pass
                    
            for i in izip(states):
                self.Chains=i
            return self.Chains
            
            
            """Now I need to create something that goes through each row counts the number of values
            (which would depend on the number of simulations) and associates their self.statespace value 
            as it moves to the next index in the row with it's respective self.statespace value"""
