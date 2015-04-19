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
import scipy
import random
import numpy as np
class ContinMarkov(object):
    def __init__(self,v=0.6,MargProbs=None,statespace=["a","t","c","g"],freqlist=[0.27,0.34,0.20,0.19],R=[0.40,0.54,0.43,0.39,0.30,0.43],QMatrix=None,NormQMatrix=None,NumSimulations=2):
        self.v=v #branch length
        self.statespace=statespace #my state space of possible states
        self.NumSimulations=NumSimulations #number of simulations I will do
        self.states=[] #empty list to hold states as they change
        self.WaitTimes=[]#empty list to hold waiting times between state changes
        #note:make sure the length of the total chain is the length of the # of wait times
        self.freqlist=freqlist
        self.R=R
        if QMatrix==None:
            self.QMatrix=[]
        if NormQMatrix==None:
            self.NormQMatrix=self.normalize()
        if MargProbs==None:
            self.MargProbs=self.MarginalProbs()
        
    def discSamp(self,freqlist,statespace):
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
        cumulProbs.extend([self.freqlist[0]])
        for i in range(1,len(self.freqlist)):
            cumulProbs.extend([self.freqlist[i]+cumulProbs[-1]])
        for i in range(0,len(self.freqlist)):
            if ranNum < cumulProbs[i]:
                return self.statespace[i]
        return None

    def normalize(self):
        """I am creating a function to normalize my matrix before simulating transtions"""
        if sum(self.freqlist)==1:
            print "Your equilibrium frequency rates equal 1"
            print "this is your frequecy list for A,C,G and T", self.freqlist
        else:
            print "Your frequencies DO NOT EQUAL 1...better start over!"

        self.QMatrix=[[-1*(self.R[0]*self.freqlist[1]+self.R[1]*self.freqlist[2]+self.R[2]*self.freqlist[3]),self.R[0]*self.freqlist[1],self.R[1]*self.freqlist[2],self.R[2]*self.freqlist[3]],
         [self.R[0]*self.freqlist[0],-1*(self.R[0]*self.freqlist[0]+self.R[3]*self.freqlist[2]+self.R[4]*self.freqlist[3]),self.R[3]*self.freqlist[2],self.R[4]*self.freqlist[3]],
         [self.R[1]*self.freqlist[0],self.R[3]*self.freqlist[1],-1*(self.R[1]*self.freqlist[0]+self.R[3]*self.freqlist[1]+self.R[5]*self.freqlist[3]),self.R[5]*self.freqlist[3]],
         [self.R[2]*self.freqlist[0],self.R[4]*self.freqlist[1],self.R[5]*self.freqlist[2],-1*(self.R[2]*self.freqlist[0]+self.R[4]*self.freqlist[1]+self.R[5]*self.freqlist[2])]]
        #this is my lists of lists that contains my q matrix values
        print "The QMatrix that will be normalized is below. The normlized matrix is returned from this method"
        for l in self.QMatrix:
            #map(a function,a sequence)---applies your function over the entire sequence
            print ', '.join(map(str, l))
            #my Qmatrix will be displayed in a "matrix" like format(each row on it's own line)
            #now I need to calculate the weighted average of rates equation
        waverage=(self.freqlist[0]*-self.QMatrix[0][0])+(self.freqlist[1]*-self.QMatrix[1][1])+(self.freqlist[2]*-self.QMatrix[2][2])+(self.freqlist[3]*-self.QMatrix[3][3])
        print "My weighted average of rates is:" ,waverage
            #now I will use the array method of numpy so I can divide all the values of my matrix by the weighted average
        QArray=np.array(self.QMatrix)
        self.NormQMatrix=QArray/waverage

        return self.NormQMatrix
    def MarginalProbs(self):
        #(v)=e^Qv
        #import scipy
        self.MargProbs=scipy.linalg.expm(NormQMatrix*self.v)

         #self.MargProbs=[self.NormQMatrix[0][1]/-self.NormQMatrix[0][0],self.NormQMatrix[0][2]/-self.NormQMatrix[0][0],self.NormQMatrix[0][3]/-self.NormQMatrix[0][0],
                   #self.NormQMatrix[1][0]/-self.NormQMatrix[1][1],self.NormQMatrix[1][2]/-self.NormQMatrix[1][1],self.NormQMatrix[1][3]/-self.NormQMatrix[1][1],
                   #self.NormQMatrix[2][0]/-self.NormQMatrix[2][2],self.NormQMatrix[2][1]/-self.NormQMatrix[2][2],self.NormQMatrix[2][3]/-self.NormQMatrix[2][2],
                   #self.NormQMatrix[3][0]/-self.NormQMatrix[3][3],self.NormQMatrix[3][1]/-self.NormQMatrix[3][3],self.NormQMatrix[3][2]/-self.NormQMatrix[3][3]]
         return self.MargProbs
    def simulate(self):
       
        # Draw a STARTING STATE from the equilibrium frequencies
        randomNum = scipy.random.random()
        if randomNum <= self.freqlist[0]:
            self.states.extend(self.statespace[0]) #should return state a
        elif self.freqlist[0]< randomNum <=(self.freqlist[1]+self.freqlist[0]):
            self.states.extend(self.statespace[1]) #should return state c
        elif self.freqlist[3]> randomNum >self.freqlist[1]:
            self.states.extend(self.statespace[2]) #should return state g
        else:
            self.states.extend(self.statespace[3]) #should return t
        # Draw a waiting time from the appropriate exponential distribution.
        for n in range(self.NumSimulations):
            while sum(self.WaitTimes)<self.v:
                 for i in self.states[-1]:#this will pull the last element of the states list
                    if self.states[-1] ==self.statespace[0]: #we are in self.statespace a
                        waitTime1=random.expovariate(self.NormQMatrix[0][0])*-1
                        self.WaitTimes.extend([waitTime1])
                          # Draw a new state from the marginal probabilities associated with the current state.
                        newstate=self.discSamp(self.MargProbs[0],self.statespace) #draw from row 1 of list
                        if newstate==self.NormQMatrix[0][1]/-self.NormQMatrix[0][0]:
                            newstate=self.statespace[1]
                            if sum(self.WaitTimes)<self.v:
                                self.states.extend(newstate)
                        elif newstate==self.NormQMatrix[0][2]/self.NormQMatrix[0][0]:
                            newstate=self.statespace[2]
                            if sum(self.WaitTimes)<self.v:
                                self.states.extend(newstate)
                        else:
                            newstate=self.statespace[3]
                            if sum(self.WaitTimes)<self.v:
                                self.states.extend(newstate)
                          
            
                    if self.states[-1]==self.statespace[1]: #we are in statespace c
                          waitTime2=random.expovariate(self.NormQMatrix[1][1])*-1
                          self.WaitTimes.extend([waitTime2])
                          # Draw a new state from the marginal probabilities associated with the current state.
                          newstate=self.discSamp(self.MargProbs[1],self.statespace) #draw from row 2 of list
                          if newstate==self.NormQMatrix[1][0]/-self.NormQMatrix[1][1]:
                              newstate=self.statespace[0]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate)
                          elif newstate==self.NormQMatrix[1][2]/-self.NormQMatrix[1][1]:
                              newstate=self.statespace[2]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate)
                          else:
                              newstate=self.statespace[3]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate) 
        
                    if self.states[-1] == self.statespace[2]: #we are in statespace g
                          waitTime3=random.expovariate(self.NormQMatrix[2][2])*-1
                          self.WaitTimes.extend([waitTime3])
                          # Draw a new state from the marginal probabilities associated with the current state.
                          newstate=self.discSamp(self.MargProbs[2],self.statespace)
                          if newstate==self.NormQMatrix[2][0]/-self.NormQMatrix[2][2]:
                              newstate=self.statespace[0]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate)
                          elif newstate==self.NormQMatrix[2][1]/-self.NormQMatrix[2][2]:
                              newstate=self.statespace[1]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate)
                          else:
                              newstate=self.statespace[3]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate)
        
                    if self.states[-1] ==self.statespace[3]: #this will give us the wait time for "t"
                          waitTime4=random.expovariate(self.NormQMatrix[3][3])*-1
                          self.WaitTimes.extend([waitTime4])
                          
                          # Draw a new state from the marginal probabilities associated with the current state.
                          newstate=self.discSamp(self.MargProbs[3],self.statespace)
                          if newstate ==self.NormQMatrix[3][0]/-self.NormQMatrix[3][3]:
                              newstate=self.statespace[0]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate)
                          elif newstate ==self.NormQMatrix[3][1]/-self.NormQMatrix[3][3]:
                              newstate=self.statespace[1]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate)
                          else:
                              newstate=self.statespace[2]
                              if sum(self.WaitTimes)<self.v:
                                  self.states.extend(newstate)
                        
                                                           
                    return (self.states , self.WaitTimes)
            
       def IndivStates(self):
            from itertools import zip

            for i in zip(self.states):#list must be same list of zip will take length of shorter list
                self.Chains=i
            return self.Chains

            from structshape import structshape

            return structshape(self.Chains) #this should tell me how many bases I have per row and how many rows
            """
             Now I am creating a function that will determine the probability of sites between
            n number of sites. You are comparing the individual sites of the chains/sequences
            created during the continuous markov chain simulation
            # Calculate the starting likelihood for currBrl. Multiply likelihoods from
            # different sites.
            """

        def estBrl (self,currBrl=self.v, diff=0.1, thresh=0.0001): 
            
                while self.diff>self.thresh :
                    """CHANGE THE ABOVE TO BE AN ARGUMENT"""
                    likeCurr=self.MargProbs
                    vUp=currBrl+diff
                    vDown=currBrl-diff
                    if (vDown < 0):
                        vDown = 0
                    likeUp=self.MargProbs
                    likeDown=self.MargProbs
                    if likeDown>likeCurr :
                        currBrl=vDown
                        vUp=currBrl+diff
                        vDown=currBrl-diff
                    elif likeUp>likeCurr :
                        vCurr=vUp
                        vUp=currBrl+diff
                        vDown=currBrl-diff
                    else :
                        self.diff *= 0.5           
                return vCurr  
