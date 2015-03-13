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
    def __init__(self,v=0.6,statespace=None,freqlist=None,R=None,QMatrix=None,NormQMatrix=None,NumSimulations=None):
    
          if freqlist==None:
             freqlist=[]
          if statespace ==None:
              statespace=[]
            
          if R==None: #now I am going to define individual r values for transition
            #user can set these values 
            R=[]
           # R=[rAC,rAG,rAT,rCG,rCT,rGT]
          if QMatrix==None:
            self.QMatrix=[]
          if NormQMatrix==None:
              self.NormQMatrix=self.normalize()
          if NumSimulations==None:
            self.NumSimulations=input("Please enter the number of simulations:")         
                    
          self.v=v #branch length
          self.statespace=statespace #my state space of possible states
          self.NumSimulations=NumSimulations #number of simulations I will do
          self.ChainStates=[] #empty list to hold states as they change
          self.WaitingTimes=[]#empty list to hold waiting times between state changes
        #note:make sure the length of the total chain is the length of the # of wait times
          self.freqlist=freqlist
          self.R=R
         #making the values a list so I can index them
          self.QMatrix=QMatrix
          self.NormQMatrix=self.normalize()
          
          
       
    def discSamp(self,statespace=["a","t","c","g"],probs=[0.27,0.34,0.20,.19]):
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
                return statespace[i]
        return None     
          
    def normalize(self,R=[0.45,0.67,0.76,0.39,0.59,0.43],freqlist=[0.27,0.34,0.20,.19],statespace=["a","t","c","g"]):
                
        """I am creating a function to normalize my matrix before simulating transtions""" 
        import numpy as np 
        if sum(freqlist) ==1:
            print "Your equilibrium frequency rates equal 1"
            print "this is your frequecy list for A,C,G and T", freqlist
        else:
            print "Your frequencies DO NOT EQUAL 1...better start over!"
        QMatrix=[[-1*(R[0]*freqlist[1]+R[1]*freqlist[2]+R[2]*freqlist[3]),R[0]*freqlist[1],R[1]*freqlist[2],R[2]*freqlist[3]],
         [R[0]*freqlist[0],-1*(R[0]*freqlist[0]+R[3]*freqlist[2]+R[4]*freqlist[3]),R[3]*freqlist[2],R[4]*freqlist[3]],
         [R[1]*freqlist[0],R[3]*freqlist[1],-1*(R[1]*freqlist[0]+R[3]*freqlist[1]+R[5]*freqlist[3]),R[5]*freqlist[3]],
         [R[2]*freqlist[0],R[4]*freqlist[1],R[5]*freqlist[2],-1*(R[2]*freqlist[0]+R[4]*freqlist[1]+R[5]*freqlist[2])]]
        #this is my lists of lists that contains my q matrix values
        #for l in self.QMatrix:
            #map(a function,a sequence)---applies your function over the entire sequence 
            #print ', '.join(map(str, l))
            #my Qmatrix will be displayed in a "matrix" like format(each row on it's own line)
            #now I need to calculate the weighted average of rates equation
        waverage=(freqlist[0]*-QMatrix[0][0])+(freqlist[1]*-QMatrix[1][1])+(freqlist[2]*-QMatrix[2][2])+(freqlist[3]*-QMatrix[3][3])
        print "My weighted average of rates is:" ,waverage 
            #now I will use the array method of numpy so I can divide all the values of my matrix by the weighted average
        QArray=np.array(QMatrix)
        NormQMatrix=QArray/waverage 
        MargProbs=[NormQMatrix[0][1]/-NormQMatrix[0][0],NormQMatrix[0][2]/-NormQMatrix[0][0],NormQMatrix[0][3]/-NormQMatrix[0][0],
                   NormQMatrix[1][0]/-NormQMatrix[1][1],NormQMatrix[1][2]/-NormQMatrix[1][1],NormQMatrix[1][3]/-NormQMatrix[1][1],
                   NormQMatrix[2][0]/-NormQMatrix[2][2],NormQMatrix[2][1]/-NormQMatrix[2][2],NormQMatrix[2][3]/-NormQMatrix[2][2],
                   NormQMatrix[3][0]/-NormQMatrix[3][3],NormQMatrix[3][1]/-NormQMatrix[3][3],NormQMatrix[3][2]/-NormQMatrix[3][3]]

        return NormQMatrix, MargProbs
        

    def simulate(self,freqlist=[0.27,0.34,0.20,.19],statespace=["a","t","c","g"],v=0.06):
        import scipy
        import random
        #create a list of the marginal probabilties associated with each nucleotides
        states=[] #this is a list to hold that states of the chain 
        waitTimes=[]
        # Draw a STARTING STATE from the equilibrium frequencies              
        randomNum = scipy.random.random()
        if randomNum <= freqlist[0]: 
            states.extend(statespace[0]) #should return state a
        elif freqlist[0]< randomNum <=(freqlist[1]+freqlist[0]): 
            states.extend(statespace[1]) #should return state c
        elif freqlist[3]> randomNum >freqlist[1]:
            states.extend(statespace[2]) #should return state g
        else:
            states.extend(statespace[3]) #should return t
        # Draw a waiting time from the appropriate exponential distribution.
        while sum(waitTimes)<v:
            for i in states[-1]: #this will pull the last element of the states list
              if states[-1] ==statespace[0]: #we are in statespace a
                  waitTime=random.expovariate(self.NormQMatrix[0][0])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[0]) #draw from row 1 of list
                  if newstate==self.NormQMatrix[0][1]/-self.NormQMatrix[0][0]:
                      newstate=statespace[1]
                      states.extend(newstate)
                  elif newstate==self.NormQMatrix[0][2]/self.NormQMatrix[0][0]:
                      newstate=statespace[2]
                      states.extend(newstate)
                  else:
                      newstate=statespace[3]
                      states.extend(newstate)   
                      
              elif states[-1]==statespace[1]: #we are in statespace c
                  waitTime=random.expovariate(-self.NormQMatrix[1][1])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[1]) #draw from row 2 of list 
                  if newstate==self.NormQMatrix[1][0]/-self.NormQMatrix[1][1]:
                      newstate=statespace[0]
                      states.extend(newstate)
                  elif newstate==self.NormQMatrix[1][2]/-self.NormQMatrix[1][1]:
                      newstate=statespace[2]
                      states.extend(newstate)
                  else:
                      newstate=statespace[3]
                      states.extend(newstate)
                      
              elif states[-1] == statespace[2]: #we are in statespace g
                  waitTime=random.expovariate(-self.NormQMatrix[2][2])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[2])
                  if newstate==self.NormQMatrix[2][0]/-self.NormQMatrix[2][2]:
                      newstate=statespace[0]
                      states.extend(newstate)
                  elif newstate==self.NormQMatrix[2][1]/-self.NormQMatrix[2][2]:
                      newstate=statespace[1]
                      states.extend(newstate)
                  else:
                      newstate=statespace[3]
                      states.extend(newstate)
                      
              else: #this will give us the wait time for "t"
                  waitTime=random.expovariate(-self.NormQMatrix[3][3])
                  waitTimes.extend(waitTime) 
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[3])
                  if newstate ==self.NormQMatrix[3][0]/-self.NormQMatrix[3][3]:
                      newstate=statespace[0]
                      states.extend(newstate)
                  elif newstate ==self.NormQMatrix[3][1]/-self.NormQMatrix[3][3]:
                      newstate=statespace[1]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[2]
                      states.extend(newstate)
            return states
            return waitTimes
