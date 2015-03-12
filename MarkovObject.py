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
    def __init__(self,v=0.6,statespace=None,freqlist=None,R=None,QMatrix=None,NumSimulations=None,NormQMatrix=None):
    
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
              self.NormQMatrix=[]
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
          self.NormQMatrix=NormQMatrix
          
    def normalize(self,NormQMatrix=None,R=[0.45,0.67,0.76,0.39,0.59,0.43],freqlist=[0.27,0.34,0.20,.19]):
                
        """I am creating a function to normalize my matrix before simulating transtions""" 
        
        if sum(self.freqlist) ==1:
            print "Your equilibrium frequency rates equal 1"
            print "this is your frequecy list for A,C,G and T", self.freqlist
        else:
            print "Your frequencies DO NOT EQUAL 1...better start over!"
        #this is my lists of lists that contains my q matrix values
        self.QMatrix=[[-1*(R[0]*self.freqlist[1]+R[1]*self.freqlist[2]+R[2]*self.freqlist[3]),R[0]*self.freqlist[1],R[1]*self.freqlist[2],R[2]*self.freqlist[3]],
                      [R[0]*self.freqlist[0],-1*(R[0]*self.freqlist[0]+R[3]*self.freqlist[2]+R[4]*self.freqlist[3]),R[3]*self.freqlist[2],R[4]*self.freqlist[3]],
                      [R[1]*self.freqlist[0],R[3]*self.freqlist[1],-1*(R[1]*self.freqlist[0]+R[3]*self.freqlist[1]+R[5]*self.freqlist[3]),R[5]*self.freqlist[3]],
                      [R[2]*self.freqlist[0],R[4]*self.freqlist[1],R[5]*self.freqlist[2],-1*(R[2]*self.freqlist[0]+R[4]*self.freqlist[1]+R[5]*self.freqlist[2])]] 
                         
        for l in self.QMatrix:
            #map(a function,a sequence)---applies your function over the entire sequence 
            print ', '.join(map(str, l))
            #my Qmatrix will be displayed in a "matrix" like format(each row on it's own line)
            #now I need to calculate the weighted average of rates equation
            waverage=(self.freqlist[0]*-self.QMatrix[0][0])+(self.freqlist[1]*-self.QMatrix[1][1])+(self.freqlist[2]*-self.QMatrix[2][2])+(self.freqlist[3]*-self.QMatrix[3][3])
            print "My weighted average of rates is:" ,waverage 
            
            import numpy as np    
            if self.NormQMatrix==None:
            #now I will use the array method of numpy so I can divide all the values of my matrix by the weighted average
                QArray=np.array(self.QMatrix)
                self.NormQMatrix=QArray/waverage 
        return self.NormQMatrix
    def discSamp(self,statespace=("a","t","c","g"),probs=[0.27,0.34,0.20,.19]):
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
    def simulate(self):
        import scipy
        import random
        #create a list of the marginal probabilties associated with each nucleotide
        self.MargProbs=[self.NormQMatrix[0][1]/-self.NormQMatrix[0][0],self.NormQMatrix[0][2]/-self.NormQMatrix[0][0],self.NormQMatrix[0][3]/-self.NormQMatrix[0][0],
                   self.NormQMatrix[1][0]/-self.NormQMatrix[1][1],self.NormQMatrix[1][2]/-self.NormQMatrix[1][1],self.NormQMatrix[1][3]/-self.NormQMatrix[1][1],
                   self.NormQMatrix[2][0]/-self.NormQMatrix[2][2],self.NormQMatrix[2][1]/-self.NormQMatrix[2][2],self.NormQMatrix[2][3]/-self.NormQMatrix[2][2],
                   self.NormQMatrix[3][0]/-self.NormQMatrix[3][3],self.NormQMatrix[3][1]/-self.NormQMatrix[3][3],self.NormQMatrix[3][2]/-self.NormQMatrix[3][3]]
        states=[] #this is a list to hold that states of the chain 
        waitTimes=[]
        # Draw a STARTING STATE from the equilibrium frequencies              
        randomNum = scipy.random.random()
        if randomNum <= self.freqlist[0]: 
            states.extend(self.statespace[0]) #should return state a
        elif self.freqlist[0]< randomNum <=sum(self.freqlist[1]+self.freqlist[0]): 
            states.extend(self.statespace[1]) #should return state c
        elif self.freqlist[3]> randomNum >self.freqlist[1]:
            states.extend( self.statespace[2]) #should return state g
        else:
            states.extend(self.statespace[3]) #should return t
        # Draw a waiting time from the appropriate exponential distribution.
        while sum(waitTimes)<self.v:
            for i in states[-1]: #this will pull the last element of the states list
              if states[-1] == self.statespace[0]: #we are in statespace a
                  waitTime=random.expovariate(-self.NormQMatrix[0][0])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[0]) #draw from row 1 of list
                  if newstate==self.NormQMatrix[0][1]/-self.NormQMatrix[0][0]:
                      newstate=self.statespace[1]
                      states.extend(newstate)
                  elif newstate==self.NormQMatrix[0][2]/-self.NormQMatrix[0][0]:
                      newstate=self.statespace[2]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[3]
                      states.extend(newstate)   
                      
              elif states[-1] ==self.statespace[1]: #we are in statespace c
                  waitTime=random.expovariate(-self.NormQMatrix[1][1])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[1]) #draw from row 2 of list 
                  if newstate==self.NormQMatrix[1][0]/-self.NormQMatrix[1][1]:
                      newstate=self.statespace[0]
                      states.extend(newstate)
                  elif newstate==self.NormQMatrix[1][2]/-self.NormQMatrix[1][1]:
                      newstate=self.statespace[2]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[3]
                      states.extend(newstate)
                      
              elif states[-1] == self.statespace[2]: #we are in statespace g
                  waitTime=random.expovariate(-self.NormQMatrix[2][2])
                  waitTimes.extend(waitTime)
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[2])
                  if newstate==self.NormQMatrix[2][0]/-self.NormQMatrix[2][2]:
                      newstate=self.statespace[0]
                      states.extend(newstate)
                  elif newstate==self.NormQMatrix[2][1]/-self.NormQMatrix[2][2]:
                      newstate=self.statespace[1]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[3]
                      states.extend(newstate)
                      
              else: #this will give us the wait time for "t"
                  waitTime=random.expovariate(-self.NormQMatrix[3][3])
                  waitTimes.extend(waitTime) 
                  # Draw a new state from the marginal probabilities associated with the current state.
                  newstate=self.discSamp(self.MargProbs[3])
                  if newstate ==self.NormQMatrix[3][0]/-self.NormQMatrix[3][3]:
                      newstate=self.statespace[0]
                      states.extend(newstate)
                  elif newstate ==self.NormQMatrix[3][1]/-self.NormQMatrix[3][3]:
                      newstate=self.statespace[1]
                      states.extend(newstate)
                  else:
                      newstate=self.statespace[2]
                      states.extend(newstate)
chain1=ContinMarkov()
chain1.simulate()      
   
   """           
       
        def NumSimulations(self):
            for i in range(self.NumSimulations): 
                simulate()   
            # Keep drawing waiting times and new states until the total time exceeds
            # the branch length.
            #Now I need estimate current likelihood of branchlength 
              
        def ProbMatrix():
            First I need to create a matrix of probabilities from the normalized QMatrix
            (v)=e^Qv
            import scipy
            ProbM=scipy.linalg.expm(NormQMatrix*self.v)
            
            return ProbM
            
        def IndivStates(self,Chains):
            from itertools import *

            for i in izip([1, 2, 3], ['a', 'b', 'c']):
                print i
            
            (1, 'a')
            (2, 'b')
            (3, 'c')
            
                  
            try:
                from itertools import izip as zip
            except ImportError: # will be 3.x series
                pass
                    
            for i in zip(states):#list must be same list of zip will take length of shorter list
                self.Chains=i
            return self.Chains
            
            from structshape import structshape
            
            return structshape(self.Chains) #this should tell me how many bases I have per row and how many rows
   """        
c1=ContinMarkov()

           
"""
for i, var in enumerate(self.Chains):
if i == len(self.Chains) - 1:
return 'last element:', var

            
  
Now I need to create something that goes through each row counts the number of values
(which would depend on the number of simulations) and associates their self.statespace value 
as it moves to the next index in the row with it's respective self.statespace value"""
"""
def estBrl(self,currBrl,diff,thresh):
 Now I am creating a function that will determine the probability of sites between 
n number of sites. You are comparing the individual sites of the chains/sequences
created during the continuous markov chain simulation
  # Calculate the starting likelihood for currBrl. Multiply likelihoods from
# different sites. 
   


Currlike=1
Currlike *=likelihood
		  
Downlike=itBrl-diff
Uplikel=
			Downlikelihood=likeCalc for initBrl-difference:  make sure this doesn't go below zero ;must be >0
			Uplikelihood= likeCalc for initBrl+difference 
			(check to see which direction is better...going lower or higher)
			(make sure you don't get a negative branch length)
			
			Assign your branch length(re set value) to the largest value you get. 
			This will be the initBrl
			
			initBrl=max( currentlikelihood,downliklihood,uplikelihood)
			
			*Keep doing this until your current is not larger or smaller than your down or up values
			You need a while loop;
			
   
# Calculate the likelihood for upBrl = currBrl+diff
# Calculate the likelihood for downBrl = currBrl-diff. Check to make sure
# this doesn't go below 0
While(diff > threshhold):
				If (upLike >current likelihood):
					***you can use recursion here 
					initiBrl=estBrl(upbrl,diff, Thresh) …..use this if the up value is better and it will keep calling until it finds the best branch length
				elif(downlike>current likelihood):
					initBrl=estBrl(downbrl,diff,thresh) 
				**make sure you change the difference to make sure it keeps incrementing 
					Else:
					Diff *=0.5 
					Currlikeli=likelicalc(initBrl) **now we have a new best estimate and then you can start over to calculate the down and up like likelihood branch lengths 
					Calcnew down Brl & likelihood 
					Calcl new up Brl & likelihood
					
		Return initBrl ...this returns your best estimate of branch length 
# NOTE: The ML branch-length estimate can sometimes be infinity. To avoid
# the function trying to reach infinity, you should add a cutoff value for
# the maximum possible branch-length. Something like 20 should be plenty
# large.
if #upBrl likelihood > currBrl likelihood :
currBrl = upBrl
elif #downBrl likelihood > currBrl likelihood :
currBrl = downBrl
else:
diff *= 0.5
# Recalculate likelihoods for new currBrl, upBrl, and downBrl
return currBrl
"""
"""
The continuous markov chain class and methods are defined above. They are used and tested below.
"""
# Simulate one site. Estimate branch length.
# Simulate several hundred sites. Estimate branch length.
# Run the above simulations repeatedly and examine variation in the estimated branch lengths.



