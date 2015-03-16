class ContinMarkov(object):
    def __init__(self,v=0.6,MargProbs=None,statespace=["a","t","c","g"],freqlist=[0.27,0.34,0.20,0.19],R=[0.45,0.67,0.76,0.39,0.59,0.43],QMatrix=None,NormQMatrix=None,NumSimulations=1):
        if QMatrix==None:
            self.QMatrix=[]
        if NormQMatrix==None:
            self.NormQMatrix=self.normalize() 
        if MargProbs==None:
            self.MargProbs=self.MargProbs()
        self.v=v #branch length
        self.statespace=statespace #my state space of possible states
        self.NumSimulations=NumSimulations #number of simulations I will do
        self.states=[] #empty list to hold states as they change
        self.WaitingTimes=[]#empty list to hold waiting times between state changes
        #note:make sure the length of the total chain is the length of the # of wait times
        self.freqlist=freqlist
        self.R=R
        
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
