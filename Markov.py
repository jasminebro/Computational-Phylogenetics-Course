"""
Exercise 4
Discrete-time Markov chains
@author: jasmine brown 
"""
"""
In this exercise, we will explore Markov chains that have discrete state spaces
and occur in discrete time steps. To set up a Markov chain, we first need to
define the states that the chain can take over time, known as its state space.
To start, let's restrict ourselves to the case where our chain takes only two
states. We'll call them A and B.
"""
# Create a tuple that contains the names of the chain's states

Chain=("A","B") 
"""
The behavior of the chain with respect to these states will be determined by
the probabilities of taking state A or B, given that the chain is currently in
A and B. Remember that these are called conditional probabilities (e.g., the
probability of going to B, given that the chain is currently in state A is
P(B|A).)
We record all of these probabilities in a transition matrix. Each row
of the matrix records the conditional probabilities of moving to the other
states, given that we're in the state associated with that row. In our example
row 1 will be A and row 2 will be B. So, row 1, column 1 is P(A|A); row 1,
column 2 is P(B|A); row 2, column 1 is P(A|B); and row 2, column 2 is P(B|B).
All of the probabilities in a ROW need to sum to 1 (i.e., the total probability
associated with all possibilities for the next step must sum to 1, conditional
on the chain's current state).
In Python, we often store matrices as "lists of lists". So, one list will be
the container for the whole matrix and each element of that list will be
another list corresponding to a row, like this: mat = [[r1c1,r1c2],[r2c1,r2c2]].
We can then access individual elements use two indices in a row. For instance,
mat[0][0] would return r1c1. Using just one index returns the whole row, like
this: mat[0] would return [r1c1,r1c2].
Define a transition matrix for your chain below. For now, keep the probabilties
moderate (between 0.2 and 0.8).
"""
# Define a transition probability matrix for the chain with states A and B
Matrix=[[0.6,0.4], [0.3,0.7]]
"""
I needed a visual to help me see what I was really doing 

     A              B

A   0.6            0.4


B   0.3            0.7

"""
# Try accessing a individual element or an individual row
# Element
#I'm going to grab the element from row 2 and the first element
print Matrix[1][0]
# Row
#lets grab the second row 
print Matrix[1]
"""
Now, write a function that simulates the behavior of this chain over n time
steps. To do this, you'll need to return to our earlier exercise on drawing
values from a discrete distribution. You'll need to be able to draw a random
number between 0 and 1 (built in to scipy), then use your discrete sampling
function to draw one of your states based on this random number.
"""
# Import scipy U(0,1) random number generator

# Paste or import your discrete sampling function

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
	
# Write your Markov chain simulator below. Record the states of your chain in
# a list. Draw a random state to initiate the chain.
# Run a simulation of 10 steps and print the output.

#I tried to do it as the instructions went step by step but this method seemed easier. Let me know if I should've done it another way. :/

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
	
	
print Markov(Chain,Matrix,10)
	
# ----> Try to finish the above lines before Tues, Feb. 10th <----
# Now try running 100 simulations of 100 steps each. How often does the chain
# end in each state? How does this change as you change the transition matrix?

print Markov(Chain,Matrix,100)

# Try defining a state space for nucleotides: A, C, G, and T. Now define a
# transition matrix with equal probabilities of change between states.


nucleotides =("a","t","g","c") #this is going  to be my state space 

NucMatrix=[[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25]]

statelist=[]
for i in range (100):
    a=Markov(nucleotides,NucMatrix,100)
    state=a[len(a)-1]
    statelist.append(state)
print statelist
    
A=statelist.count("a")    
T=statelist.count("t")
G=statelist.count("g")
C=statelist.count("c")

print(A,T,G,C)

# Again, run 100 simulations of 100 steps and look at the ending states. Then
# try changing the transition matrix.

nucleotides2 =("a","t","g","c") #this is going  to be my state space 
#I'll just use random state numbers and see what happrens 
NucMatrix2=[[0.20,0.40,0.30,0.10],[0.40,0.30,0.20,0.10],[0.25,0.25,0.25,0.25],[0.50,0.20,0.10,0.20]]

statelist2=[]
for i in range (100):
    a=Markov(nucleotides2,NucMatrix2,100)
    state=a[len(a)-1] #tried this to see if it would make a difference as far as the letters I got and to make it more general
    statelist2.append(state)
print statelist2
    
A2=statelist2.count("a")    
T2=statelist2.count("t")
G2=statelist2.count("g")
C2=statelist2.count("c")

print(A2,T2,G2,C2)
  
