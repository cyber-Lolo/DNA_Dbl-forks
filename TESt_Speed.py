#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 19:03:40 2022

@author: lopid
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 13:03:51 2022

@author: lopid
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 13:47:28 2022

@author: laurianepierrotdeseilligny
"""

#############################################################################

"""
Libraries
"""

import numpy as np
import random
import matplotlib.pyplot as plt

#############################################################################
"""
DNA 
legend : b - hole  * island
Exemplen(origin is arbitrary because of circurlarity, we set it on an island)

**----***--***--****----

I   2  3  .. 
H   4  2  ...

"""

"""
Function that gives the cummulated distances 
"""

def cummu(A): #A is a l*3 matrice
	l = len(A)
	A[0,2]=0
	for i in range(1,l):
		A[i,2] = np.sum(A[0:i,1])+np.sum(A[0:i,0])
	return(A)

#############################################################################

# DNA is a slope with n sites (island size, hole size, cummulated hole size)

#Column 0 : Islands
#Column 1 : Holes sizes
#Column 2 : Cummulated holes' sizes


# For convenience we will consider it to be a circle : the n+1 th site corresponds to the 1st site

# Therefore, N_sites newly activated origins lead to N_sites holes

# At the very first step of replication, the sum of holes's size = DNA length

#############################################################################

"""
Nucleation Algorithm
"""

def DoubleL(
	 dt,  # time step 
	 v,   # fork veolocity (speed of dilatation)
	 alpha, # poisson law parameter->  
	 L_DNA, # size of dna 
	 N # number of initial sites
	):
	
	
#### Initialization ####
	Nme=10

	#At time zero, nucleation
	t=0
	#Number of activated sites at firsts step is chosen
	#Nucleation stars with the activation of at least 1 site
	#The first site serve as genome positon reference
	#N_s=np.array([ 41, 425, 529, 591, 634, 767, 873, 954])


	N_s = np.array(random.sample(range(1,L_DNA),N))
	N_s = np.sort(N_s)
	#N_s=np.array([1,20,25,27,50,70,80,97])
	#We want the distance between holes
	L1 = np.insert(N_s, 0, 0) 
	L2 = np.append(N_s,L_DNA)
	N_s = L2 - L1


	To_a = np.zeros((len(N_s),2))
	To_a[:,1] = N_s
	
	l = len(N_s)
	Db_lst = np.zeros((l,3))
	
	Db_lst[:,0:2] = To_a
	
	Db_lst = cummu(Db_lst)
	Db_lst =np.delete(Db_lst, 0,axis=0)
	Dl=np.shape(Db_lst)[0]
	# Island[:,0] coordinates
	# Island[:,1] time
	Island = t*np.ones((Dl,2))
	#At initialization, island positions are the ones of cummulated holes size
	#And the firts island is placed at zero
	X = Db_lst[:,2]
	Island[:,0] = X
	

	
	# holes[:,0] coordinates
	# holes[:,1] time
	holes = np.zeros((1,2))

	S=N
	P=0
	Copy = np.array(Db_lst)
#### Nulceation ####

	while np.sum(Copy[:,1])>0:
		#print([0,0,0,0],np.sum(Db_lst[:,1]))
		t = t+dt 
		I = alpha*t

		
		 #rate of newly activated sites at time t per unit of time and space
		V=0.75*np.exp(-Db_lst[:,1]/L_DNA)

		Db_lst[:,0]=Db_lst[:,0]+2*V*dt  #islands grow
		Db_lst[:,1]=Db_lst[:,1]-2*V*dt  #holes shrink
		
		
		Db_lst[:,1]=Db_lst[:,1]*(Db_lst[:,1]>0) #ReLu on holes : if <0, then becomes 0
		

		Db_lst = cummu(Db_lst) # Updating the column of cummulated holes' sizes
		
		Copy = np.array(Db_lst) #Copy 

		#A = np.zeros((1,3)) #empty the double list so that it will carry the new values

		H = Copy[:,1] #holes'sizes
		
		#Tracking holes that are 0
		Coal_ind = np.nonzero(H==0)[0] #Index of holes that are disappearing (coalescence)
		
		
		P=P+np.sum(H==0)
		m = len(Coal_ind) #number of points that coalesce
		ho = t*np.ones((m,2)) # Preparing to get their genome position
		
#### Coalescence ####

		for j in range(m): #for all the coalescing holes
			

			c = Coal_ind[j] #get their index
			#compute their genome position
			
			#fusion with the adjacent island (on the right)

			if c <= len(Copy)-2:  #if coalescence does not happen on the last point
				#print(Copy[c+1,:],Copy[c,:])
				ho[j,0] = Copy[c,2]+Copy[c,0]*0.5#-Copy[0,0]*0.5 #compute their genome position
				Copy[c+1,0] = Copy[c,0]+Copy[c+1,0]
				#print(c,Copy[c+1,:])
			elif c == len(Copy)-1:
				Copy[c,0]=L_DNA-Copy[c,2]
				ho[j,0]=L_DNA
			
		holes = np.append(holes , ho, axis=0)	  
		Coal_ind_del=Coal_ind[Coal_ind <=len(Copy)-2]
		#Delete the holes that have disappeared		
		
		Copy = np.delete(Copy,Coal_ind_del,axis=0)
		
		#print(Copy)
		H=Copy[:,1]
		
		n = np.shape(Copy)[0]	# number of (island, holes, cummulated holes)
		Copy=cummu(Copy)
		Db_lst=Copy
	return(Island,holes,Db_lst,S,P)
	"""	
#### New sites coordinates ####
	
		Db_lst = np.ones((1,3))
	
		for i in range(n) :
			
			h = Copy[i,1] # hole number i
			N_sites = I*dt*h  #number of new sites to be activated
			
			
			if N_sites<1 or int(N_sites)>int(h):  #if no new sites, Db_list appends the unchanged island and holes
				x = np.array([Copy[i,:]])
				Db_lst = np.append(Db_lst,x, axis=0)

			if N_sites>=1: #if there is at least 1 new size
				S=S+int(N_sites)
				#Place them randomly on the current hole
				
				New_sites = np.array(random.sample(range(1,int(h)),int(N_sites)))
				New_sites = np.sort(New_sites)
				
				L1 = np.insert(New_sites, 0, 0) 
				L2 = np.append(New_sites,h)
				#We want the distance between the holes
				
				New_sites = L2 - L1
				#print(np.sum(New_sites),h)
				#New  sites means islands of length zero				
				To_add = np.zeros((len(New_sites),3))
				
				#The first island of the new sub ensemble is the one corresponding to
				#the hole that is being split
				To_add[0,0]=Copy[i,0]
				
				#Add the new holes
				To_add[:,1] = New_sites
				Db_lst = np.append(Db_lst, To_add, axis=0)
			
				o = len(New_sites)
				isl = t*np.ones((o,2))
				isl[:,0]=Copy[i,2]-h+New_sites
				Island = np.append(Island, isl,axis=0)
		
		# delete first line that was null
		
		Db_lst = np.delete(Db_lst,0, axis=0)
		
		#Update the cummulated holes of the newly added sites
		Db_lst = cummu(Db_lst)
	return(Island,holes,Db_lst,S,P)
	
		#print(Db_lst[np.shape(Copy)[0]-1,2])

		Or = Db_lst[:,0]
		Or_ind = np.nonzero(Or==0)[0]
		o = len(Or_ind)
		isl = t*np.ones((o,2))

#### New sites coordinates ####
		
		for k in range(o):
			s = Or_ind[k]
			isl[k,0] = Db_lst[s,2]#-Db_lst[0,0]*0.5
		Island = np.append(Island, isl,axis=0)
	print(Copy[len(Copy)-1,2])
	return(Island,holes,Db_lst,S,P)


#############################################################################


Plotting
"""

if __name__=="__main__":
	
#### Variables' choice ####

	K=0
	#total time
	t = 0
	
	#time step
	dt = 0.1
	
	#growth speed
	v = 0.5
	
	#DNA length
	L_DNA = 1000
	
	#Rate
	alpha = 10**-5
	
	#at t=0 DNA starts to replicate
	t = t+dt
	
	# Number of activated sites 
	N_sites = 10
	
	Temps_moy=[]
	Fluc=[]
	for i in range (1,110,5):
		Temps=[]
		moy=0
		std=0
		for j in range (50):
			Island,holes,Db_lst,T,R = DoubleL(dt, v, alpha, L_DNA, i)
			M=np.max(holes[:,1])
			holes = np.delete(holes,0, axis=0)
			M=np.max(holes[:,1])
			Temps.append(M)
		moy=np.mean(np.array(Temps))
		std=np.std(np.array(Temps))
		Temps_moy.append(moy)
		Fluc.append(std)
	pass
#############################################################################
"""
	Island,holes,Db_lst,T,R = DoubleL(dt, v, alpha, L_DNA, N_sites)
	Max_ind=np.nonzero(Island[:,0]==np.max(Island[:,0]))

	holes = np.delete(holes,0, axis=0)
	A=plt.scatter(Island[:,0],Island[:,1])
	B=plt.scatter(holes[:,0],holes[:,1])
	plt.plot([0,70],[0,70/v])
	plt.xlabel('Genome position')
	plt.ylabel('Time')
	plt.legend((A,B),('Activated Islands', 'Coalescence'),scatterpoints=1,loc='upper right',ncol=1,fontsize=8)
	M=np.max(holes[:,1])
"""
	

	
