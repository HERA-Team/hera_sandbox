from omnical import arrayinfo as arr
import numpy as np

def sort_pol_reds(antpos,num_pols=2,reds_tol=0.1):
	"""
	this sorts antpos into:
	XY/YX pairs (which we want to contend are equal by appending the flip of each tuple)
	XX/YY pairs (we need to sort out which is which)
	
	
	XXX need to add support for +x+x+ orientations
	XXX need to switch to tolerance method rather than exact redundancy
	XXX need to switch to tolerance method for third dimension (I assume perfect flatness and offset polarizations by value 'boost_factor' in the z-direction)
	
	"""
	if num_pols>2:
		raise NotImplementedError('Not ready for >2 orientations yet.')
		##TODO: add support for arbitraty number of orientations##
	boost_factor=100
	boost = antpos+[0,0,boost_factor]#Ypol is +100 above Xpol in the z direction
	antpos2 = np.append(antpos,boost,axis=0) #X and Y
	
	#call the vanilla omnical on the not-vanilla (chocolate?) antpos array
	reds = arr.compute_reds(antpos2,tol=reds_tol)

	crosspol_reds, linpol_reds = [],[]
	
	#classify whether each group is crosspol or linpol by checking third dimension
	#a little fragile since it assumes an array on a perfectly flat surface
	
	#TODO: switch to a tolerance method
	for group in reds:
		tuple0 = group[0]
		i,j = tuple0[0],tuple0[1]
		#XXX is this a little cavalier? vv
		if i==j+len(antpos) or j==i+len(antpos): continue #gets rid of autocorrelations
		#XXX ^^
		if antpos2[i][2] - antpos2[j][2] != 0: crosspol_reds.append(group) #check 3rd dimension
		else: linpol_reds.append(group)

	#impose XY/YX redundancy by flipping each of these ij tuples
	#there must be a more clever way to do this
	for i,group in enumerate(crosspol_reds):
		temp = []
		for tuple in group: temp.append(tuple[::-1]) 
		#^^ if I appended straight to "group" it would loop forever
		for tuple in temp: group.append(tuple)
		crosspol_reds[i] = group
	
	#need to sort out which one is X and which is Y in the non-pol case
	XX_reds, YY_reds = [], []
	for i, group in enumerate(linpol_reds):
		temp_group_XX,temp_group_YY = [],[]
		for tuple in group:
			i,j = tuple
			if antpos2[i][2]==0: temp_group_XX.append(tuple)
			elif antpos2[i][2]==boost_factor: temp_group_YY.append(tuple)
			else:
				raise ValueError('Something has gone wrong with number of orientations')
				
				#TODO as above, this is a little fragile as stands. Should shift to a tolerance
				
		if len(temp_group_XX) < 2: continue #it's not actually redundant with anything
		XX_reds.append(temp_group_XX)
		YY_reds.append(temp_group_YY)
	
	return {'xx':XX_reds,'yy':YY_reds,'crosspol':crosspol_reds}

def antnum_nopol2pol(antnum,antpos,num_pols=2,return_pos=False):
	"""
	maps the output of sort_pol_reds back to physical antenna number 
	(on a per antenna basis)
	"""
	antpos = np.array(antpos)
	N = antpos.shape[0]
	
	if antnum > num_pols*N: raise ValueError('Invalid antenna index. No mapping for this number of polarizations.')
	
	if num_pols==1: num = antnum
	elif num_pols==2:
		if antnum<N: num = antnum
		if antnum>=N: num = antnum-N
	else:
		raise NotImplementedError('Not ready > 2 orientations yet.')
	
	if not return_pos: return num
	else: return [num,antpos[num]]

def red_dict_nopol2pol(red_dict,antpos,num_pols=2,return_pos=False):
	"""
	maps the output of sort_pol_reds back to physical antenna number (in batch)
	
	XXX all these for loops are inefficient and should be vectorized
	"""
	for key in red_dict.keys():
		for i,red in enumerate(red_dict[key]):
			for j,tuple in enumerate(red):
				#tuples are immutable, lists are not
				red_dict[key][i][j] = list(red_dict[key][i][j])
				for k,num in enumerate(tuple):
					red_dict[key][i][j][k] = antnum_nopol2pol(num,antpos,num_pols=num_pols,return_pos=return_pos)	
	return red_dict


def red_dict_2polstring(red_dict):
	"""
	append polarization keys directly to red_dict object which requires a parsing to string-type. Don't know how useful this is in the long run but it's good for checking on things.
	"""
	for key in red_dict.keys():
		if key != 'crosspol':
			pol1=key[0]
			pol2=key[1]
		else:
			pol1='x'
			pol2='y'
		for i,red in enumerate(red_dict[key]):
			for j,tuple in enumerate(red):
				"""
				tuples are immutable, lists are not
				"double" parsing to list (ie once here and once in red_dict_nopol2pol) does 
				not affect functionality and allows for backwards functionality
				"""
				red_dict[key][i][j] = list(red_dict[key][i][j])
			 
				b = str(red_dict[key][i][j][0])
				c = str(red_dict[key][i][j][1])
			
				red_dict[key][i][j] = list((b+pol1,c+pol2))
	return red_dict

			
			
			
				


############ TESTING ############
antpos =  np.array([[0,0,0],[0,1,0],[0,2,0],[1,0,0],[1,1,0],[1,2,0]])

red_dict = sort_pol_reds(antpos)

print 'ANTENNAE | POSITION'
for a, pos in enumerate(antpos):
	print a, pos

print ''
print 'RED_DICT (sort_pol_reds output)'
for key in red_dict.keys(): 
	print key,red_dict[key]
	print ''

print ''
print 'RED_DICT (mapped back to physical antennae)'
print red_dict_nopol2pol(red_dict,antpos)

print ''
print 'RED_DICT (with pol keys)'
R = red_dict_2polstring(red_dict)
for key in R.keys(): 
	print key,R[key]
	print ''
####################################
	
