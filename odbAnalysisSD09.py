#These scripts are for analysing output database stresses in the rail
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random	
from odbAccess import *	
from  abaqus import session
import numpy as np
nt=20 # number of ties
nl=nt-1 #number of load positions
ki=[0+190*i for i in range(10)]
for p in range(10):
	h=open('SD_stresses_'+str(p)+'.txt','w')
	s=open('SD_Stiff_'+str(p)+'.txt','w')
	pi=open('SD_Poiss_'+str(p)+'.txt','w')
	o=open('stiff_av_'+str(p)+'.txt','w')
	t=open('stresses_av_'+str(p)+'','w')
	q1=open('Concrete_E_'+str(p)+'.txt','w')
	r1=open('Concrete_P_'+str(p)+'.txt','w')
	SD_s33_list=[]
	Seff=open('Eff_Stress'+str(p)+'.txt','w')
	Meff=open('Mean_Stress'+str(p)+'.txt','w')
	#kf=189
	for k in range(1):
	# for k in [0,1]:
		l0=ki[(p)]+nl*k
		lf=ki[(p)]+nl*(k+1)
		#Create a file that contain maximum stresses for every every standard deviation given to random stiffness in the ties
		g=open('RailMaxStress_'+str(p)+'_'+str(k)+'.txt','w')
		gp=open('RailMaxPrinc_'+str(p)+''+str(k)+'.txt','w')
		# q1=open('Concrete_E'+str(k)+'.txt','w')
		# r1=open('Concrete_P'+str(k)+'.txt','w')
		maxs33=[]
		MDirs=[]
		MPrs=[]
		n=open('RailMaxDisp_'+str(p)+''+str(k)+'.txt','w')
		maxs33=[]
		MaxDisp=[]
		for j in range(l0,lf):
		# lf1=l0+2
		# for j in range(l0,lf1):
			# session.upgradeOdb("G:/EDissertProject/odbs04-19-14/JobAuto-"+str(j)+".odb", 
				# "G:/EDissertProject/odbfiles-1/JobAuto-"+str(j)+"-1.odb",)
			#Open output database model
			odb=openOdb('JobAuto-'+str(j)+'.odb')
			#Stresses:
			stresses=odb.steps['Step-2'].frames[-1].fieldOutputs['S']
			rail=odb.rootAssembly.instances['RAILINST']
			railStress=stresses.getSubset(region=rail,position=INTEGRATION_POINT)
			s33=[]
			MaxPrinc=[]
			DirVec=[]
			f=open('RailbendingStress-'+str(j)+'.txt','w')
			dir=open('Princ-Direction-'+str(j)+'.txt','w')
			#
			#Displacements:
			disp=odb.steps['Step-2'].frames[-1].fieldOutputs['U']		
			#Designate the location of displacements:
			#-------------------------
			u1=disp.locations[0]
			u=disp.getSubset(location=u1)
			#Maximum displacements array initialization and file creation:
			#----------------------------
			Dis=[]
			princip0=[]
			Princip1=[]
			Princip2=[]
			for i in range(0,len(railStress.values)):
				str0=railStress.values[i].data[0] 	#SGMxx
				str1=railStress.values[i].data[1]	#SGMyy
				str2=railStress.values[i].data[2]	#SGMzz
				str3=railStress.values[i].data[3]	#SGMxy
				str4=railStress.values[i].data[4]	#SGMxz
				str5=railStress.values[i].data[5]	#SGMyz	
				ST=[[str0,str3,str4],[str3,str1,str5],[str4,str5,str2]] #Stress tensor at integration point in an element
				#Compute Eigenvalues of the stress tensor (principal stresses) and eigenvectors of the stress tensors (directions)
				evals,evecs=np.linalg.eig(ST)
				maxpr=max(evals)
				Princip0.append([evals[0])
				Princip1.append([evals[1])
				Princip2.append([evals[2])
				indp=list(evals).index(maxpr)
				astress=abs(railStress.values[i].data[2])
				DirVec.append([evecs[0][indp],evecs[1][indp],evecs[2][indp]])
				# PrincStr[i]=dict({"Max. princ.":maxpr,"Direction":DirVec})
				s33.append(astress)
				MaxPrinc.append(maxpr)
				f.write(str(astress)+'\n')
				dir.write(str(MaxPrinc[i])+'\t')
				[dir.write(str(DirVec[i][k])+'\t') for k in range(3)]
				dir.write('\n')
				#Extract vertical displacements
			for i in range(0,len(u.values)):
				d=u.values[i].data[1]
				Dis.append(d)
			#computing effective stress
			effstr0max=max(Princip0)
			effstr0min=0.0
			sigm0a=(effstr0max-effstr0min)/2.0
			sigm0m=(effstr0max+effstr0min)/2.0
			effctr1max=max(Princip1)
			if effctr1max<0.0:
				effctr1max=0.0
				effctr1min=min(Princip1)
			else:
			effctr1min=0.0
				continue
			sigm1a=(effstr1max-effstr1min)/2.0
			sigm1a=(effstr1max+effstr1min)/2.0
			effctr2max=max(Princip2)
			if effctr2max<0.0:
				effctr2max=0.0
				effctr2min=min(Princip2)
			else:
				effctr2min=0.0
				continue
			sigm2a=(effstr2max-effstr2min)/2.0
			sigm2a=(effstr2max+effstr2min)/2.0
			sigmae=1/sqrt(2)*((sigma0a-sigm1a)**2+(sigma1a-sigm2a)**2+(sigma2a-sigm0a)**2)**0.5
			sigmam=1/sqrt(2)*((sigma0m-sigm1m)**2+(sigma1m-sigm2m)**2+(sigma2m-sigm0m)**2)**0.5
			#Maximum displacement for a given load application
			
			maxDis=max(np.absolute(np.array(Dis)))
			MaxDisp.append(maxDis)
			#Find maximum bending stress
			Ms33=max(s33)
			MPr=max(MaxPrinc) #highest maximum principle in the rail
			Ind_MPr=MaxPrinc.index(MPr)
			MDir=DirVec[Ind_MPr]	
			#Append maximum bending stress on the rail in a different files
			maxs33.append(Ms33)
			MPrs.append(MPr)
			MDirs.append(MDir)
			g.write(str(Ms33)+'\n')
			gp.write(str(MPr)+'\t')
			[gp.write(str(MDir[m])+'\t') for m in range(3)]
			gp.write('\n')
			f.close()
			dir.close
			n.write(str(maxDis)+'\n')
		g.close()
		n.close()
		gp.close()
		#Find average maximum stress when the load has swept different position of the rail
		#Find also the standard deviation of the maximum bending stress
		avms33=(sum(maxs33))/nl
		sum_s33=0
		l=0
		while l<nl:
			sum_s33=sum_s33+(maxs33[l]-avms33)**2
			l=l+1
		SD_s33=sqrt(sum_s33/nl)
		#List these standard deviation for different standard deviations of the ties stiffness
		SD_s33_list.append([SD_s33])
		h.write(str(SD_s33)+'\n')
		t.write(str(avms33)+'\n')
		#Find the materials 
		#	
		odbm=openOdb('JobAuto-'+str(l0)+'.odb')
		Concr_E=[odbm.materials['CONCRETE'+str(m)].elastic.table[0][0] for m in range(nt)]
		Concr_P=[odbm.materials['CONCRETE'+str(m)].elastic.table[0][1] for m in range(nt)]
		[q1.write(str(Concr_E[n])+'\n') for n in range(nt)]
		[r1.write(str(Concr_P[n])+'\n') for n in range(nt)]
		q1.write('\n')
		r1.write('\n')
		
		#Calculate the standard deviations of Young' s modulus and stiffness's
		sumVE=0
		sumVP=0
		AveragE=(sum(Concr_E)/nt)
		AveragPr=(sum(Concr_P)/nt)
		#Compute new standards deviation for Young's modulus and Poisson ratio
		n=0
		while n<nt:
			sumVE=sumVE+(Concr_E[n]-AveragE)**2
			sumVP=sumVP+(Concr_P[n]-AveragPr)**2
			n=n+1
		stdE=sqrt(sumVE/nt)
		stdPr=sqrt(sumVP/nt)
		s.write(str(stdE)+'\n')
		pi.write(str(stdPr)+'\n')
		o.write(str(AveragE)+'\n')
	h.close()
	s.close()
	pi.close()
	o.close()
	t.close()
	q1.close()
	r1.close()