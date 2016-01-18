
## Importing modules
import abaqus
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
#P=Fd/Area1
a=9.2
b=7.8
PLoad=1360.700002
P=PLoad*2.0/3.0
for SDi in range(10): #SDi define incremental standard deviation index of tie Young's moduli
	##INFORMATION ON TIES
	##Input dimensional parameters, metric units, mm
	#Tie spacings
	tsp=21*25.4 # distance between centrelines of two successive ties; units in mm
	wt=228.6 # Tie width on z axis
	st=tsp-wt # gap distance between two successive ties
	##Information on ties
	nt=20#number of ties
	ht=177.8 #height of the tie on y axis
	Lt=457.2 #length of the tie on x axis
	#Longitudinal starting point of ties
	l0=st/2.0
	f = open('tieModuli_SD_Incr_'+str(SDi)+'.txt', 'w')
	for Nrand in range(1):
		#Generating random ties Young's Modulus and poisson ratio for concrete ties
		stE0=60.0
		muE0=100.0
		dmuE=50.0
		muE=muE0+SDi*dmuE #mean of Young's modulus,
		stE=stE0  #Standard deviation  of Young's modulus,
		muPr=0.18
		stPr=0.05
		Pr=0.15
		## This function generate a sample of random numbers nt from of Gaussian distribution 
		# such that their average and standard deviation are the same as the pool fixed mean muE and standard deviation stE
		ConcrE1=[]
		ConcrE=[]
		# def randomstiffness(stE,muE,nt,ConcrE1):
		keepsampling=0.0
		while keepsampling<nt:
			ConcrEi=[random.gauss(muE,stE) for i in range(nt)]
			# for l in range(nt):
				# if ConcrE[l]<0:
					# ConcrE[l]=muE0+abs(ConcrE[l])
				# else:
					# continue
			ConcrPr=[random.gauss(muPr,stPr) for i in range(nt)]
			#Force the mean to be the same as requested
			sumstiff=sum(ConcrEi)
			newmean=sumstiff/nt
			Xdev=(nt*muE-sumstiff)/nt
			ConcrEc=[ConcrEi[i]+Xdev for i in range(nt)] 
			
			#Compute the standard deviation of the new sample
			ConcrEsq=[(ConcrEc[i])**2 for i in range(nt)]
			StD=sqrt(abs(sum(ConcrEsq)/nt-muE**2))
			
			#Deviation from requested variances
			dStD2=stE**2-StD**2
			
			#Force the standard deviation of the sample to be equal to the standard deviation of the whole population
			#dstD is subtracted from numbers less than the median and and added to numbers greater than the median , that is:
			#Compute the median
			
			ListE=ConcrEc
			ListE.sort() #sort numbers in a list from the lowest to greatest
			Xmed=1.0/2.0*(ListE[nt/2-1]+ListE[nt/2])
			##Divide the list into two lists separated by the mean
			ListE1=[ListE[i] for i in range(nt/2)]
			ListE2=[ListE[i] for i in range(nt/2,nt)]
			muE1=(sum(ListE1))/(nt/2.0)
			muE2=(sum(ListE2))/(nt/2.0)
			dmuE=muE2-muE1
			Determ=dmuE**2+4.0*dStD2
			
			if Determ<0.0:
				keepsampling=0.0
			else:
				dStD=0.5*dmuE-0.5*sqrt(dmuE**2+4.0*dStD2)
				#perturb the numbers to have requested standard dev and mean
				ConcrE=[]
				for i in range(nt):
					if ConcrEc[i]<Xmed:
						Ec=ConcrEc[i]+dStD
						ConcrE.append(Ec)
					elif ConcrEc[i]>Xmed:
						Ec=ConcrEc[i]-dStD
						ConcrE.append(Ec)
					else:	
						Ec=ConcrEc[i]
						ConcrE.append(Ec)
				random.shuffle(ConcrE)
				keepsample=[]
				for i in range(nt):
					if ConcrE[i]<=20.0:
						#print (tie modulus error)
						keepsample.append(0)
					else:
						keepsample.append(1)
				keepsampling=sum(keepsample)
			# saving:
		#dataModulus=ConcrE
		#ConcrE=randomstiffness(stE,muE,nt,ConcrE1)
		#print ConcrE
		[f.write(str(ConcrE[k])+'\n') for k in range(nt)]        # column names
		 # numpy.savetxt(f, numpy.array([ConcrE, ConcrPr]).T)
		#Compute new averages 
		sumVE=0.0
		sumVP=0.0
		AveragE=(sum(ConcrE)/nt)
		AveragPr=(sum(ConcrPr)/nt)
		#Compute new standards deviation for Young's modulus and Poisson ratio
		n=0
		while n<3:
			sumVE=sumVE+(ConcrE[n]-AveragE)**2
			sumVP=sumVP+(ConcrPr[n]-AveragPr)**2
			n=n+1
		stdE=sumVE/nt
		stdPr=sumVP/nt
		#INFORMATION ON THE RAIL
		hr=186.0 #height on y axis
		wr=73.0	#width on x axis
		Lr=nt*tsp	#length on z axis
		re=20.0
		te=20.0
		
		#base coordinate of the rail
		R_vx=((wr/2.0,hr/2.0),(-wr/2.0,hr/2.0),(-wr/2.0,-hr/2.0),(wr/2.0,-hr/2.0)) #Define rail base vertices
		#base coordinate of the tie
		T_vx=((Lt/2.0,-hr/2.0),(-Lt/2.0,-hr/2.0),(-Lt/2.0,(-hr/2.0-ht)),(Lt/2.0,(-hr/2.0-ht))) #Define Tie base vertices
		#create surfaces for the rail loads
		#Datum plane for partition
		nl=nt-1 #number of load surfaces
		ll=pi*a*b/wr #longitudinal length of load surface in mm
		#to=(wt-ll)/2.0 #Offset distance from the load surface above the tie
		#dl=st/2.0+to-ll/2.0;#offset distance between the beginning of two successive load surface
		#(Wu & Thompson, 1999, the wheel load of 75kN) #traction on the rail  #My ORIGINAL TRACTION = 19.2
		##BEGIN SKETCHING AND CREATING THE MODEL
		# Area1=wr*ll
		# F=142.3*10**3 #Newtons
		# theta=0.25
		# Fd=(1+theta)*F #dynamic load
		# P=Fd/Area1
		##Base coordinates
		for j in range(nl):
			##Sketching the rail base
			l=19*(10*SDi+Nrand)+j
			myModel=mdb.Model(name='SimTrack'+str(l))
			sketch1=myModel.ConstrainedSketch(name='Rail section', sheetSize=400.0)
			sketch1.Line(point1=R_vx[0],point2=R_vx[1])	
			sketch1.Line(point1=R_vx[1],point2=R_vx[2])
			sketch1.Line(point1=R_vx[2],point2=R_vx[3])	
			sketch1.Line(point1=R_vx[3],point2=R_vx[0])
			##Create Rail part by extrusion
			myrailPart=myModel.Part(name='Rail',dimensionality=THREE_D,
				type=DEFORMABLE_BODY)
			myrailPart.BaseSolidExtrude(depth=Lr, sketch=sketch1)
			#
			#CREATE THE STEEL MATERIAL AND THE RAIL SECTION
			myMaterial=myModel.Material(name='Steel')
			myMaterial.Elastic(table=((200000.0, 0.3),))
			myModel.HomogeneousSolidSection(material='Steel', name='RailSteel')
			myrailPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
					cells=myrailPart.cells.findAt(((wr/4.0,(hr/4.0),Lr/2.0), 
						), )), sectionName='RailSteel', thicknessAssignment=
					FROM_SECTION)	
			#Creating a rail Instance
			myModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
			RailInst=myModel.rootAssembly.Instance(dependent=ON, name='RailInst', 
						part=myrailPart)
			#Create surfaces for which node sets will be created
			myModel.rootAssembly.Set(faces= RailInst.faces.findAt(((
				0, 0, 0), ), ), name='Set-1'+str(j))
			myModel.rootAssembly.Set(faces= RailInst.faces.findAt(((
				0, 0, Lr), ), ), name='Set-2'+str(j))
			#Create step 2 for the load
			myModel.StaticStep(initialInc=0.1, name='Step-2', previous='Initial')
			##APPLYING LOADS ON THE RAIL
			#Start by creating datum planes for partitioning
			dlp11=tsp*(j+0.5)-ll/2.0 #starting distance of load surface above the tie
			dlp12=dlp11+ll #ending distance of load surface above the ties
			dlp21=tsp*(j+1.0)-ll/2.0 #starting distance of load surface above the tie
			dlp22=dlp21+ll #ending distance of load surface above the ties
			dlp=((dlp11,dlp12),(dlp21,dlp22));
			v=1
			planep1=myrailPart.DatumPlaneByPrincipalPlane(offset=dlp[v][0], 
				principalPlane=XYPLANE)
			planep2=myrailPart.DatumPlaneByPrincipalPlane(offset=dlp[v][1], 
				principalPlane=XYPLANE)
			#Do partitions of the rail
			partitionp1=myrailPart.PartitionFaceByDatumPlane(datumPlane=
				myrailPart.datums[3], faces=myrailPart.faces.findAt(((wr/4.0, hr/2.0, dlp[v][0]+ll/2.0), ), ))
			partitionp2=myrailPart.PartitionFaceByDatumPlane(datumPlane=myrailPart.datums[4], faces=myrailPart.faces.findAt(((wr/4, hr/2, 
				dlp[v][1]+ll/2), ), ))
			#Apply loads at partitioned surface
			myModel.SurfaceTraction(createStepName='Step-2', 
				directionVector=((0.0, 0.0, 0.0), (0.0, -1.0, 0.0)), distributionType=
				UNIFORM, field='', localCsys=None, magnitude=P, name='Load-1'+str(j), region=
				Region( side1Faces=RailInst.faces.findAt(((wr/4.0, hr/2.0, (dlp[v][0]+ll/2.0)), ), )), resultant=ON, traction=GENERAL)
			##Mesh the rail part
			myrailPart.setElementType(elemTypes=(ElemType(
				elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15, 
					elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD)), 
				regions=(myrailPart.cells.findAt(((-wr/4.0, 
					hr/4.0, Lr/4.0), ), ), ))
			myrailPart.seedPart(deviationFactor=0.4, minSizeFactor=0.4, size=re)		
			myrailPart.generateMesh()
			##create a dummy element
			sketchd=myModel.ConstrainedSketch(name='dummy section', sheetSize=400.0)
			sketchd.Line(point1=(0,0),point2=(re,0))
			sketchd.Line(point1=(re,0),point2=(re,re))
			sketchd.Line(point1=(re,re),point2=(0,re))	
			sketchd.Line(point1=(0,re),point2=(0,0))
			mydummyPart=myModel.Part(name='dummy',dimensionality=THREE_D,
				type=DEFORMABLE_BODY)
			mydummyPart.BaseSolidExtrude(depth=re, sketch=sketchd)
			mydummyPart.SectionAssignment(offset=0.0, 
						offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
						cells=mydummyPart.cells.findAt(((0,0,Lr), 
						), )), sectionName='RailSteel', thicknessAssignment=
							FROM_SECTION)
			#Assigning dammy element material
			mydummyPart.SectionAssignment(offset=0.0, 
				offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
				cells=mydummyPart.cells.findAt(((re/2.0,(re/2.0),re/2.0), 
				), )), sectionName='RailSteel', thicknessAssignment= FROM_SECTION)				
			#Creating a dummy Instance
			myModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
			DummyInst=myModel.rootAssembly.Instance(dependent=ON, name='dmyInst', part=mydummyPart)
			#Translate the the dummy instance at a distance Lr
			myModel.rootAssembly.translate(instanceList=('dmyInst', ), 
					vector=(0.0, 0.0, Lr))	
			#Applying Encastred BC
			myModel.EncastreBC(createStepName='Initial', name='BC-dummy'+str(i), region=Region(faces=DummyInst.faces.findAt(((re/2, re/2, Lr+re),),)))
			#Mesh the dummy part
			mydummyPart.setElementType(elemTypes=(ElemType(elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15, 
				elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD)), 
				regions=(mydummyPart.cells.findAt(((re/2.0, re/2.0, re/e), ), ), ))
			mydummyPart.seedPart(deviationFactor=0.5, minSizeFactor=0.5, size=re)
			mydummyPart.generateMesh()
			#create dummy node
			myModel.rootAssembly.Set(name='Set-3'+str(j), vertices=DummyInst.vertices.findAt(((0.0, 0.0, Lr), ), ))
			##CREATING TIE PARTS
			for i in range(nt): 
			## Start sketching the ties
				zi=l0+i*(wt+st)
				Tiesketch2=myModel.ConstrainedSketch(name='Tiesection'+str(i), sheetSize=800.0)
				Tiesketch2.Line(point1=T_vx[0],point2=T_vx[1])	
				Tiesketch2.Line(point1=T_vx[1],point2=T_vx[2])
				Tiesketch2.Line(point1=T_vx[2],point2=T_vx[3])	
				Tiesketch2.Line(point1=T_vx[3],point2=T_vx[0])
			#Create Tie part by extrusion
				mytiePart=myModel.Part(name='Tie'+str(i),dimensionality=THREE_D,
				type=DEFORMABLE_BODY)
				mytiePart.BaseSolidExtrude(depth=wt, sketch=Tiesketch2)
			#CREATE THE TIE CONCRETE  MATERIALS 
				myModel.Material(name='Concrete'+str(i))
				myModel.materials['Concrete'+str(i)].Elastic(table=((ConcrE[i], Pr),))
				#Create section
				myModel.HomogeneousSolidSection(material='Concrete'+str(i), name='ConcreteTie'+str(i))
				#Assign section
				mytiePart.SectionAssignment(offset=0.0, 
					offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
					cells=mytiePart.cells.findAt(((-Lt/4.0, (-hr/2.0-ht/2.0), wt/2.0), 
					), )), sectionName='ConcreteTie'+str(i), thicknessAssignment=FROM_SECTION)
			
				#Instancing each tie
				myModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
				mytieInst=myModel.rootAssembly.Instance(dependent=ON, name='Ties'+str(i), 
					part=mytiePart)
				#Translate the tie at a distance Tr along the rail
				Tr=st/2+i*tsp
				myModel.rootAssembly.translate(instanceList=('Ties'+str(i), ), 
					vector=(0.0, 0.0, Tr))
				##Applying Encastred BC
				myModel.EncastreBC(createStepName='Initial', name='BC-Tie'+str(i), region=Region(
					faces=mytieInst.faces.findAt(((Lt/4.0, (-hr/2.0-ht), (Tr+wt/2.0)), ), )))
				
				##Apply tie constraints to each tie and the rail
				myModel.Tie(adjust=ON, constraintEnforcement=NODE_TO_SURFACE, 
					master=Region(side1Faces=RailInst.faces.findAt(((wr/4.0, -hr/2.0, st/4.0), ), )), name='Constraint-1'+str(i), 
					positionToleranceMethod=COMPUTED, slave=Region(side1Faces=mytieInst.faces.findAt(
					((-3*Lt/8.0, -hr/2.0, (Tr+wt/2.0)), ), )), thickness=ON, tieRotations=ON)
					
				##Mesh each tie part
				mytiePart.setElementType(elemTypes=(ElemType(elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15, 
					elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD)), 
					regions=(mytiePart.cells.findAt(((Lt/4.0, (-hr-ht)/2.0, wt/2.0), ), ), ))
				mytiePart.seedPart(deviationFactor=0.2, minSizeFactor=0.2, size=te)
				mytiePart.generateMesh()
			#Set Field output requests
			myModel.FieldOutputRequest(createStepName='Step-2', name=
				'F-Output-'+str(l), variables=('S', 'MISES', 'MISESMAX',
				'E', 'PE', 'PEEQ', 'PEEQT', 'EE','ER','U', 'UT', 'TF', 'NFORC', 'BF', 'P'))
			#Apply Periodic boundary conditions

			set1=myModel.rootAssembly.sets['Set-1'+str(j)].nodes
			set2=myModel.rootAssembly.sets['Set-2'+str(j)].nodes
			set3=myModel.rootAssembly.sets['Set-3'+str(j)].nodes
			# node1=mdb.models['Model-1'].rootAssembly.instances['fullrail-1'].nodes
			node2=DummyInst.nodes
			set2new=[]
			set1new=[]
			x2s=[]
			y2s=[]
			z2s=[]
			z1s=[]
			# cst=152.4
			#compute the minimum distance in xy plane between one node from 
			#one face and each node from the other face and create a list of matched nodes
			x=[]
			y=[]
			z=[]
			def ArrNodeSet(set2,set2new):
				for n in range(len(set2)):
					x1=set2[n].coordinates[0]
					y1=set2[n].coordinates[1]
					z1=set2[1].coordinates[2]
					labe=set2[n].label
					x2s.append(x)
					y2s.append(y)
					z2s.append(z)
					set2new.append([labe,x1, y1, z1])
				z2cst=max(z2s)
			ArrNodeSet(set2,set2new)
			ArrNodeSet(set1,set1new)
			Tm=1 #Condition to match coordinates in set 2 to coordinates in set 1 else vice versa.
			NodeEqxy=[]
			set2n=[]
			def editingNodes(set1new,set2new,Tm):
				# set2n=[]
				for i in range(len(set1new)):
					d12xyp=[((set2new[n][1]-set1new[i][1])**2+(set2new[n][2]-set1new[i][2])**2)**0.5 for n in range(len(set2new))]
					dnmin=min(d12xyp) 
					mlab=d12xyp.index(dnmin) # Index of matching node in set2 to a node of index i in set1
					if (Tm==1):
						nameOfset='Set-2'+str(j)
						#lab1 and lab2 are labels of matching nodes.
						lab1=set1new[i][0]
						lab2=set2new[mlab][0]
						no2nx=set1new[i][1]
						no2ny=set1new[i][2]
						no2nz=set1new[1][3]
						set2n.append([no2nx,no2ny,no2nz,i,mlab])
						NodeEqxy.append([set1new[i],lab1,lab2])
						myModel.rootAssembly.editNode(coordinate1=set2n[i][0], 
						coordinate2=set2n[i][1], coordinate3=set2n[i][2],
						nodes=myModel.rootAssembly.sets[nameOfset].nodes[set2n[i][4]])	
					elif (Tm==2):
						nameOfset='Set-1'+str(j)
						labn=set1new[nj][0]
						no2nx=set2new[i][1]
						no2ny=set2new[i][2]
						no2nz=set2new[1][3]
						set2n.append([no2nx,no2ny,no2nz,mlab,lab2])
						NodeEqxy.append([set1new[i],lab1,lab2])
			#start editing nodes with new coordinates keeping the node number the same
						myModel.rootAssembly.editNode(coordinate1=set2n[i][0], 
						coordinate2=set2n[i][1], coordinate3=set2n[i][2],
						nodes=myModel.rootAssembly.sets[nameOfset].nodes[set2n[i][3]])		
			##find extra nodes on the other face that has more nodes
			#set2new-set2n
			editingNodes(set1new,set2new,Tm)
			if (range(len(set2))>range(len(set1))):
				set2node=[set2new[n][0] for n in range(len(set2new))]
				set2noden=[set2n[n][0] for n in range(len(set2n))]
				hangnode=list(set(set2node)-set(set2noden))
				#Edit these nodes 
				SetToEdit=[]
				def indexOfnodes(hangnode,set2,SetToEdit):
			#This function returns the set of nodes with their coordinates and their index in the nodes set)
					for k in range(len(hangnode)):
						NoToEdit=hangnode[k]
						for l in range(len(set2)):
							Nodenum=set2[l].label
							if (Nodenum==NoToEdit):
								xEdit=set2[l].coordinates[0]
								yEdit=set2[l].coordinates[1]
								zEdit=zcst
								SetToEdit.append([Nodenum,xEdit,yEdit,ZEdit,l])
								break
							elif (Nodenum!=NoToEdit):
								continue
				Tmm=2
				editingNodes(SetToEdit,Set1new,Tmm,'Set-2'+str(j))
				#Print nodes numbers pairs of the same coordinates
				print 'Equal x y coordinates node pairs are', NodeEqxy
				print 'hanging nodes are', SetToEdit
			#	Start defining boundary conditions
			elif (range(len(set2))==range(len(set1))):
				for i in range(len(set1)):
			#create the names of node sets
					setname1='set1a-'+str(i)
					setname2='set2a-'+str(i)
					setname3='set3'+str(j)
			##creates the node sets
					# node1=mdb.models['Model-1'].rootAssembly.sets['Set-1'].nodes
					# node2=mdb.models['Model-1'].rootAssembly.sets['Set-2'].nodes
					# node3=mdb.models['Model-1'].rootAssembly.sets['Set-3'].nodes
					st1=set1[set2n[i][3]:set2n[i][3]+1]
					st2=set2[set2n[i][4]:set2n[i][4]+1]
					# nodelabel1=nst1.label
					# nodelabel2=nst2.label
					# st1=set1[nodelabel1-1:nodelabel1]
					# st2=set1[nodelabel2-1:nodelabel2]
					#use global index sets to create node sets
			#Create a dummy node
					st3=set3;
			#Create the node sets
					myModel.rootAssembly.Set(name=setname1, nodes=st1)
					myModel.rootAssembly.Set(name=setname2, nodes=st2)	
					myModel.rootAssembly.Set(name=setname3, nodes=st3)	
			#Apply constraints equations for periodic boundary conditions	
					myModel.Equation(name='ConstraintX1-'+str(i), terms=((1.0, 
						setname1, 1), (-1.0, setname2, 1)))
					myModel.Equation(name='ConstraintY1-'+str(i), terms=((1.0, 
						setname1, 2), (-1.0, setname2, 2)))
					myModel.Equation(name='ConstraintZ1-'+str(i), terms=((1.0, 
						setname1, 3), (-1.0, setname2, 3),(1, setname3, 3)))
			## Create jobs
			# mytrackJob=mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
			# explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
			# memory=90, memoryUnits=PERCENTAGE, model='SimTrack'+str(j), modelPrint=OFF, 
			# multiprocessingMode=DEFAULT, name='JobAuto'+str(j), nodalOutputPrecision=SINGLE, 
			# numCpus=1, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
			# userSubroutine='', waitHours=0, waitMinutes=0)
			# mytrackJob.submit(consistencyChecking=OFF)
			
			mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
			explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
			memory=90, memoryUnits=PERCENTAGE, model='SimTrack'+str(l), modelPrint=OFF, 
			multiprocessingMode=DEFAULT, name='JobAuto-'+str(l), nodalOutputPrecision=SINGLE, 
			numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
			waitHours=0, waitMinutes=0)	
			mdb.jobs['JobAuto-'+str(l)].writeInput(consistencyChecking=OFF)
		f.close()