import os
import shutil
from os import path
import subprocess
#import numpy as np
startpath1=r"C:\Users\celestink\Documents\TestingStation\changAv"
newpath = path.join(startpath1,"Rerun1")
batchfile="\ProAn.bat"
if not os.path.exists(newpath): os.makedirs(newpath)
g=open('FilesToRerun.txt','w')
#loop over all files to see if they exist
checkfiles=[0+i for i in range(19)]+[190+j for j in range(19)]+[380+j for j in range(19)]+[570+j for j in range(19)]+[760+j for j in range(19)]+[950+j for j in range(19)]
for k in checkfiles:
	if  os.path.exists(startpath1+"\JobAuto-"+str(k)+".odb"): 
		continue	
	else:
		#write names of missing files in a different file
		g.write('JobAuto-'+str(k)+'.inp'+'\n')
		#copy corresponding inp file in a different folder
		shutil.copy2(startpath1+"\JobAuto-"+str(k)+".inp",newpath+"\JobAuto-"+str(k)+".inp")
		
g.close()
batchpath=newpath+batchfile
#run again the inp files the didn't run before 
shutil.copy2(startpath1+batchfile,newpath+batchfile)
# p = subprocess.Popen(batchpath, shell=True)
# p.communicate()
# print (p.returncode) # is 0 if success

