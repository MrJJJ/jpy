#!/usr/bin/env python
# encoding: utf-8

import sys
import glob
import os

if len(sys.argv)==1:
	print '***************************'
	print '*** Welcome in multish ****'
	print '***************************'
	print "Parameters MUST be between SIMPLE quote"
	print "-c or --command : The command you want to execute with $f as target file (example : './mapnh param=$f'"
	print "-f --files : Your target files (default = '*.fasta')"
	print "-n --number : The rough number of slices to divide your target files (default '100')"
	print "-N --name : The basename of sub-sh files (default = 'z')"
	sys.exit()

command = ''
files = '*.fasta'
number = '100'
name = 'z'
submit = 0

for i in range(len(sys.argv)):
	if sys.argv[i]=='-c' or sys.argv[i]=='--command':
		command = sys.argv[i+1]
		print 'Your command :',command
	if sys.argv[i]=='-f' or sys.argv[i]=='--files':
		files = sys.argv[i+1]
		print 'Your files :',files
	if sys.argv[i]=='-n' or sys.argv[i]=='--number':
		number = sys.argv[i+1]
		print 'Number of sub-sh files to approach :',number
	if sys.argv[i]=='-N' or sys.argv[i]=='--name':
		name = sys.argv[i+1]
		print 'Name of su-sh files :',name
	if sys.argv[i]=='-s' or sys.argv[i]=='--submit':
		submit = 1
		print 'The sh-files will be sumbitted to SGE at the end of the script'

print '-----------------------------------------------------------------'

fileList = glob.glob(files)
numberOfFiles = len(fileList)

numberOfFilesPerSh = int(numberOfFiles/int(number))+1

listOfFileLists = []
fileChunks = range(0,numberOfFiles,int(numberOfFilesPerSh))

for i in range(len(fileChunks)-1):
	listOfFileLists.append(fileList[fileChunks[i]:fileChunks[i+1]])
listOfFileLists.append(fileList[fileChunks[-1]:])

#Verifying that all files are in listOfFileLists
for f in fileList:
	test=0
	for l in listOfFileLists:
		if f in l:
			test = 1
	if test == 1:
		pass
	elif test == 0:
		print '!!!!! Big problem,',f,'is not in the sh-files created !!!!!'

for i in range(len(listOfFileLists)):
	ofile = open(name+str(i)+'.sh', 'w')
	ofile.write('for f in ')
	ofile.write(' '.join(listOfFileLists[i]))
	ofile.write("\n")
	ofile.write('do ')
	ofile.write(command)
	ofile.write("\n")
	ofile.write('done')
	
print str(len(listOfFileLists))+' '+name+'.sh files created ('+str(len(listOfFileLists)-1)+' with '+str(len(listOfFileLists[0]))+' files, the last with '+str(len(listOfFileLists[-1]))+' files).'

print '---> You can submit these '+str(len(listOfFileLists))+' '+name+'.sh files in sge with a bash script'
print '---> (or relaunch the script with -s or --submit)'
print('')
print 'for f in '+name+'*.sh ; do qsub -cwd -q long.q -b y ./$f ; done'

if submit == 1:
	os.system('for f in '+name+'*.sh ; do qsub -cwd -q long.q -b y ./$f ; done')
