#!/usr/bin/python3
# -*- coding: utf-8 -*-

##################################################
## Modules
##################################################
## Python modules
import sys, os, argparse
from time import time, localtime, strftime, sleep
import pprint
##################################################
## Variables Globales
nomsEspeces = []
sequences = []
balises = []
version="3.4"

##################################################
## Functions
def checkParameters (arg_list):
	## Check input related options

	if (not arg_list.alignfile):
		print ('Error: No input file defined via option -i/--input !' + "\n")
		parser.print_help()
		exit(1)

#def checkVersion ():
	#lacmd=os.popen("grep -E ^version= /home/popphyl/PROGRAMME/GapCleaner.py | cut -f2 -d'\"'")
	#vp =lacmd.read().replace("\n","")
	#if float(version) < float(vp):
		#print('\nCurrent version '+version+' is not more recent ('+vp+')\nDo you  want to start programme? (y/n)')
		#answer = input()
		#while answer != "y" and answer != "n":
			#print("Please answer y or n !!!!!!!!!\n")
			#answer = input()
		#if answer == "n":
			#print("Exit programme\n")
			#exit(1)
		#if answer == "y":
			#print("Start programme with version "+version+"\n")


##################################################
## Main code
##################################################

if __name__=="__main__":
	
	# Initializations
	#checkVersion()
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog='Gapcleaner.py', description='''This Programme remove Gap with seuil from alignement''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help='display Gapcleaner version number and exit')
	parser.add_argument('-d', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")
	
	files = parser.add_argument_group('Input info for running')
	files.add_argument('-s', '--seuil',type=float, metavar="<float>", dest = 'seuil', default=0.25, help = '''Seuil ?\n\
																											0 = no gap in cleaned alignment,\n\
																											1 = all gap in cleaned alignment.\n\
																											(default = 0.25)''')

	files.add_argument('-S', '--removeStop', choices=('False','True'), dest = 'removeStop', default = 'False', help = 'Remove Stop codons (default = False)')
	files.add_argument('-t', '--type', type=str, metavar='<str>', dest='type', default='dna', help = 'Type of sequences : dna or aa ? (default = dna')

	files.add_argument('-i', '--input', metavar="<filename>", dest = 'alignfile', help = 'Name of file containing alignement')
	files.add_argument('-o', '--output', metavar="<filename>", dest = 'alignfileout', help = 'Name of output file containing alignement remove (default= nomAlign_clean.fasta)')


	args = parser.parse_args()
	# Check parameters
	checkParameters(args)
	
	# Welcome message
	print("#################################################################")
	print("#             Welcome in Gapcleaner (Version " + version + ")               #")
	print("#################################################################")
	print("\nStart time: ", start_time, "\n")

	#Lecture du fichier
	nomAlign = args.alignfile
	print(nomAlign, "en cours de lecture...\n")
	fichierAlign = open(nomAlign, "r")
	align = fichierAlign.readlines() 	# align est une liste avec dans chaque case une ligne
	fichierAlign.close()
	# Récupération seuil
	seuil = args.seuil

	for i in range(0,len(align)):
		if align[i][0] == '>' :
			nomsEspeces.append(align[i])
			balises.append(i)

	balises.append(len(align))

	for i in range(0,len(balises) - 1):
		sequences.append(align [ balises[i] +1  :  balises[i+1] ])

	# Enlever les \n des sequences
	for i in range(0, len(sequences)):
		for j in range(0, len(sequences[i])):
			sequences[i][j]=sequences[i][j].rstrip('\n')

	#Rassembler les sites de chaque séquences dans un seul str
	for i in range(0, len(sequences)):
		sequences[i] = ''.join(sequences[i])

	nbSitesInitial = len(sequences[0])

	if args.type == 'dna' and args.removeStop == True:
		# Vérification présence codon stop ['TAA', 'TAG'] et remplacement par NNN
		j=0
		for i in range(0, len(sequences)):
			while j <= len(sequences[0]):
				#if sequences[i][j:j+3] in ['TAA', 'TAG']:
				if sequences[i][j:j+3] in ['TAA', 'TGA','TAG']:
					print("Présence d'un codon stop dans la sequence",nomsEspeces[i].rstrip(),"au site",j+1,". Il sera remplacer par NNN")
					sequences[i] = sequences[i][:j]+"NNN"+sequences[i][j+3:]
				j += 3
			j=0
	if args.type == 'dna':
		# Vérification présence des lettre ACTG N
		j=0
		for i in range(0, len(sequences)):
			while j < len(sequences[0]):
				sequences[i][j].replace("a","A").replace("c","C").replace("t","T").replace("g","G")
				if sequences[i][j] not in ['A', 'C', 'T', 'G','-','N','n','a', 't', 'c', 'g']:
					print("Présence de la lettre",sequences[i][j], "dans la sequence",nomsEspeces[i].rstrip(),"au site",j,". Il sera remplacer par N")
					sequences[i] = sequences[i][:j]+"N"+sequences[i][j+1:]
				j += 1
			j=0
		
		
		
		
	if args.debug == "True":
		pprint.pprint(sequences)

	gapScore = [0]
	gapScore *= nbSitesInitial

	#Compter les gap dans chaque sites
	for i in range(0,nbSitesInitial):
		for seq in sequences:
			if seq[i] == "-" or seq[i] == "N" or seq[i] == "?" :
				gapScore[i] += 1
	
	#Separer les sites de chaque sequences dans des cases d'une liste
	for i in range(0, len(nomsEspeces)):
		sequences[i] = list(sequences[i])
	nogap=[""]
	nogap*= len(nomsEspeces)
	
	
	#Changer les gapScore en proportions
	if args.type=='dna':
		for i in range(0,nbSitesInitial,3):
			gapScore[i] = gapScore[i]/len(nomsEspeces)
			gapScore[i+1] = gapScore[i+1]/len(nomsEspeces)
			gapScore[i+2] = gapScore[i+2]/len(nomsEspeces)
			
			if gapScore[i]  <= seuil and gapScore[i+1]  <= seuil and gapScore[i+2]  <= seuil:
				for j in range (0, len(nomsEspeces)):
					#print(sequences[j])
					nogap[j] = nogap[j]+sequences[j][i]
					nogap[j] = nogap[j]+sequences[j][i+1]
					nogap[j] = nogap[j]+sequences[j][i+2]
		nbSitesFinal = len(nogap[0])

	elif args.type=='aa':
		for i in range(0,nbSitesInitial):
			gapScore[i] = gapScore[i]/len(nomsEspeces)
			
			if gapScore[i]  <= seuil:
				for j in range (0, len(nomsEspeces)):
					#print(sequences[j])
					nogap[j] = nogap[j]+sequences[j][i]
		nbSitesFinal = len(nogap[0])
	else:
		print('ERROR : Type is not dna or aa...')
		exit(0)


	#Réecrire le fichier
	if args.alignfileout == None:
		nomFichierSortie = nomAlign.split(".")[0]+"_clean.fasta" 
	else:
		nomFichierSortie = args.alignfileout
	
	fichierSortie = open(nomFichierSortie, "w")

	for i in range(0, len(nomsEspeces)):
		if nogap[i] != "-"*len(nogap[i]) and nogap[i] != "N"*len(nogap[i]) and nogap[i] != "n"*len(nogap[i]):
			fichierSortie.write(nomsEspeces[i])
			fichierSortie.write(nogap[i])
			fichierSortie.write("\n")
	fichierSortie.close()


	# Display a summary of the execution
	print('\nExecution summary:' + "")
	print('\t-',nomAlign, "rassemble",len(nomsEspeces),"espèces.")
	if len(nomsEspeces) == 4:
		print('\t- Seuil utilisé : 0.25')
	else:
		print('\t- Seuil utilisé :',seuil)
	print("\t- Il y a", nbSitesInitial, "sites dans l'alignement de départ.")
	print("\t- Il y a", nbSitesFinal, "sites dans l'alignement nettoyé --->", (1-(nbSitesFinal / nbSitesInitial)) *100 , "% des sites ont étés éliminés")
	print("\t- Fichier \"",nomFichierSortie,"\" créer")
	
		
	print("\nStop time: ", strftime("%d-%m-%Y_%H:%M:%S", localtime()))
	print("#################################################################")
	print("#                        End of execution                       #")
	print("#################################################################")
