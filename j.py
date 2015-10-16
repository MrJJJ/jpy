# encoding utf-8

import math

########################
###   text parsing   ###
########################
def findall(string,sub): #Retourne la liste de tous les index d'un substring dans un string
	start=0
	indexlist=[]
	while True:
		index=string.find(sub,start)
		if index is -1: break
		indexlist.append(index)
		start = index+1
	return indexlist

def findall_re(s,substring): #renvoie les positions de tous les substrings du string s
	import re
	return [m.start() for m in re.finditer('substring', 's')]

#######################
###   list tricks   ###
#######################

def multisort(listOfLists): #sort a list of lists according to the first list
	sorted_listOfLists = []
	for i in listOfLists:
		a,b = zip(*sorted(zip(listOfLists[0],i)))
		sorted_listOfLists.append(b)
	return [list(x) for x in sorted_listOfLists]

def test_multisort(): 
	a=[2,1,3]
	b=['Deux','Un','Trois']
	c=['D','U','T']
	a,b,c = multisort([a,b,c])
	print(a,b,c)
	
class Csv(list): #Extension of class list to make cool stuff with csv
	def __init__(self, *args):
		list.__init__(self, *args)
		self.headers = self[0] #headers attribute
		for i,header in enumerate(self[0]):
			setattr(self,header,[x[i] for x in self[1:]]) #attribute columnName return the column as a list
			setattr(self,header+'_id',self[0].index(header))#attribute columnName_id return the id of the column
	def sort(self, column, reverse=False): #sort by a column
		return [y for (x,y) in sorted(zip(column,self[1:]),reverse=reverse)]


###############	
###   seq   ###
###############
def countACGT(seq):
    return sum([float(seq.count(x)) for x in 'ACGT'])

def pos(seq,posNumber):
	s=seq[posNumber-1::3]
	return(s)

def gc(s): #GC% d'un string
	return sum([float(s.count(x)) for x in 'GC']) / sum([float(s.count(x)) for x in 'ACGT'])

def g(s): #GC% d'un string
	return sum([float(s.count(x)) for x in 'G']) / sum([float(s.count(x)) for x in 'ACGT'])
def c(s): #GC% d'un string
	return sum([float(s.count(x)) for x in 'C']) / sum([float(s.count(x)) for x in 'ACGT'])
def a(s): #GC% d'un string
	return sum([float(s.count(x)) for x in 'A']) / sum([float(s.count(x)) for x in 'ACGT'])
def t(s): #GC% d'un string
	return sum([float(s.count(x)) for x in 'T']) / sum([float(s.count(x)) for x in 'ACGT'])

def gc1(s):
	s=s[0::3]
	return gc(s)
def gc2(s):
	s=s[1::3]
	return gc(s)
def gc3(s):
	s=s[2::3]
	return gc(s)


def hamming(seq1,seq2): #Compte le nombre de diff entre 2 string (distance de Hamming)
	return sum( [a != b and a != 'N' and b != 'N' and a != '-' and b != '-' for a,b in zip(seq1,seq2)] )




#################
###   align   ###
#################
def fasta2str(f): #Un alignement fasta en un string python
	onestring=''.join([x.rstrip() for x in open(f).readlines() if x[0] != '>'])
	return onestring

def fasta2tuple(f):
	fasta=open(f).readlines()
	seqName = [x.rstrip().replace('>','') for x in fasta if x[0] == '>']
	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	return (seqName,seq)

def fasta2dic(f):
	fastaStream = open(f)
	fasta = fastaStream.readlines()
	fastaStream.close()
	sp = [x.replace('\n','').rstrip().replace('>','') for x in fasta if x[0] == '>']
	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	dic = {}
	for i in range(len(sp)):
		dic[sp[i]] = seq[i]
	return dic 

def fasta2richdic(f):
	fastaStream = open(f)
	fasta = fastaStream.readlines()
	fastaStream.close()
	sp = [x.replace('\n','').rstrip() for x in fasta if x[0] == '>']
	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	GC3 = [ gc3(x) for x in seq ]
	length = [ len(x) for x in seq ]
	ngap = [ x.count('N')+x.count('-') for x in seq ]
	return {'sp':sp,'seq':seq, 'gc3':GC3, 'len':length, 'ngap':ngap}

def phy2dic(f):
    phyStream = open(f)
    phy = phyStream.readlines()
    phyStream.close()
    phyDic={}
    for l in phy:
        seqName = l.split()[0]
        seq = l.split()[1]
        if seqName in phyDic:
            phyDic[seqName] = phyDic[seqName]+seq
        else:
            phyDic[seqName] = seq
    return phyDic



def allFastaHaveSameLength(globstr): #Test if all seq of an alignment have same length
	import glob
	f=glob.glob(globstr)
	tab=[]
	for i in f:
		b=[x for x in open(i).readlines() if x[0]!='>']
		tab.append(all(len(x)==len(b[0]) for x in b))
	return tab
def allFastaHaveSameNumberOfSpecies(globstr): #Test if all alignments have same number of species
	import glob
	f=glob.glob(globstr)
	tab=[]
	for i in f:
		b=[x for x in open(i).readlines() if x[0]=='>']
		print(len(b))

def gapHomogenization(seqList): #Take a list of seq and put N when not ACGT in all positions
    seqListHomo = [['N']*len(seqList[0]) for i in range(len(seqList))]
    for i in range(len(seqList[0])):
        posIsToEliminate = 0
        for j in range(len(seqList)):
            if seqList[j][i] in 'ACGT':
                pass
            else:
                posIsToEliminate = 1
                break
        if posIsToEliminate == 1:
            pass
        else:
            for s in range(len(seqList)):
                seqListHomo[s][i] = seqList[s][i]
    return [''.join(x) for x in seqListHomo]

################
###   tree   ###
################


def treeToTaxa(tree): # treeTotaxa : renvoie une liste des noms des noeuds du gene mais ou les noeuds internes sont designes par l'ensemble des especes qu'ils contiennent (sp1_sp2_..._spN)
	nbNode = tree.count(":") #compte le nombre de ':' (=noeuds)
#	print("Il y a",nbNode,"noeuds")
	taxas=[] #Tableau qui va acceuillir tous les noms de taxons
	indexNode = [] # Tableau qui va acceuillir toutes les positions des ':' (=noeuds)
	for i in range(0,nbNode):
		indexNode.append(tree.replace(":","X",i).find(":")) #trouve la position du ieme ':' et la stocke dans le tableau indexNode
	for i in indexNode:
		j=1
		char=tree[i-j]
		while char.isdigit() or char==' ' or char=='.': #Remonte depuis le ':' jusqu'a la premiere lettre ou premiere parenthese qu'il rencontre
			j+=1
			char=tree[i-j]
		k=j # sauvegarde de la position j
		if char.isalpha():  # Si c'est une lettre (ou underscore pour 2 especes de tortues)
			while char.isalpha() or char=='_': #Tant que c'est une lettre on remontre pour chopper tout le nom de l'espece
				k+=1
				char=tree[i-k]
			taxa=tree[i-k+1:i-j+1] # =le nom de l'espece
			taxas.append(taxa) #Ajout du nom de l'sp dans le tableau taxas
		if char==')':       # Si c'est une parenthese
			parenthesisSum=-1
			while parenthesisSum!=0: # on cherche la parenthese ouvrante correspondante pour avoir toutes les especes qui sont contenues dans la branche ancestrale correspondante. 
				k+=1
				char=tree[i-k]
				if char==')':
					parenthesisSum-=1
				if char=='(':
					parenthesisSum+=1
			taxa=tree[i-k+1:i-j+1] # = portion de l'arbre auquel correspond le ':'
			taxa=taxa.replace('1','')
			taxa=taxa.replace('2','')
			taxa=taxa.replace('3','')
			taxa=taxa.replace('4','')
			taxa=taxa.replace('5','')
			taxa=taxa.replace('6','')
			taxa=taxa.replace('7','')
			taxa=taxa.replace('8','')
			taxa=taxa.replace('9','')
			taxa=taxa.replace('0','')
			taxa=taxa.replace('.','')
			taxa=taxa.replace('(','')
			taxa=taxa.replace(')','')
			taxa=taxa.replace(':','')
			taxa=taxa.replace(' ','')
			taxa=taxa.replace(',','+')
			taxas.append(taxa) # Ajout de sp1_sp2_sp3... dans le taleau taxas
	return taxas



def taxaConversion(taxas,spphylo): # Prend tableau taxas cree avec fonction treeToTaxa et fichier spphylo (=1 ligne par especes avec toute leur phylogenie, ex : Mus,Rodentia,...,Mammal) -> renvoie une liste des noms des noeuds du gene mais ou les noeuds internes sont designes par l'ensemble des especes qu'ils contiennent (sp1_sp2_..._spN)
	taxasConverted=[] # Resultat final de la fonction = Liste qui va contenir toutes les conversion
	for taxa in taxas: #Pour chaque noeuds non convertits
		taxa=taxa.split(',')#On saucissone tous les '_'
		if len(taxa) > 1: #Si on est dans un noeud interne
			phylodico = {} #Dico qui va associer chaque taxon au nombre de fois ou rencontrera dans spphylo
			for sp in taxa: #Pour chaque espece dans taxa -> On va creer un dico
				for pline in spphylo: #On parourt chaque ligne de spphylo
					pline=pline.split(',')#On saucissone selon les ','
					if pline[0] == sp: #On cherche ligne de l'espece orrespondante dans spphylo
						phyloorder=pline #Sauvegarde d'une ligne de spphylo avec bonne espece
						for p in pline: #Pour chaque case de la ligne de spphylo
							if p in phylodico: 
								phylodico[p]+=1 #Si dans dico, on ajoute 1 au nombre de fois ou on l'a rencontree
							else:
								phylodico[p]=1 #Si pas dans le dico on l'initialise a 1.
			position=9999999999999999 #On initialise position en derniere position
			for cle in phylodico: #Pour chaque cle du dico
				if phylodico[cle]==max(phylodico.values()): #On selectionne toutes les cles qui ont valeur maximale
					if phyloorder.index(cle) < position:
						position=phyloorder.index(cle)
			taxaConverted=phyloorder[position]
			taxasConverted.append(taxaConverted)
		if len(taxa) == 1: #Si on est dans un noeud terminal
			taxasConverted.append(taxa[0])
	return taxasConverted



def treeToGC(tree): # arbre -> liste des labels (bootstraps, GC de nhml)
	nbNode = tree.count(":")   # Compte le nombre de noeuds donc d'especes (ancetre commun compris)
#	print("Il y a",nbNode,"noeuds")
	GC=[]                    # Tableau qui va acceuillir toutes les valeurs de GC
	indexNode = []				# Tableau qui va acceuillir toutes les positions des ':' (=noeuds)
	for i in range(0,nbNode):
		indexNode.append(tree.replace(":","X",i).find(":"))  # Trouve la position du ieme ':' et la stocke dans le tableau indexNode
	for i in indexNode:
		j=1
		char=tree[i-j]
		while char.isdigit() or char=='.':
			j+=1
			char=tree[i-j]
		GCsp=tree[i-j+1:i]       # On garde la partie correspondante a la valeur de GC
		GC.append(GCsp)
		#print(char)
	return GC



def treeToBL(tree): # arbre -> liste des longueurs de branches
	nbNode = tree.count(":")   # Compte le nombre de noeuds donc d'especes (ancetre commun compris)
#	print("Il y a",nbNode,"noeuds")
	BL=[]                    # Tableau qui va acceuillir toutes les valeurs de longueur de branches (BL)
	indexNode = []				# Tableau qui va acceuillir toutes les positions des ':' (=noeuds)
	for i in range(0,nbNode):
		indexNode.append(tree.replace(":","X",i).find(":"))  # Trouve la position du ieme ':' et la stocke dans le tableau indexNode
	for i in indexNode:
		j=1
		char=tree[i+j]
		while char.isdigit() or char=='.':
			j+=1
			char=tree[i+j]
		BLsp=tree[i+1:i+j]       # On garde la partie correspondante a la valeur de longueur de branhe
		BL.append(BLsp)
		#print(char)
	return BL


#################
###   Stats   ###
#################

