import os 
import re

def read_file(directory, pro_name, op=True):
	f = open(directory+pro_name, "r")
	lines=[]
	if op==True:
		next(f)##----for not read the first line!!! be careful with this!!
	#first_line = f.readline()
	for line in f:
		lines.append(re.split(',|;',line))
	f.close()
	#print(lines)
	return lines

def add_interaction(number):
	if number==15:
		return 'aromatic stacking'
	elif number==16:
		return 'hydrogen bond'
	elif number==17:
		return "\n" #'hydrogen bond water'
	elif number==18:
		return 'hydrophobic'
	elif number==19:
		return 'repulsive'
	elif number==20:
		return 'salt bridge'
	else:
		return 'not definied!'

def replace_quotation_marks(file_n, new_file):
	f=open(file_n, "r")

	#lines=f.readlines()
	line=f.readline()
	f.close()
	#print(len(line))
	#print(line)
	new_line=line.replace("\'", "\"")
	#new_line=line
	##"(HETATM|ATOM)(\_\d+\_\w+\_\w+\_\d+\_\w)\"/\'$1$2\'/", "\""
	#re.sub("\'\s[HETATM]|[ATOM](_d+_w+_w+_d+_\w)\s\'\1\2 \'", "\"", new_line)
	
	#print(new_line)
	#test= "h-->[\'1\':key]"
	#teste_new=test.replace("\'", "\"")
	#print(test)
	#print(teste_new)

	fo = open(new_file, "w")
	fo.write(new_line)
	fo.close()
	print("quotation marks replaced!")
	print(new_file+" created!")

def construct_file_output(dir_set,protein_name):
	# Reading file
	lines=read_file(dir_set, protein_name) ## ricin directory 
	#print (lines)
	#creating the ARFF file
	header="%\n%CDK2 ligand-protein interactions\n%Author: Alexandre Fassio\n%Pre-processed by: Susana Medina\n%Date: November 2014\n"
	header=header+"@RELATION cdk2\n"
	header=header+"@ATTRIBUTE Atom1 NOMINAL\n"
	header=header+"@ATTRIBUTE Type1 NOMINAL\n"
	header=header+"@ATTRIBUTE Atom2 NOMINAL\n"
	header=header+"@ATTRIBUTE Type2 NOMINAL\n"
	header=header+"@ATTRIBUTE Distance NUMERICAL\n"
	header=header+"@ATTRIBUTE Interaction NOMINAL\n"
	header=header+"@DATA\n"

	data=""
	att="protein_id,atom1,type1,atom2,type2,distance,interaction\n"


	#choose name of file??
	#f.write(header)
	#f.write(att)
	l_line=""
	inter_length=6

	for l in lines:
		l_line=protein_name.split('.')[0]+","
		for i in range(0,len(l)-1):
			l_line=l_line+l[i]+","
		data=data+l_line+l[len(l)-1]
		#if l[len(l)-1] != '\n': print(l_line+l[len(l)-1])
		l_line=""

	return data

def construct_output(dir_set, protein_name, delim=':'):
	#delimitator is : for contact new format files
	# Reading file
	lines=read_file(dir_set, protein_name,False)
	data=""
	l_line=""
	fl=""
	for i in range(0,len(lines[0])-7):
		fl=fl+lines[0][i]+","
	fl=fl+lines[0][len(lines[0])-7]+",interaction"+"\n"
	print("fl----->"+fl+"\n")
	contacts=lines[1:]
	for l in contacts:
		l_line=protein_name.split(delim)[0]+","

		for i in range(0,len(l)):
			if i >14:
				if int(l[i])==1:
					l_line=l_line+add_interaction(i)+"\n"
					#print(add_interaction(i))
			else:
				l_line=l_line+l[i]+","
		#data=data+l_line+l[len(l)-1]
		#print(l)
		data=data+l_line
		#if l[len(l)-1] != '\n': print(l_line+l[len(l)-1])
		l_line=""
	first_l=fl
	print(fl+data)
	return fl, data

def parse_files(file_pro_lig, file_out, dir_set):
	print ("Begining parsing...")
	#dir set --> directory of output file
	#Reading list of contact files
	#list_proteins= input() ##automatic process
	att="protein_id,Atom1,Type1,Atom2,Type2,Distance,Interaction\n"
	##redefine attributes read from file! first line!

	f = open(file_pro_lig, 'r')
	list_pro_lig = f.readlines()
	f.close()
	#print(list_pro_lig)
	f = open(str(file_out+".csv"), 'w')
	tmp=""
	fl=""
	co=""
	#print list_pro_lig
	for protein in list_pro_lig:
		#tmp=tmp+construct_file_output( dir_set,protein.rstrip('\n')) ## old format csv files
		if protein!="\n":
			fl, co=construct_output( dir_set,protein.rstrip('\n'))
			tmp=tmp+co
			print (str(protein.rstrip('\n')+ " processed \n"))
	fl="protein_id,"+fl
	print("fl value--->"+fl)
	f.write(fl+tmp)
	f.close()

	print ("File created!")
	print ("End!")

def create_model_pro_lig(dir_data, list_pro_lig ):
	
	print ("Creating model interactions protein-ligand...")

	#Reading list of contact files
	#file_pro_lig= "list_files_inter.txt" #"list_test.txt" # "list_files.txt"

	f = open(file_pro_lig, 'r')
	list_pro_lig = f.readlines()
	f.close()

	att="Atom1,Type1,Atom2,Type2,Distance,Interaction\n"
	data=""
	for protein in list_pro_lig:
		#HETATM  - ATOM / ATOM - HETATM
		
		lines=read_file(dir_data, protein.rstrip('\n')) ## directory of ricin ### not tested!!
		f = open(str(protein.rstrip('\n')+".csv"), 'w')
		f.write(att)
		l_line=""
		inter_length=7

		for l in lines:
			#ATOM1
			if l[0][0:6]=='HETATM' and l[2][0:4]=='ATOM':
				for i in range(0,inter_length-2):
					l_line=l_line+l[i]+","
				print (l_line)
				data=data+l_line+l[inter_length-2]+"\n"
				l_line=""
			#ATOM2
			elif l[2][0:6]=='HETATM' and l[0][0:4]=='ATOM':
				for i in range(0,inter_length-2):
					l_line=l_line+l[i]+","
				print (l_line)
				data=data+l_line+l[inter_length-2]+"\n"
				l_line=""
		f.write(data)
		f.close()	
		data=""
		print (str(protein.rstrip('\n')+ " processed \n"))

	print ("Model created!")
	print ("End!")

def create_model_lig_lig(dir_data, file_pro_lig):
	
	print ("Creating model interactions ligand-ligand...")

	#Reading list of contact files
	#file_pro_lig= "list_files.txt"

	f = open(file_pro_lig, 'r')
	list_pro_lig = f.readlines()
	f.close()

	att="Atom1,Type1,Atom2,Type2,Distance,Interaction\n"
	data=""
	for protein in list_pro_lig:
		#HETATM  - ATOM / ATOM - HETATM
		
		lines=read_file(dir_data, protein.rstrip('\n'))
		f = open(str(protein.rstrip('\n')+".csv"), 'w')
		f.write(att)
		l_line=""
		inter_length=7

		for l in lines:
			#ATOM1 and ATOM2 are HETATM
			if l[0][0:6]=='HETATM'and l[2][0:6]=='HETATM':
				for i in range(0,inter_length-2):
					l_line=l_line+l[i]+","
				data=data+l_line+l[inter_length-2]+"\n"
				l_line=""

		f.write(data)
		f.close()	
		
		print (str(protein.rstrip('\n')+ " processed \n"))

	print ("Model created!")
	print ("Bye!")


#Ligplot - test1
#function test
#construct_output("inter-mapped/cluster1/", "1M52:A:P17:119.map.interaction.csv")
#parse_files("inter-mapped/list_cluster1", "cluster1_t", "inter-mapped/cluster1/")
#parse_files("inter-mapped/list_cluster2", "cluster2_t", "inter-mapped/cluster2/")
#parse_files("inter-mapped/list_cluster3", "cluster3_t", "inter-mapped/cluster3/")
#replace_quotation_marks("cluster1_t_G2.json", "cluster1_G2.json")
#replace_quotation_marks("cluster2_t_G2.json", "cluster2_G2.json")
#replace_quotation_marks("cluster3_t_G2.json", "cluster3_G2.json")

##Ligplot -test2
#parse_files("comparacao2/figura1/fig1_1m52-2hyy.txt", "fig1_1m52-2hyy_test", "./comparacao2/figura1/1m52-2hyy/")
#parse_files("comparacao2/figura1/fig1_2hyy.txt", "fig1_2hyy_test", "./comparacao2/figura1/2hyy/")
#parse_files("comparacao2/figura1/fig1_2hyy-2hzi.txt", "fig1_2hyy-2hzi_test", "./comparacao2/figura1/2hzi-2hyy/")
#parse_files("comparacao2/figura1/fig1_2hyy-3cs9.txt", "fig1_2hyy-3cs9_test", "./comparacao2/figura1/3cs9-2hyy/")
#parse_files("comparacao2/figura2/fig2.txt", "fig2_test", "./comparacao2/figura2/")
#parse_files("comparacao2/figura3/fig3.txt", "fig3_test", "./comparacao2/figura3/")
#parse_files("comparacao2/figura4/fig4_1byq-3d36.txt", "fig4_1byq-3d36_test", "./comparacao2/figura4/1byq-3d36/")
#parse_files("comparacao2/figura4/fig4_1z5a-3d36.txt", "fig4_1z5a-3d36_test", "./comparacao2/figura4/1z5a-3d36/")


##Ricin
#parse_files("list_ricin.txt", "ricin_test", "./ricin_data_contacts/")
##new format
#parse_files("./datasets_contacts/list_ricin.txt", "ricin_complex", "./datasets_contacts/ricin_new_format/")

##CDK2
#parse_files("list_files_inter.txt", "cdk2_test_inter", "./cdk2_data_contacts/cdk2_interactions/")
##new format
#parse_files("./datasets_contacts/list_cdk2.txt", "cdk2_complex", "./datasets_contacts/cdk2_new_format/")

##Use this for files that have contacts of any type
#create_model_pro_lig("./ricin/", "list_files_inter.txt")
#create_model_lig_lig("../data_graphs_cdk2/", "list_files.txt")
