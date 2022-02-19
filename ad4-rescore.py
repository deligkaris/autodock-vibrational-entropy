#!/usr/bin/env python
import re,os,sys,math, decimal, subprocess
from copy import deepcopy
from openeye.oechem import *
from decimal import *

#######################################################################################################################
## Author: GW McElfresh and Christos Deligkaris
## Date: 5-10-2015
## Program: AD4_rescore.py
## This work is published on Comput Biol Chem 2018 Jun;74:286-293. doi: 10.1016/j.compbiolchem.2018.03.027
##
#Function: AD4_rescore.py performs 4 main functions:
#          1: Creation of .pdb files from AutoDock4's docking output file (.dlg)
#          2: Creates .pdb files for AutoDock4's proposed ligands that adhere to the symmetry of the macromolecule if desired.
#          3: Calculates the RMSD for each conformation according to OpenEye's OERMSD function (AD4 has symmetry problems that make it untrustworthy)
#          4: Reclustering based on the OERMSD between conformations
#          5: Rescoring based on a method from "Empirical Entropic Contributions in Computational   Docking: Evaluation in APS Reductase Complexes" by Chang et al in J Comput Chem 29: 1753-1761, 2008.
#Usage: Using a .dlg file and the macromolecule and ligand .pdb files used to create the dlg file, appropriately use the arguments below to accomplish reclustering/rescoring 
#
#ARGUMENTS AND FLAGS:
# -ref (example: reference_molecule.pdb) [1 argument is required after the flag] (OPTIONAL)
# -dlg (example: dlg_file.dlg) [1 argument is required after the flag] (REQUIRED
# -outdir is the ouput directory (Example: output/ ) [1 argument is required after the flag] (REQUIRED)
# -flip is a flag that notifies the script that you have macromolecular symmetry [no arguments are required after the flag] (OPTIONAL)
# -conect is a flag that results in using chimera to open and write the PDB files including CONECT information because sometimes the openeye connectivity is too sensitive
#		to the small coordinate changes of the different binding sites
# -rescore is a flag that notifies the script that you want to rescore the clusters [no arguments are required after the flag] (OPTIONAL)
# -macro is the macromolecule or complex that was used in the calculations.  (example: macromolecule.pdb) [1 argument is required after the flag] (REQUIRED)
# -res1 and -res2 are the user defined residues used for matching  (example: 17-18.B) [1 argument is required after each flag] (REQUIRED IF -flip IS USED)
# -rmsd sets the RMSD tolerance for clustering (example: 2.0) [1 argument is required after the flag] (REQUIRED)
#
#This is an example command line:
#AD4_rescore.py -ref 127D_ligand_no_H_oe_fixed.pdbqt -outdir output/ -dlg 127D.dlg -flip -rescore -macro 127D.pdb -res1 17-18.B -res2 5-6.A -rmsd 2.0
#
#######################################################################################################################

#read the command line arguments
i=0
flipped_ligand_flag="no" 	#initalize the flipped ligand flag
rescoring_flag= "no" 		#initalize the rescoring flag
chimera_conect_flag="no"	#initialize the chimera conect flag
reference_path= None
output_path= None
PDB_ID= None
rmsd_tolerance = None
dlg_path= None
residue_1= None
residue_2 = None
while i < len(sys.argv):
	if sys.argv[i]=="-ref":
            	reference_path=sys.argv[i+1]
        elif sys.argv[i]=="-outdir":
            	output_path=sys.argv[i+1]
        elif sys.argv[i]=="-dlg":
            	dlg_path=sys.argv[i+1]
        elif sys.argv[i]=="-flip":
            	flipped_ligand_flag="yes"
        elif sys.argv[i]=="-rescore": 
            	rescoring_flag="yes"
        elif sys.argv[i]=="-macro":
            	PDB_ID= sys.argv[i+1]
        elif sys.argv[i]=="-res1":
            	residue_1 = sys.argv[i+1]
        elif sys.argv[i]=="-res2":
            	residue_2 = sys.argv[i+1]
        elif sys.argv[i]=="-rmsd":
            	rmsd_tolerance = sys.argv[i+1]
	elif sys.argv[i]=="-conect":
		chimera_conect_flag="yes"
        #elif sys.argv[i]=="-help":
        #        print ("\nUsage: Using a .dlg file and ligand .pdb files used to create the dlg file, appropriately use the arguments below to accomplish reclustering/rescoring. Information relating to reclustering/rescoring will be found in output/Cluster-information.txt. PDB files can be found in the output directory \n")
        #        print ("Example of script being used properly:" + "\n" + "AD4_rescore.py -ref ligand.pdb  -outdir output/ -dlg 1D32.dlg -macro 1D32.pdb -res1 1-4.A -res2 5-8.B -rmsd 2.0  -rescore -flip \n \nARGUMENTS: \n -ref is the ligand used in the docking calculation (example: ligand.pdb) [1 argument is required after the flag] (REQUIRED) \n -dlg is the DLG file created by the docking calculation (example: dlg_file.dlg) [1 argument is required after the flag] (REQUIRED) \n -outdir is the ouput directory (Example: output/ ) [1 argument is required after the flag] (REQUIRED) \n -flip is a flag that notifies the script that you have macromolecular symmetry [no arguments are required after the flag] (OPTIONAL) \n -rescore is a flag that notifies the script that you want to rescore the clusters [no arguments are required after the flag] (OPTIONAL) \n -macro is the macromolecule or complex that was used in the calculations.  (example: macromolecule.pdb) [1 argument is required after the flag] (REQUIRED IF -flip IS USED) \n -res1 and -res2 are the user defined residues used for matching  (example: 17-18.B) [1 argument is required after each flag] (REQUIRED IF -flip IS USED) \n -rmsd sets the RMSD tolerance for clustering (example: 2.0) [1 argument is required after the flag] (REQUIRED)\n")
        #        print("Questions and concerns about this code's function may be posted on the Deligkaris group mailing list or emailed to GW McElfresh at mcelfreshgw@gmail.com")
        #        sys.exit()
        i=i+1

##########################################################################################
             
#provide an output file name for the ligand files
output_name="ligand-run-"

##########################################################################################
#ERROR reports depending on misuse of flags

if rmsd_tolerance==None or dlg_path== None or output_path==None:
        if rmsd_tolerance== None:
                print "ERROR: The RMSD tolerance for clustering was not provided. Please type 'ad4_rescore.py -help' for more information."
        if dlg_path==None:
                print "ERROR: The .dlg file was not provided. Please type 'ad4_rescore.py -help' for more information."
        if output_path==None:
                print "ERROR: A directory to store ligand files was not defined. Please type 'ad4_rescore.py -help' for more information."
        sys.exit()
if flipped_ligand_flag=="yes" and PDB_ID.endswith('.pdb')==False: #if a macromolecule is not provided when the user asks to detect macromolecular symmetry
        print "ERROR: The pdb file for the macromolecule used in the .dlg file was not provided. This error was caused because the '-flip' flag was used and either the '-macro' flag was not used, or the argument supplied after -macro was not a pdb file. Please type 'AD4_rescore.py -help' for more information."
        sys.exit()
if flipped_ligand_flag=="yes" and PDB_ID.endswith('.pdb')==True:
        if residue_1 == None:
                print "ERROR: Residue 1 was not provided. When using -flip, you must specify residues that Chimera will match to allow the ligand to use the macromolecule's symmetry. Please type 'ad4_rescore.py -help' for more information."
                sys.exit()
        if residue_2 == None:
                print "ERROR: Residue 2 was not provided. When using -flip, you must specify residues that Chimera will match to allow the ligand to use the macromolecule's symmetry. Please type 'ad4_rescore.py -help' for more information."
                sys.exit()
if not output_path.endswith("/"):
    	print "ERROR: the output destination is a directory. Please put a / character after your desired output folder name."
    	sys.exit()

##########################################################################################

#open the file with the data
inputfile=open(dlg_path,'r')
#read the file
input_line=inputfile.readlines()
pdb_info=dict() 		#maps the run number (j) to the conformation's atom's information
pdb_info_no_hydrogens=dict() 	#has the same purpose as above, but for pdb files without hydrogens
LE_dict=dict() 			#maps the run number (j) to the run's lowest energy
pdb_and_LE=dict() 		#combines atomic info and lowest energy mapped to j, such that it can be written.
pdb_and_LE_no_h=dict() 		#has the same purpose as above, but for pdb files without hydrogens
cluster_population=dict() 	#maps the number of conformations in a cluster to the cluster number (k)
directory=output_path 		#path where the "ligand-run-#" files are written
total_conformations=0 		#counter for the number of conformations generated by AD4
entropic_term=0 		#intializing a variable for rescoring
reranked_energy=0 		#initializing a variable for rescoring for the clusters
rescored_ligand_energy=0 	#initializing a variable for the ranking of ligands within their cluster during reranking
automorph_flag = True           #flag for automorphism during the OERMSD calculations

##########################################################################################

rmsd_dict=dict() 		#initialize a dictionary such that RMSD values can be associated with their run number (in the same fashion as pdb_info)
cluster=[] 			#initialize an empty cluster
cluster_dict=dict() 		#initialize a dictionary such that the contents (ligands) of a cluster can be associated with a cluster's ranking relative to other clusters
flipped_reference_dict=dict() 	#initalize a dictionary such that we can tell if a ligand (by its run number) in a cluster is in its flipped geometry or not.

##########################################################################################

#Physical and Emperical Constants
AD4_correction= 1/0.577 	#a scaling constant used in Chang et. al. that scales their AD4 results. We can use this constant to "un-scale" their constants for our use.  
cluster_weight= AD4_correction*1.148 #an emperical weight defined by Chang et. al. 
AD4_weight=AD4_correction*0.658 #user defined weight to appropriately score molecules
kcal_to_joule = 4184		#conversion from kcal to joules, defined, see front page of Physical Chemistry by McQuarrie and Simon
R = 8.3145 / kcal_to_joule	#page 694 of Mohr2008: CODATA recommended values of the fundamental physical constants: 2006, gas contstant used in the rescoring equation
T=298.15 			#temperature in kelvin used in the rescoring equation, set by autodock

##########################################################################################

#Paths
pdb_no_h_path=output_path + "PDB-without-hydrogens/" #where the unflipped ligand files without hydrogens will be written
flipped_no_h_path=pdb_no_h_path[0:-1] + "-flipped/" #where the flipped ligand files without hydrogens will be written
pdb_with_h_path= output_path + "PDB-with-hydrogens/" #where the unflipped ligand files with hydrogens will be written
flipped_with_h_path=pdb_with_h_path[0:-1] + "-flipped/" #where the flipped ligand files with hydrogens will be written

##########################################################################################

#write the output directories to store the ligands
if not os.path.isdir(output_path):
	print "OUTPUT FOLDER CREATED!"
	os.mkdir(output_path, 0755 )
if not os.path.isdir(pdb_no_h_path):
	os.mkdir(pdb_no_h_path, 0755)
if flipped_ligand_flag == "yes":
        if not os.path.isdir(flipped_no_h_path):
	    	os.mkdir(flipped_no_h_path, 0755)
if not os.path.isdir(pdb_with_h_path):
	os.mkdir(pdb_with_h_path, 0755)
if flipped_ligand_flag == "yes":
        if not os.path.isdir(flipped_with_h_path):
	    	os.mkdir(flipped_with_h_path, 0755)

##########################################################################################

i=0 			#line number counter
imax=len(input_line)
j=0 			#run number counter
k=0 			#cluster number counter

##########################################################################################
# REGULAR EXPRESSIONS

#This regular expression recognizes the run number so the runs can be mapped to output files
run_regex = re.compile(r"""
                        ^DOCKED: \s MODEL
                        \s* (\d*)                #run number out of 100 runs
                        """, re.X)
#This regex recognizes the lowest energy for the best conformation in each run.
LE_regex = re.compile(r"""
                      ^DOCKED: \s USER \s* Estimated \s Free \s Energy \s of \s Binding \s* =
                      \s* (-?\+?\d{1,5}\.\d{1,3})   #lowest energy for the current run
                      \s* kcal/mol \s* .*
                      """, re.X)
#This regex captures information about the conformation's atoms for use in OERMSD
atom_info_regex = re.compile(r"""
                     ^DOCKED: \s ATOM
                     \s* (\d+)                      #atom serial number
                     \s* ([a-zA-Z]+)             #Atomic symbol
                     (\d*)                      #Atom's number
                     \s* (\w*\d*)                      #Residue Name (AltLoc is ignored) 
                     \s* (\w*)                      #ChainID
                     \s\s?\s? (-?\d*?)                      #resSeq
                     \s* (-?\d{1,5}\.\d{1,3})       #x coordinates
                     \s* (-?\d{1,5}\.\d{1,3})       #y coordinates
                     \s* (-?\d{1,5}\.\d{1,3})       #z coordinates
                     \s* \-?\+?\d{1,4}\.\d{1,5}?        #van der Waal interaction
                     \s* \-?\+?\d{1,4}\.\d{1,5}?       #Electrical interaction             
                     \s* \-?\+?\d{1,4}\.\d{1,5}?        #charge on the atom
                     \s* [a-zA-Z]*                        #Type
                     """, re.X)
#This regex captures RMSD values from the output files for clustering
output_regex= re.compile(r"""
			^REMARK\ 300\ OERMSD\ =\ (\d*\.?\d*)
			""", re.X)
file_name_regex= re.compile(r"""
			    ligand-run-(\d*)-without-h.pdb
				
			    """, re.X)

#this regular expression determines if the string/text line is part of a PDB file header
pdb_header_regex = 	re.compile(r"""
				^REMARK	#is there a REMARK at the beginning of the line?
			""", re.X)

##########################################################################################

#This loop iterates over the dlg file (where i is the line number) 
print "READING DLG FILE"
while i<imax:
    #application of the regexs to find information
    run_number= re.search(run_regex, input_line[i])
    LE_number= re.search(LE_regex, input_line[i])
    atom_info= re.search(atom_info_regex, input_line[i])
    #If a regex is satisfied, the appropriate if statement will pull information from the line
    if run_number is not None:
        j=j+1
    if LE_number is not None:
        LE = float(LE_number.group(1))
        if not j in LE_dict: 
            LE_dict[j] = LE
    if atom_info is not None:
        serial_number = atom_info.group(1)
        atomic_symbol = atom_info.group(2)
        atom_ID_number = atom_info.group(3)
        Residue = atom_info.group(4)
        chainID= " " + atom_info.group(5) + " "
        resSeq= atom_info.group(6)
        x_coordinates = atom_info.group(7)
        y_coordinates = atom_info.group(8)
        z_coordinates = atom_info.group(9)
        #These while loops are for formatting the lines in the PDB
        while len(serial_number)==1:
            serial_number = str(" " + serial_number)
        while len(serial_number) < 5:     
            serial_number = " " + serial_number
        while (len(atomic_symbol) + len(atom_ID_number))< 3:
            atom_ID_number=atom_ID_number + " "
        while len(Residue)<3:
            Residue= Residue + " "
        while len(resSeq)<3:
            resSeq= " " + resSeq 
        while len(x_coordinates) < 8:
            x_coordinates = " " + x_coordinates
        while len(y_coordinates) < 8:
            y_coordinates = " " + y_coordinates
        while len(z_coordinates) < 8:
            z_coordinates = " " + z_coordinates
	#store the coordinates & other useful information into a dictionary
        if not j in pdb_info:
            pdb_line = (("ATOM  " + serial_number+ "  " + atomic_symbol+ atom_ID_number + " " + Residue + chainID + resSeq + "    " + x_coordinates + y_coordinates + z_coordinates + "  0.00  0.00     0.000 " + atomic_symbol + "\n"))
            pdb_info[j] = pdb_line
            #if the atom isn't a hydrogen, also store it to create the PDB files without hydrogen
	    if not input_line[i][85] =="H":
		pdb_line_no_h=pdb_line
	        pdb_info_no_hydrogens[j]=pdb_line
        else:
            pdb_line = pdb_line + (("ATOM  " + serial_number+ "  " + atomic_symbol+ atom_ID_number + " " + Residue + chainID + resSeq + "    " + x_coordinates + y_coordinates + z_coordinates + "  0.00  0.00     0.000"+ atomic_symbol.rjust(2) + "\n"))
            pdb_info[j]= pdb_line
	    if not input_line[i][85] =="H":
		pdb_line_no_h= pdb_line_no_h + (("ATOM  " + serial_number+ "  " + atomic_symbol+ atom_ID_number + " " + Residue + chainID + resSeq + "    " + x_coordinates + y_coordinates + z_coordinates + "  0.00  0.00     0.000"+ atomic_symbol.rjust(2) + "\n"))
	        pdb_info_no_hydrogens[j]= pdb_line_no_h
    i=i+1
    jmax=j
    kmax=k

print "READING DLG FILE: DONE"

##########################################################################################

#this variable stores the size of a dictionary that corresponds to the total number of conformations from AD4
total_conformations= float(len(LE_dict))

#this loop joins the coordinate data (pdb_info) and lowest energies (LE_dict) to create the PDB files for each binding site
print ("CREATING PDB FILES")
#reset the run number back to 1, so the the run number can be properly mapped to output files
j=1
while j<=jmax:
    	pdb_and_LE[j]=("REMARK 300 LIGAND RUN NUMBER " + str(j) + "\n" + "REMARK 300 AD4 Energy Score = " + str(LE_dict[j])+ "\n" + pdb_info[j])
    	outfile = open(pdb_with_h_path + output_name + str(j)+ "-with-h.pdb", 'w')
    	text=str(pdb_and_LE[j]+ "END" + "\n")
    	outfile.write(text)
    	outfile.flush()
    	os.fsync(outfile.fileno())
    	if j in pdb_info_no_hydrogens.keys():
		pdb_and_LE_no_h[j]=("REMARK 300 LIGAND RUN NUMBER " + str(j) + "\n" + "REMARK 300 AD4 Energy Score = " + str(LE_dict[j])+ "\n" + pdb_info_no_hydrogens[j])
		outfile = open(pdb_no_h_path+ output_name + str(j)+ "-without-h.pdb", 'w')
        	text=str(pdb_and_LE_no_h[j]+ "END" + "\n")
        	outfile.write(text)
        	outfile.flush()
        	os.fsync(outfile.fileno())
    	j=j+1

print ("CREATING PDB FILES: DONE")

##########################################################################################

#sometimes openeye gets the connectivity (bond orders to be precise) to be different even though
#the molecule is the same with only tiny differences in bond distances
#the distances are different because of rounding when the geometries are written to the DLG file
if chimera_conect_flag=="yes":
	print ("NOW USING CHIMERA TO INCLUDE CONECT INFORMATION IN THE PDB FILES")
	for ligand in os.listdir(pdb_no_h_path):
		#make the chimera command file for opening and saving the PDB files with CONECT information this time
		script_file=open("chimera_script_for_" + ligand[0:-4] +".com", 'w')
        	script_file.write("open "+ pdb_no_h_path + ligand+"\n")
        	#write the PDB file without hydrogens including the CONECT information
        	script_file.write("write format pdb 0 "+ pdb_no_h_path + ligand + "\n")
        	script_file.write("stop")
        	script_file.close()
        	chimera=subprocess.call("chimera" +" --silent" + " --nogui "+"chimera_script_for_" + ligand[0:-4] +".com", shell=True)
        	#delete the script
        	os.remove(script_file.name)
	print ("NOW USING CHIMERA TO INCLUDE CONECT INFORMATION IN THE PDB FILES: DONE")

##########################################################################################

#Make the flipped ligands if the molecule has symmetry
if flipped_ligand_flag=="yes":
    	print "WRITING SYMMETRY-RELATED LIGAND MOLECULES"
    	for ligand in os.listdir(pdb_no_h_path):
        	#make the chimera command file for obtaining a flipped ligand
        	script_file=open("chimera_script_for_" + ligand[0:-4] +".com", 'w')
        	script_file.write("open "+ PDB_ID+"\n")
        	script_file.write("open " + PDB_ID+"\n")
        	script_file.write("open "+ pdb_no_h_path + ligand+"\n")
        	script_file.write("open "+ pdb_with_h_path + ligand[0:-13] + "with-h.pdb" +"\n")
        	script_file.write("match #0:"+residue_1 + " " + "#1:"+ residue_2+"\n")
        	#write the flipped PDB file without hydrogens
        	script_file.write("write format pdb relative 0 2 "+ flipped_no_h_path + "flipped-" + ligand + "\n")
        	#write the flipped PDB file with hydrogens
        	script_file.write("write format pdb relative 0 3 " + flipped_with_h_path+ "flipped-" + ligand[0:-13] + "with-h.pdb" + "\n")
        	script_file.write("stop")
        	script_file.close()
        	chimera=subprocess.call("chimera" +" --silent" + " --nogui "+"chimera_script_for_" + ligand[0:-4] +".com", shell=True)       
        	#delete the script, leaving just the flipped ligand
        	os.remove(script_file.name)    
    	print "WRITING SYMMETRY-RELATED LIGAND MOLECULES: DONE"
    
#This for loop iterates over the created PDB files WITHOUT HYDROGENS and calculates the OERMSD from the reference structure
if not reference_path == None:
	print ("CALCULATING RMSDs FROM LIGAND-REFERENCE")
	for conformation in os.listdir(pdb_no_h_path):
            	#This skips any file in the output folder that isn't a conformation
            	if output_name in conformation:
                	#A regex is applied to find out which run number was iterated on
	        	file_number= re.search(file_name_regex, conformation)
	        	j= int(file_number.group(1))

        		refmolecule=reference_path
        		molecule=pdb_no_h_path+conformation

            		fitfs=oemolistream()
    	        	fitfs.SetFormat(OEFormat_PDB)
    
         		reffs=oemolistream()
          		reffs.SetFormat(OEFormat_PDB)
        		#Modification of the OERMSD function. 
        		if fitfs.open(molecule):
        	    		if reffs.open(refmolecule):
        	        		for fitmol in fitfs.GetOEGraphMols():
						fitc = OEDoubleArray(3* fitmol.GetMaxAtomIdx())
        					fitmol.GetCoords(fitc)
                            			for refmol in reffs.GetOEGraphMols():    
        		
                       	        			rot  = OEDoubleArray(9)
        	        				trans = OEDoubleArray(3)

                					automorph = automorph_flag
                                			heavyOnly = True
                                			overlay = False
			
							refc = OEDoubleArray(3* refmol.GetMaxAtomIdx())
    							refmol.GetCoords(refc)

							#rmsd = OERMSD(refc, fitc, len(refc)/3, False) #use coordinates, not molecules
                					rmsd  = OERMSD(refmol, fitmol, automorph, heavyOnly, overlay, rot, trans)
							if rmsd < 0.:
								print ("ERROR! RMSD calculated has a negative value! Check your PDB files")
								sys.exit()
                                			#Apply the OERMSD function to the flipped ligands if the -flip flag was used, indicating macromolecular symmetry
        						if flipped_ligand_flag=="yes":
                                        			flipped_molecule=flipped_no_h_path+ "flipped-" + conformation
        							fitfs=oemolistream()
            	                        			fitfs.SetFormat(OEFormat_PDB)
    
            	                        			reffs=oemolistream()
            	                        			reffs.SetFormat(OEFormat_PDB)

        							if fitfs.open(flipped_molecule):
        	   			    				if reffs.open(refmolecule):
            	        		        				for fitmol in fitfs.GetOEGraphMols():
											fitc = OEDoubleArray(3* fitmol.GetMaxAtomIdx())
                                                					fitmol.GetCoords(fitc)
                            			    					for refmol in reffs.GetOEGraphMols():    
        		
                										rot  = OEDoubleArray(9)
                										trans = OEDoubleArray(3)
        
                										automorph = automorph_flag
                										heavyOnly = True
                										overlay = False
												
												refc = OEDoubleArray(3* refmol.GetMaxAtomIdx())
                                                        					refmol.GetCoords(refc)

                                                        					#rmsd_flipped = OERMSD(refc, fitc, len(refc)/3, False) #use coordinates, not molecules
                										rmsd_flipped  = OERMSD(refmol, fitmol, automorph, heavyOnly, overlay, rot, trans)
												if rmsd_flipped < 0.:
													print ("ERROR! RMSD calculated has a negative value! Check your PDB files!")
													sys.exit()
                                                       						#Whether the default AD4 ligand or the flipped ligand's rmsd is used is then associated with its run number
                                                        					if rmsd > rmsd_flipped:
                                                                					flipped_reference_dict[j]="yes"
                                                        					else:
                                                                					flipped_reference_dict[j]="no"
                                                        					
                                                        					rmsd=min(rmsd,rmsd_flipped)
                                                        
			        			#store a binding site's current text, such that REMARK lines are all at the top of a file  
			        			conformation_file= open(molecule,'r')
		                			conformation_text= conformation_file.read()
	                        			conformation_file.flush()
			        			os.fsync(conformation_file.fileno())

        						#The rmsd with respect to the reference is written to the conformation's text file
                					outfile=open(molecule, 'w')
                					text=("REMARK 300 OERMSD TO REFERENCE = " + str(rmsd)+ "\n")
                					outfile.write(text)
	        					outfile.write(conformation_text)
                					outfile.flush()
	        					os.fsync(outfile.fileno())

	        					#The OERMSD from the reference is then associated with its run number
	        					rmsd_dict[j]=rmsd

	        					#These are functions from OERMSD that are necessary for alternative RMSD calculations (possibly future work)
                					#OERotate(fitmol, rot)
                					#OETranslate(fitmol, trans)
	print ("CALCULATING RMSDs FROM LIGAND-REFERENCE: DONE")
else:
        for conformation in os.listdir(pdb_no_h_path):
                #This skips any file in the output folder that isn't a conformation
                if output_name in conformation:
                        #A regex is applied to find out which run number was iterated on
                        file_number= re.search(file_name_regex, conformation)
                        j= int(file_number.group(1))
			rmsd_dict[j]=999 #use 999 as a dummy RMSD from the non-given reference

##########################################################################################

k=1
j=1
print "CLUSTERING"

#clustering
deepcopy_LE_dict=dict()
deepcopy_LE_dict=deepcopy(LE_dict)

#Deepcopy is used such that we can iterate over LE_dict and remove entries into it. 
for key in deepcopy_LE_dict.keys():
    #intialize a new cluster
    cluster=[]
    #ensure that the conformation is not already in a cluster. If so, skip it.
    if LE_dict.keys():
        #find the minimum energy across all the remaining conformations
        lowest_run = min(LE_dict, key=LE_dict.get)
        #iterate over the remaining conformations and calculate OERMSD with respect to the minimum energy to find clusters
        j=1
	#intialize the lowest energy ligand in a cluster (flipped if specified)
	refmolecule=pdb_no_h_path+output_name+str(lowest_run)+"-without-h.pdb"
        while j<=jmax:
	    if j in LE_dict.keys():
		if int(j) == int(lowest_run):
		    rmsd=0.000
                #check to see if the molecule has symmetry, if so, calculate the flipped rmsd
                elif(flipped_ligand_flag)=="yes": #perform both OERMSD calculations
                    #OERMSD calculations
            	    molecule=pdb_no_h_path+output_name+str(j)+"-without-h.pdb"
    		
    	    	    fitfs=oemolistream()
    	    	    fitfs.SetFormat(OEFormat_PDB)
    
        	    reffs=oemolistream()
                    reffs.SetFormat(OEFormat_PDB)

                    if fitfs.open(molecule):
	                if reffs.open(refmolecule):
			    for refmol in reffs.GetOEGraphMols():
				refc = OEDoubleArray(3* refmol.GetMaxAtomIdx())
    				refmol.GetCoords(refc)
    	    	                for fitmol in fitfs.GetOEGraphMols():
                            
        		            rot  = OEDoubleArray(9)
                                    trans = OEDoubleArray(3)

         		            automorph = automorph_flag
        		            heavyOnly = True
        		            overlay = False

				    fitc = OEDoubleArray(3* fitmol.GetMaxAtomIdx())
        			    fitmol.GetCoords(fitc)
				    #rmsd = OERMSD(refc, fitc, len(refc)/3, False) #use coordinates, not molecules
         		    	    rmsd  = OERMSD(refmol, fitmol, automorph, heavyOnly, overlay, rot, trans)
				    if rmsd < 0.:
				    	print ("ERROR! RMSD calculated has a negative value! Check your PDB files.")
				    	sys.exit()
		   #Flipped OERMSD calculations
            	    molecule=flipped_no_h_path+"flipped-"+output_name+str(j)+"-without-h.pdb"
    		
    	    	    fitfs=oemolistream()
    	    	    fitfs.SetFormat(OEFormat_PDB)
    
        	    reffs=oemolistream()
                    reffs.SetFormat(OEFormat_PDB)

                    if fitfs.open(molecule):
	                if reffs.open(refmolecule):
			    for refmol in reffs.GetOEGraphMols():
				refc = OEDoubleArray(3* refmol.GetMaxAtomIdx())
                                refmol.GetCoords(refc)
    	    	                for fitmol in fitfs.GetOEGraphMols():
                                                    
        		            rot  = OEDoubleArray(9)
                                    trans = OEDoubleArray(3)

         		            automorph = automorph_flag
        		            heavyOnly = True
        		            overlay = False
				
   	   			    fitc = OEDoubleArray(3* fitmol.GetMaxAtomIdx())
                                    fitmol.GetCoords(fitc)
                                    #rmsd_flipped = OERMSD(refc, fitc, len(refc)/3, False) #use coordinates, not molecules
         		    	    rmsd_flipped  = OERMSD(refmol, fitmol, automorph, heavyOnly, overlay, rot, trans)
				    if rmsd_flipped < 0.:
			            	print ("ERROR! RMSD calculated has a negative value! Check your PDB files.")
				    	sys.exit()
                else:
	            #OERMSD calculations
            	    molecule=pdb_no_h_path+output_name+str(j)+"-without-h.pdb"
    		
    	    	    fitfs=oemolistream()
    	    	    fitfs.SetFormat(OEFormat_PDB)
    
        	    reffs=oemolistream()
                    reffs.SetFormat(OEFormat_PDB)

                    if fitfs.open(molecule):
	                if reffs.open(refmolecule):
			    for refmol in reffs.GetOEGraphMols():
				refc = OEDoubleArray(3* refmol.GetMaxAtomIdx())
                                refmol.GetCoords(refc)
    	    	                for fitmol in fitfs.GetOEGraphMols():
                                                    
        		            rot  = OEDoubleArray(9)
                                    trans = OEDoubleArray(3)

         		            automorph = automorph_flag
        		            heavyOnly = True
        		            overlay = False

		   		    fitc = OEDoubleArray(3* fitmol.GetMaxAtomIdx())
                                    fitmol.GetCoords(fitc)
		 		    #rmsd = OERMSD(refc, fitc, len(refc)/3, False) #use coordinates, not molecules
         		    	    rmsd  = OERMSD(refmol, fitmol, automorph, heavyOnly, overlay, rot, trans)
				    if rmsd < 0.:
				    	print ("ERROR! RMSD calculated has a negative value! Check your PDB files.")
				    	sys.exit()

	        #Determine if the OERMSD with respect to the conformation with the lowest energy satisfies the RMSD tolerance. If not, also check if the flipped OERMSD satisfies the rmsd tolerance. (if -flip was used)
	        if float(rmsd) <= float(rmsd_tolerance):
		    molecule=pdb_no_h_path+output_name+str(j)+"-without-h.pdb"
		    #rounding for formatting
		    relative_rmsd= ("%.2f" % Decimal(str(round(float(rmsd),2))))
		    #if not reference_path == None:
	            #	rmsd_reference=("%.2f" % Decimal(str(round(float(rmsd_dict[j]),2))))
		    #else:
		    #	rmsd_reference = "N/A (no reference provided)"
		    rmsd_reference=("%.2f" % Decimal(str(round(float(rmsd_dict[j]),2))))
                    #Default the status of "Flipped?" to no in case of no macromolecular symmetry
		    flipped = "no "

	            #store a binding site's current text, such that REMARK lines are all at the top of a file  
                    conformation_file= open(molecule,'r')
		    conformation_text= conformation_file.read()
	            conformation_file.flush()
		    os.fsync(conformation_file.fileno())
		    #write the OERMSD and cluster number to the conformation's file
		    
		    outfile=open(molecule, 'w')
                    text=("REMARK 300 OERMSD TO LOWEST IN CLUSTER = " + str(relative_rmsd)+ "\n")
                    outfile.write(text)
                    text=("REMARK 300 CONFORMATION IS IN CLUSTER " + str(k)+ "\n")
	            outfile.write(text)
		    #write the rest of the binding site's text after the remarks
		    outfile.write(conformation_text)
 		    outfile.flush()
		    os.fsync(outfile.fileno())
			    
	            #Store useful information for the conformation as a tuple (change the stored information based on macromolecular symmetry
                    if flipped_ligand_flag == "yes":
	                information_tuple= (j, ("%.2f" % Decimal(str(LE_dict[int(j)]))), rmsd_reference, relative_rmsd, flipped, flipped_reference_dict[j])
                    else:
                    	information_tuple= (j, ("%.2f" % Decimal(str(LE_dict[int(j)]))), rmsd_reference, relative_rmsd)
	            #append the tuple to a list of tuples 
                    cluster.append(information_tuple)
	            #remove the conformation that was just put into a cluster
                    LE_dict.pop(j, None)
                #compare the flipped version of the ligand with the cluster's lowest deltaG ligand
                elif flipped_ligand_flag=="yes":
		    if float(rmsd_flipped) <= float(rmsd_tolerance):
			molecule=pdb_no_h_path+output_name+str(j)+"-without-h.pdb"
		       	#rounding for formatting
	       	        relative_rmsd= ("%.2f" % Decimal(str(round(float(rmsd_flipped),2))))
			rmsd_reference=("%.2f" % Decimal(str(round(float(rmsd_dict[j]),2))))
                        #Store whether or not the RMSD of the flipped ligand is lower than the default ligand or not
		        flipped = "yes"

	                #store a binding site's current text, such that REMARK lines are all at the top of a file of the original pdb files with no hydrogens
                        conformation_file= open(molecule,'r')
		        conformation_text= conformation_file.read()
	                conformation_file.flush()
		        os.fsync(conformation_file.fileno())
		        #write the OERMSD and cluster number to the conformation's file
		    
		        outfile=open(molecule, 'w')
                        text=("REMARK 300 OERMSD TO LOWEST IN CLUSTER = " + str(relative_rmsd)+ "\n")
                        outfile.write(text)
                        text=("REMARK 300 CONFORMATION IS IN CLUSTER " + str(k)+ "\n")
	                outfile.write(text)
		        #write the rest of the binding site's text after the remarks
		        outfile.write(conformation_text)
 		        outfile.flush()
		        os.fsync(outfile.fileno())
			    
	                #Store useful information for the conformation as a tuple (change the stored information based on macromolecular symmetry
                        if flipped_ligand_flag == "yes":
                        	information_tuple= (j, ("%.2f" % Decimal(str(LE_dict[int(j)]))), rmsd_reference, relative_rmsd, flipped, flipped_reference_dict[j])
                        else:
                            	information_tuple= (j, ("%.2f" % Decimal(str(LE_dict[int(j)]))), rmsd_reference, relative_rmsd)
	                #append the tuple to a list of tuples 
                        cluster.append(information_tuple)
	                #remove the conformation that was just put into a cluster
                        LE_dict.pop(j, None)
		    
                j=j+1
            
	    #skip the conformation if it already belongs to a different cluster
	    else:
	        j=j+1
        #This statement stores the list of tuples defined in the RMSD tolerance if statement as a dictionary, such that cluster_dict[i] returns the ith cluster's information 
        cluster_dict[k]=cluster
        k=k+1
        #store kmax for iterating over the clusters
        kmax=k
print "CLUSTERING: DONE!"

##########################################################################################

#initialize a file for the clusters to be written into
outfile=open(output_path+"Cluster-information.txt", 'w')
outfile.flush()
os.fsync(outfile.fileno())

k=1
#This loop writes the cluster information into a single file
print ("WRITING CLUSTER INFORMATION TO OUTPUT FILE")
while k<kmax:
	cluster=cluster_dict[k]
	
	#This skips empty clusters created by excess iterations during the cluster creation
	if cluster == []:
	    k=k+1
	    continue
	#sort the clusters according to lowest energy
        sorted_cluster=cluster_dict[k]
    	sorted_cluster.sort(key=lambda x: float(x[1]))
	#check for rescoring flag
	if str(rescoring_flag) == str("yes"):
		#begin rescoring according to the RMSD method in literature
		conformations_in_cluster=len(cluster_dict[k])
		prob= float((conformations_in_cluster)/(total_conformations))
		LE_in_cluster=sorted_cluster[0][1]
		entropic_term=- float(cluster_weight)*float(R)*float(T)*float((math.log(prob)))
		New_energy=float(AD4_weight)*float(LE_in_cluster) + float(entropic_term)
	else:
	    #this variable allows for new energies to be written as numbers if the Rescoring flag is set to yes, but tells the user if rescoring wasn't performed.
	    New_energy="N/A"	
	
	#write the cluster information for cluster[k]
    	outfile=open(output_path+"Cluster-information.txt", 'a')
    	text=("CLUSTER # " + str(k)+ "\n")
    	outfile.write(text)
    	text=("Rescored Lowest Energy for cluster  = " + str(New_energy)+ "\n")
    	outfile.write(text)
	imax=len(sorted_cluster)
	i=0
	text=(str("Number of conformations in this cluster= ")+ str(imax) + "\n")
	outfile.write(text)
        #Check if the macromolecule has symmetry, and change the columns of the Cluster_information file accordingly
        if flipped_ligand_flag == "yes":
		text=(str("Run number      AD4 Free Energy Score	Rescored Free Energy W*DGvib  RMSD to reference          RMSD to lowest in cluster           Flipped for cluster RMSD?                Flipped for Reference RMSD?")+ "\n")
        else:
        	text=(str("Run number      AD4 Free Energy Score	Rescored Free Energy W*DGvib  RMSD to reference          RMSD to lowest in cluster")+ "\n")
	outfile.write(text)
        #write each of the cluster's members as the tuples created previously (now sorted by energy)
    	while i<imax:
	    if rescoring_flag == "yes":
		rescored_ligand_energy=str("%.2f" % Decimal(str(round(float(AD4_weight)*float(sorted_cluster[i][1])+entropic_term,2)))) + str(" %.2f" % Decimal(str(round(entropic_term,2))))
	    else:
		rescored_ligand_energy= "N/A"
            #if the macromolecule has symmetry, include information about if the ligand is flipped or not (stored in sorted_cluster[i][4])
            if flipped_ligand_flag == "yes":
    	        text=(str(Decimal(sorted_cluster[i][0])) + "	           " + str(Decimal(sorted_cluster[i][1]))+ "		   	" + rescored_ligand_energy + "		        	" + str(Decimal(sorted_cluster[i][2]))+ "		   	"+ str(Decimal(sorted_cluster[i][3]))+ "                              "+ str(sorted_cluster[i][4])+ "                                           " + str(sorted_cluster[i][5]) + "\n")
            else:
		if not reference_path == None:
                	text=(str(Decimal(sorted_cluster[i][0])) + "	           " + str(Decimal(sorted_cluster[i][1]))+ "		   	" + rescored_ligand_energy + "		        	" + str(Decimal(sorted_cluster[i][2]))+ "		   	"+ str(Decimal(sorted_cluster[i][3]))+ "\n")
		else:
			text=(str(Decimal(sorted_cluster[i][0])) + "               " + str(Decimal(sorted_cluster[i][1]))+ "                    " + rescored_ligand_energy + "                          " + "N/A" + "                       "+ str(Decimal(sorted_cluster[i][3]))+ "\n")
      	    outfile.write(text)
	    outfile.flush()
	    os.fsync(outfile.fileno())
	    i=i+1
	k=k+1
print ("WRITING CLUSTER INFORMATION TO OUTPUT FILE: DONE")
cluster_reranking_regex= re.compile(r"""
                          ^Rescored \s Lowest \s Energy \s for \s cluster \s \d* \s = \s 	#the text that is in the line with the rescored energy
			  (-?\d*\.\d*)								#the rescored energy                         
                          """, re.X)

##########################################################################################

print "RERANKING CLUSTERS BASED ON RESCORED FREE ENERGY"

#The file created from the previous loop needs to be reranked according to lowest cluster energy
cluster_file= open(output_path+"Cluster-information.txt",'r')
cluster_line= cluster_file.readlines()
imax=len(cluster_line)
#this list stores the conformation lists that contain text from Cluster_information.txt
clusters_list=[]

i=0
j=0

#This loop iterates over the written Cluster_information.txt file to sort the clusters by their new energy (Thus reranking the clusters)
while i<imax: #the counter i iterates over Cluster_information.txt until it finds the start of a cluster
    #Once a cluster is identified, lists are initalized to store the cluster's text
    if "CLUSTER #" in cluster_line[i]:
	conformation_list=[]
	conformation_and_energy_list=[]
	#the statement then iterates one more line (while using a different counter) 
	j=i+1
	while j<imax:
	    #Each line is the checked for the start of a new cluster
	    #if a new cluster is found, the loop breaks
	    cluster_energy= re.search(cluster_reranking_regex, cluster_line[j])
	    if "CLUSTER #" in cluster_line[j]:
		i=i+1
	        break
	    #otherwise, the line of text is added to the list of binding sites
	    else:
	        conformation_list.append(cluster_line[j])
	    #additionally, the rescored energy from the cluster is stored
            if cluster_energy is not None:
	        reranked_energy=cluster_energy.group(1)
	    j=j+1
	#then, the text containing the cluster's binding sites and the cluster's rescored energy are joined in a 2-d list
	conformation_and_energy_list=(conformation_list, reranked_energy)
	#which is then stored into a list that can be sorted
	clusters_list.append(conformation_and_energy_list)
    i=j

i=0
k=1
imax=len(clusters_list)
#initalize an empty list for sorting the clusters's energies
sorted_cluster_list=[]
cluster_file.flush()
os.fsync(cluster_file.fileno())
#sort the clusters according to lowest energy
sorted_cluster_list=clusters_list
sorted_cluster_list.sort(key=lambda x: float(x[1]))

#truncate the old cluster_information.txt file and write the reranked file
outfile=open(output_path+"Cluster-information.txt",'w')

#This loop writes the clusters that the previous loop reranked according to rescored energy
while i<imax: #iterate over the sorted clusters (i) in cluster list
    jmax=len(sorted_cluster_list[i][0])
    j=0
    #write the current cluster into Cluster_information.txt
    text=("CLUSTER # " + str(k) + "\n")	
    outfile.write(text)
    #iterate over the binding sites of cluster[i] 
    while j<jmax:
	#write the individual binding sites and relevant headers for each cluster
	text=(sorted_cluster_list[i][0][j])
	outfile.write(text)
	j=j+1
    i=i+1
    k=k+1

outfile.flush()
os.fsync(outfile.fileno())

print "RERANKING CLUSTERS BASED ON RESCORED FREE ENERGY: DONE"

##########################################################################################

print ("FORMATTING PDB FILES")

#formatting all other PDB files to have the same header as the PDBs without hydrogens
jmax=len(os.listdir(pdb_no_h_path))
#determine if you need to iterate over the flipped ligands 
if flipped_ligand_flag == "yes":
        list_of_directories=[pdb_with_h_path,flipped_no_h_path, flipped_with_h_path]
else:
        list_of_directories=[pdb_with_h_path]
for directory in list_of_directories:
        j=1
        while j<=jmax:

                #read the appropriate header from the 'without hydrogen' pdbs
	        infile= open(pdb_no_h_path+output_name+str(j) +"-without-h.pdb",'r')
	        full_text= infile.readlines()
		header_text = ""
		i=0
		while i < len(full_text):
			if re.search(pdb_header_regex,full_text[i]) is not None:
				header_text = header_text + full_text[i]
			i=i+1
	        infile.flush()
	        os.fsync(infile.fileno())

                #establish several variables such that this loop can access files regardless of hydrogens or flipping
                if "without-hydrogens" in directory: 
                        file_suffix= "-without-h.pdb"
                elif "with-hydrogens" in directory:
                        file_suffix= "-with-h.pdb"
                if "flipped" in directory:
                        flipped_identifier= "flipped-"
                else:
                        flipped_identifier= ""

                #read the current text of the pdb files in the current directory    
                outfile=open(directory + flipped_identifier + output_name+str(j) + file_suffix, 'r')
	        file_text=outfile.readlines()
	        outfile.flush()
	        os.fsync(outfile.fileno())

                #write the header from the 'without hydrogen' pdbs and then write the text from the 'with hydrogen' pdbs
	        outfile=open(directory + flipped_identifier + output_name+str(j) + file_suffix, 'w')
	        outfile.write(header_text)

                #Fix the Chimera bond issue with multiresidue ligands by ensuring that ligand residues are grouped together for Chimera
		i=0
                sorted_residues=[]
                #this loop adds the lines of the PDB files to a list that stores the atom number, residue number, and then the whole atom line in the first, second, and third element of the list respectively
                while i<(len(file_text)-1):
			if re.search(pdb_header_regex,file_text[i]) is None:      #the header should not be included
	                        sorted_residues.append([file_text[i][7:11], file_text[i][22:26], file_text[i]])
                        i=i+1
                #sort the list first by residue #, then by atom #
		if not sorted_residues[1][1]:
                	sorted_residues.sort(key= lambda x: (float(x[1])))

                #this loop writes the lines of the created list above to create the PDB file with an appropriate header
		i=0
                while i<len(sorted_residues):
		    	if re.search(pdb_header_regex,str(sorted_residues[i][2])) is None:	#the header is included from the 'without hydrogen' pdbs
	            		outfile.write(str(sorted_residues[i][2]))
                    	i=i+1
		#outfile.write("END")
	        outfile.flush()
	        os.fsync(outfile.fileno())

	        j=j+1

##########################################################################################

print ("FORMATTING PDB FILES: DONE")
print ("ALL DONE!")

##########################################################################################
