#!/usr/bin/env python
import re, os, sys
from decimal import Decimal
##########################
#Author: GW McElfresh
#Date: 3-23-2015
#Purpose: This program reads the output file of the Cluster-information.txt file and creates an appropriate latex table 
#Usage: ad4tolatex.py your_input_file 

#Error statement
if len(sys.argv) != 2:
        print("Usage: ./ad4tolatextable.py your_input_file")
	sys.exit()

#initalize the regular expressions used to read the Cluster-information.txt file
cluster_number_regex = re.compile(r"""
                        ^CLUSTER \s \# \s* (\d*)
                        """, re.X)
lowest_energy_regex = re.compile(r"""
                        ^Rescored \s Lowest \s Energy \s for \s cluster \s* = \s* (-?\d*\.\d*)
                        """, re.X)
#this will read a conformation from a non-rescored cluster-information.txt
lowest_energy_not_rescored_regex= re.compile(r"""
                        ^\d* \s*           #run number
                        (-?\d*\.\d*) \s*     #ad4 rescore energy
                        N/A \s*     #rescored energy
                        (\d*\.\d*) \s*       #rmsd to reference
                        \d*\.\d* \s*       #rmsd to cluster's lowest energy
                        \s* .*
                        """,re.X)
conformations_in_cluster_regex = re.compile(r"""
                        ^Number \s of \s conformations \s in \s this \s cluster= \s* (\d*)
                        """, re.X)
#this will read a conformation from a rescored cluster-information.txt
rmsd_regex = re.compile(r"""
                        ^\d* \s*           #run number
                        -?\d*\.\d* \s*     #ad4 rescore energy
                        -?\d*\.\d* \s*     #rescored energy
                        (\d*\.\d*) \s*       #rmsd to reference
                        \d*\.\d* \s*       #rmsd to cluster's lowest energy
                        \s* .*
                        """, re.X)


input_file=sys.argv[1]
#open the cluster-information.txt file for reading
cluster_information_txt= open(input_file, 'r')
text_file_lines=cluster_information_txt.readlines()
imax=len(text_file_lines)

i=0
#intitalize a counter for the number of clusters to show in the table
cluster=0
#initialize a flag such that the program knows when to write a line into the latex file
write_a_line="no"
latex_file=open(input_file[0:-4]+".tex",'w')
latex_file.write("\n\\begin{table}\n\\caption{}\n\\label{table:}\n\\begin{center}\n\\begin{tabular}{cccc}\n\\hline\n")
latex_file.write("CLUSTER & LOWEST BINDING& NUMBER IN     & RMSD \\\\\n")
latex_file.write("RANK    & ENERGY        & CLUSTER       &       \\\\\n")
latex_file.write("\\hline\n")

#iterate on the lines of the text file
while i<imax:
    #break after 10 clusters are printed into a table
    if int(cluster)<=10:
        #apply the regexs to the line
        cluster_number= re.search(cluster_number_regex, text_file_lines[i])
        LE_number = re.search(lowest_energy_regex, text_file_lines[i])
        conformation_number = re.search(conformations_in_cluster_regex, text_file_lines[i])
        rmsd_number = re.search(rmsd_regex, text_file_lines[i+2])
        not_rescored_energy= re.search(lowest_energy_not_rescored_regex,text_file_lines[i+2])
        #check if any of the regular expressions find a match
        if cluster_number is not None:
            cluster=cluster_number.group(1)
	    #since the energy location can vary depending on if rescoring is desired, we'll apply the LE regex at the start of each cluster
	    LE_number = re.search(lowest_energy_regex, text_file_lines[i+1])
            #if it matches a rescored file, it grabs the rescored energy
	    if LE_number is not None:
		LE=LE_number.group(1)
                LE=("%.2f" % Decimal(str(round(float(LE),2))))
	    #otherwise, it grabs the ad4 energy of the best ligand in a cluster
	    else:
                not_rescored_energy= re.search(lowest_energy_not_rescored_regex,text_file_lines[i+4])
		LE=not_rescored_energy.group(1)
        if conformation_number is not None:
            conformation=conformation_number.group(1)
            #provide a flag such that the program knows when to write a line into the latex file
            write_a_line="yes"
        if write_a_line=="yes":
            if rmsd_number is not None:  
                rmsd=rmsd_number.group(1)
                rmsd=("%.2f" % Decimal(str(round(float(rmsd),2))))
            if not_rescored_energy is not None:
                rmsd=not_rescored_energy.group(2)
                rmsd=("%.2f" % Decimal(str(round(float(rmsd),2))))
            write_a_line="no"
            latex_file.write(" $"+ str(cluster).rjust(2) + "$" + "\t& $"+ str(LE).rjust(6)+"$"+ "\t& $"+ str(conformation).rjust(4)+"$"+ "\t& $"+ str(rmsd).rjust(6)+"$\t\\\\\n")        i=i+1
    else:
        break

latex_file.write("\end{tabular}\n\end{center}\n\end{table}")

print ("Table is written in the current directory as " + input_file[0:-4]+ ".tex")
