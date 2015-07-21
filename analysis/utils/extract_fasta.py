import glob

pdblist = glob.glob('./*0001.pdb')

outfileA = open("ChainA_seqs.fasta", "w")
outfileB = open("ChainB_seqs.fasta", "w")
AA3toAA1= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

for pdb in pdblist:
	pdbfile = open(pdb, 'r')

	pdbnum = pdb.split("_")[2]

	chain_A = {}
	chain_B = {}



	for line in pdbfile:
		chain = "None"
		#print line
		if line[:4] == "ATOM":
			chain = line[21]
			#print line
		if chain in "AB":
			resnum = int(line[23:26].strip())
			threeLetRes = line[17:20]
			oneLetRes = AA3toAA1[threeLetRes]

			if chain == "A":
				this_Chain = chain_A
			elif chain == "B":
				this_Chain = chain_B

			this_Chain[resnum] =  oneLetRes

	A_string = ""
	B_string = ""


	#A_max = 150
	#B_max = 150
	try:
		A_min = min(chain_A.keys()) -1
		A_max = max(chain_A.keys())
	except ValueError:
		print "error"
		print pdb

	try:
		B_min = min(chain_B.keys()) -1
		B_max = max(chain_B.keys())
	except ValueError:
		B_max = 0
		B_min = 0


	for pos in range(A_min,A_max):
		pos = pos+1
		try:
			A_string =  A_string + chain_A[pos]
		except KeyError:
			A_string =  A_string + "-"
	A_string =  A_string + "\n"

	for pos in range(B_min,B_max):
		pos = pos+1
		try:
			B_string =  B_string + chain_B[pos]
		except KeyError:
			B_string =  B_string + "-"
	B_string =  B_string + "\n"


	accessionline = str(">%s\n") % pdbnum

	outfileA.write(accessionline)
	outfileA.write(A_string)

	outfileB.write(accessionline)
	outfileB.write(B_string)







	pdbfile.close()


outfileA.close()
outfileB.close()

