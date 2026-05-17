#!/urs/bin/LAB1 python
import sys 
pdbfile=sys.argv[1]
chain=sys.argv[2]

fh=open(pdbfile)

for line in fh:
    if line.startswith('ATOM') and line[21]==chain: print (line.rstrip()) 
