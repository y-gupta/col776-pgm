from sys import argv
import re

f0 = open(argv[1],'r') # output file
count = 1
outLines = f0.readlines() # read output file
numLines = len(outLines)
# Read each pair of lines of output file
for i in range(0,numLines,2):
    query = outLines[i].strip()
    #check if query is in correct format
    m=re.search(r'\d+\s\d+\s\[(\d+)?(,\d+)*\]$',query)
    if not m :
        print "Wrong format of query " + str(count)
        print "Exiting the program !!!"
        exit()
    # Read output of query
    result = outLines[i+1].strip()
    m=re.search(r'^yes$|^no\s\[(\d+,)+\d+\]$',result)
    if not m :
        print "Wrong format of result in query " + str(count)
        print "Exiting the program !!!"
        exit()
    count += 1

f0.close()