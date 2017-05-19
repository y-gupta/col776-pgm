from sys import argv
from os import system
f0 = open(argv[1],'r') # bayesian network file
count = 1 # counts number of queries in output file
treeString = "digraph mentions {\n" # string describing tree in DOT language
numNodes = int(f0.readline())
# add title
treeString += 'node[shape=circle,fixedsize=true,width=0.6]\n' # all nodes should be circles
edges = [] #list of edges in the tree

# Read network and write into treeString
for line in f0:
	[source,dest] = line.split()
	if dest == '[]':
		continue
	dest = dest.strip('[]').split(',')
	for d in dest:
		edges.append((source,d))
f0.close()

# If only network file is provided
if len(argv)==2:
	f1 = open('bn.dot','w')
	for j in range(len(edges)):
		treeString += '"' + str(edges[j][0]) + '"->"' + str(edges[j][1]) + '" \n'
		
	treeString += '}'		
	f1.write(treeString)
	f1.close()
	system("dot -q -Teps bn.dot -o bn.eps") # generate figure from dot file. -Teps tells format of figure. To generate pdf, give argument -Tpdf, -q suppresses warnings
	print "Figure successfully generated..."
	exit()

f1 = open(argv[2],'r') # output file

outLines = f1.readlines() # read output file
numLines = len(outLines)
# Read each pair of lines of output file
for i in range(0,numLines,2):
	colorTree = treeString # colorTree stores string for creating tree with colors in DOT language
	inputList = outLines[i].split()
	
	color_edges = ['' for _ in edges] #stores color for each edge
	f2 = open(str(count)+".dot",'w') # Create dot file for each query
	
	# Create 5 clusters for creating legends in figures, 3 clusters are joined vertically by invisible edges
	colorTree += '\
	subgraph cluster_0 {\
	style=invis\
	b1 [shape=circle,label="query",style="filled",color=white,fixedsize=true,width=0.3];\
	a1 [shape=circle,label="",style=filled,fillcolor=red,width=0.4];\
	a1->b1[constraint=false,style=invis];\
	}\
	subgraph cluster_1 {\
		style=invis\
		b2 [shape=circle,label="observed",style="filled",color=white,fixedsize=true,width=0.3];\
		a2 [shape=circle,label="",style=filled,fillcolor=green,width=0.4];\
		a2->b2[constraint=false,style=invis];\
	}\
	subgraph cluster_2 {\
		style=invis\
		b3 [shape=circle,label="active trail",style="filled",color=white,fixedsize=true,width=0.3,fontsize=10];\
		a3 [shape=rarrow,label="",color=blue,width=0.4,height=0.0];\
		a3->b3[constraint=false,style=invis];\
	}\
	a1->a2[style=invis];\
	a2->a3[style=invis];\
	'
	# Read Xi and Xj and color them red
	node1 = inputList[0] 
	colorTree += '"' + node1 + '"' + '[shape=circle, style=filled, fillcolor=red]\n'
	node2 = inputList[1]
	colorTree += '"' + node2 + '"' + '[shape=circle, style=filled, fillcolor=red]\n'

	#Read evidence and color them green
	evidence = inputList[2].strip('[]')
	if evidence != '':
		evidence = evidence.split(',')
		for e in evidence:
			colorTree += '"' + e + '"' + '[shape=circle, style=filled, fillcolor=green]\n'
	
	# Read output of query
	result = outLines[i+1].strip()
	if result == 'yes':
		colorTree += 'label="No active trail present\n";'
	else:
		activeTrail = result.split()[1].strip('[]').split(',')
		for j in range(len(activeTrail)-1):
			if (activeTrail[j],activeTrail[j+1]) in edges:
				edge_index = edges.index((activeTrail[j],activeTrail[j+1]))
			elif (activeTrail[j+1],activeTrail[j]) in edges:
				edge_index = edges.index((activeTrail[j+1],activeTrail[j]))
			else:
				print "Error !!! In query number "+str(count) + ",Invalid active trail, some edges are not present, program exiting !!!"
				exit()
			color_edges[edge_index] = '[color=blue]' # color the active trail edges blue
	
	for j in range(len(edges)):
		colorTree += '"' + str(edges[j][0]) + '"->"' + str(edges[j][1]) + '" ' + color_edges[j] + '\n'
		
	colorTree += '}'		
	f2.write(colorTree)
	f2.close()
	system("dot -q -Teps "+str(count)+".dot -o "+str(count)+".eps") # generate figure from dot file. -Teps tells format of figure. To generate pdf, give argument -Tpdf, -q suppresses warnings
	#system("rm "+str(count)+".dot") # removes dot files
	print "Figure " + str(count) + " successfully generated..."
	count += 1

