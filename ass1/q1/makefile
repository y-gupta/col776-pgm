N = 10
K = 4
BN = sample-bn.txt
Q = sample-query.txt
a.out: q1.cpp bayes.h
	g++ -std=c++11 -O2 q1.cpp -o a.out
gen: a.out
	./a.out g "$(N)" "$(K)"
query: a.out
	./a.out q "$(BN)" "$(Q)"

