# generated files
#TEX=1894945_ass1.tex
PLOTSCRIPTS=plotscript1.gpl
RESULTS1=Errors_P_numb_0.000500.txt Errors_P_numb_1.000000.txt Errors_P_numb_10.500000.txt
RESULTS=$(wildcard *.txt)
PROGRAM=adr
OBJS=main.o adr.o sparsematrix.o
PLOTS=ADR.pdf
#REPORT=1894945_ass1.pdf
# additional variables
CPPFLAGS=-std=c++11
#CPPFLAGS=-std=c++11 -g


all: $(PLOTS)

#$(REPORT): $(PLOTS) $(TEX)
	#pdflatex -interaction=batchmode 1894945_ass1.tex
	#pdflatex -interaction=batchmode 1894945_ass1.tex


$(PLOTS): $(RESULTS1) $(PLOTSCRIPTS)
	gnuplot plotscript1.gpl


$(RESULTS1): $(PROGRAM)
	./$(PROGRAM)

#g++ -O3 -Wall -Wfatal-errors -pedantic $(CPPFLAGS) $(OBJS) -o $(PROGRAM)
$(PROGRAM): $(OBJS)
	g++ -Ofast -Wall -Wfatal-errors -pedantic $(CPPFLAGS) $(OBJS) -o $(PROGRAM)


#g++ -O3 -Wall -Wfatal-errors -pedantic $(CPPFLAGS) -c $^ -o $@
$(OBJS): %.o: %.cc
	g++ -Ofast -Wall -Wfatal-errors -pedantic $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(RESULTS)
