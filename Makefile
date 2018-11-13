# generated files
#TEX=1894945_ass1.tex
#PLOTSCRIPTS=plotscript_delta.gpl plotscript_lambda.gpl
RESULTS= Test.dat
PROGRAM=adr
OBJS=main.o adr.o sparsematrix.o
#PLOTS=Gauss_Seidel_delta.pdf Gauss_Seidel_lambda.pdf
#REPORT=1894945_ass1.pdf
# additional variables
CPPFLAGS=-std=c++11
#CPPFLAGS=-std=c++11 -g


all: $(RESULTS)

#$(REPORT): $(PLOTS) $(TEX)
	#pdflatex -interaction=batchmode 1894945_ass1.tex
	#pdflatex -interaction=batchmode 1894945_ass1.tex


#$(PLOTS): $(RESULTS) $(PLOTSCRIPTS)
	#gnuplot plotscript_delta.gpl
	#gnuplot plotscript_lambda.gpl


$(RESULTS): $(PROGRAM)
	./$(PROGRAM)

#g++ -O3 -Wall -Wfatal-errors -pedantic $(CPPFLAGS) $(OBJS) -o $(PROGRAM)
$(PROGRAM): $(OBJS)
	g++ -Ofast -Wall -Wfatal-errors -pedantic $(CPPFLAGS) $(OBJS) -o $(PROGRAM)


#g++ -O3 -Wall -Wfatal-errors -pedantic $(CPPFLAGS) -c $^ -o $@
$(OBJS): %.o: %.cc
	g++ -Ofast -Wall -Wfatal-errors -pedantic $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(RESULTS)
