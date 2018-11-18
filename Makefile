# generated files
TEX=u1894945-assignment2.tex
PLOTSCRIPTS=errors.gpl sol1.gpl sol2.gpl sol3.gpl sol4.gpl sol5.gpl
RESULTS=Errors_P_numb_0.000500.txt Errors_P_numb_1.000000.txt Errors_P_numb_10.500000.txt P_numb_0.000500_Solution_J_9.txt P_numb_0.000500_Solution_J_19.txt P_numb_0.000500_Solution_J_39.txt P_numb_0.000500_Solution_J_79.txt P_numb_0.000500_Solution_J_159.txt P_numb_0.000500_Solution_J_319.txt P_numb_0.000500_Solution_J_639.txt P_numb_1.000000_Solution_J_9.txt P_numb_1.000000_Solution_J_19.txt P_numb_1.000000_Solution_J_39.txt P_numb_1.000000_Solution_J_79.txt P_numb_1.000000_Solution_J_159.txt P_numb_1.000000_Solution_J_319.txt P_numb_1.000000_Solution_J_639.txt P_numb_10.500000_Solution_J_9.txt P_numb_10.500000_Solution_J_19.txt P_numb_10.500000_Solution_J_39.txt P_numb_10.500000_Solution_J_79.txt P_numb_10.500000_Solution_J_159.txt P_numb_10.500000_Solution_J_319.txt P_numb_10.500000_Solution_J_639.txt P_numb_0.010000_Solution_J_9.txt P_numb_0.010000_Solution_J_19.txt P_numb_0.010000_Solution_J_39.txt P_numb_0.010000_Solution_J_79.txt P_numb_0.010000_Solution_J_159.txt P_numb_0.010000_Solution_J_319.txt P_numb_0.010000_Solution_J_639.txt P_numb_0.000000_Solution_J_9.txt P_numb_0.000000_Solution_J_19.txt P_numb_0.000000_Solution_J_39.txt P_numb_0.000000_Solution_J_79.txt P_numb_0.000000_Solution_J_159.txt P_numb_0.000000_Solution_J_319.txt P_numb_0.000000_Solution_J_639.txt P_numb_0.000500_Residual_J_9.txt P_numb_0.000500_Residual_J_19.txt P_numb_0.000500_Residual_J_39.txt P_numb_0.000500_Residual_J_79.txt P_numb_0.000500_Residual_J_159.txt P_numb_0.000500_Residual_J_319.txt P_numb_0.000500_Residual_J_639.txt P_numb_1.000000_Residual_J_9.txt P_numb_1.000000_Residual_J_19.txt P_numb_1.000000_Residual_J_39.txt P_numb_1.000000_Residual_J_79.txt P_numb_1.000000_Residual_J_159.txt P_numb_1.000000_Residual_J_319.txt P_numb_1.000000_Residual_J_639.txt P_numb_10.500000_Residual_J_9.txt P_numb_10.500000_Residual_J_19.txt P_numb_10.500000_Residual_J_39.txt P_numb_10.500000_Residual_J_79.txt P_numb_10.500000_Residual_J_159.txt P_numb_10.500000_Residual_J_319.txt P_numb_10.500000_Residual_J_639.txt P_numb_0.010000_Residual_J_9.txt P_numb_0.010000_Residual_J_19.txt P_numb_0.010000_Residual_J_39.txt P_numb_0.010000_Residual_J_79.txt P_numb_0.010000_Residual_J_159.txt P_numb_0.010000_Residual_J_319.txt P_numb_0.010000_Residual_J_639.txt P_numb_0.000000_Residual_J_9.txt P_numb_0.000000_Residual_J_19.txt P_numb_0.000000_Residual_J_39.txt P_numb_0.000000_Residual_J_79.txt P_numb_0.000000_Residual_J_159.txt P_numb_0.000000_Residual_J_319.txt P_numb_0.000000_Residual_J_639.txt
PROGRAM=adr
OBJS=main.o adr.o sparsematrix.o
PLOTS=ADRer.pdf ADRsol1.pdf ADRsol2.pdf ADRsol3.pdf ADRsol4.pdf ADRsol5.pdf
REPORT=u1894945-assignment2.pdf
# additional variables
CPPFLAGS=-std=c++11
#CPPFLAGS=-std=c++11 -g


all: $(REPORT)

$(REPORT): $(PLOTS) $(TEX)
	pdflatex -interaction=batchmode u1894945-assignment2.tex
	pdflatex -interaction=batchmode u1894945-assignment2.tex


$(PLOTS): $(RESULTS) $(PLOTSCRIPTS)
	gnuplot errors.gpl
	gnuplot sol1.gpl
	gnuplot sol2.gpl
	gnuplot sol3.gpl
	gnuplot sol4.gpl
	gnuplot sol5.gpl

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
