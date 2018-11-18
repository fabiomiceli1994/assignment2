################################################################################

                     C1 - ASSIGNMENT 2 - ID:1894945
                              README file

################################################################################

To execute the program type "make" (without " ") into the terminal and press
enter. The tests executed by the program, indexed through J (number of points in
the interior of the domain), alpha, beta, gamma and the Pe or Da number, will
appear printed on the terminal while the program runs. Different .txt files,
containing the residual error as function of the iterations required for
convergence, the approximate solution and the final errors for the different
cases tested will be printed. The report will be produced by the Makefile as
well. The scripts required by gnuplot for the construction of the plots are
provided and called by the Makefile: the different plots produced are saved as
pdf files.
If other tests are to be performed, it is sufficient to change the values of
the parameters vector (Jvec, alpha, beta, gamma) in the main file. A suitable
message is printed in case the parameters are not changed correctly. The Program
automatically change its output in case the paramters inserted define the other
type equation defined in the assignment instructions. 
