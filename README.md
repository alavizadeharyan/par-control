This project computes the solution of a control problem using parallel CPUs.
The control problem is governed by a Schrodinger-type equation.
Several  algorithms  can be used to solve  such  problems, one of which  is
"monotonic" algorithm.
By  dividing the time  interval into  smaller intervals, and  sending these
problems to  different  CPUs, each  CPU  solves the provided  problem using
monotonic  algorithm and returns the obtained  solution to the main CPU. By
simply concatenating these solutions, an approximation of the solution will
be achieved after several iterations.

For more details, please read the pdf  file "Using parallel CPUs to solve a
control problem".
