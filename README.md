# ALGO Code Fest #0

# Software Requirements

A modern C compiler such as a recent version of ``gcc`` or ``clang`` is required. The ability to compile against the GNU Scientific Library (v2.5 or later) is also required.

In addition, development tools such as ``make`` (build automation), a debugger, and ``valgrind`` (dynamic analysis) will be very useful, but are not strictly required.

# How to Compile a Problem to a Range of Solvers

1. Compile utilities by running the `Makefile` in the `util` directory.
```bash
[nasf4nio-cs]$ cd src/util/
[util]$ make
```
2. Compile solvers by running the `Makefile` in the `solver` directory.
```bash
[nasf4nio-cs]$ cd src/solver/
[solver]$ make
```

3. Edit the `Makefile` in the problem model directory to select (uncomment) or deselect (comment) the solvers to apply to the problem. In the following example, *GRASP* is selected and *Iterated Greedy with Biased Destruction* is not selected.
```bash
[nasf4nio-cs]$ cd problem/problemX/
[problemX]$ more Makefile
(...)
PROG = problemX-grasp1
(...)
problemX-grasp1:  problemX.o
$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/grasp1.o $(LIBS)
# problemX-bd_ig: problemX.o
# $(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/bd_ig.o $(LIBS)
(...)
```
4. Compile the problem model for the selected solvers and run them. In the following example, *GRASP* is applied to the problem.
```bash
[problemX]$ make
[problemX]$ ./problemX-grasp1 <input data> 0.3 1000
```
