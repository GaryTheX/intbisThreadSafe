# intbisThreadSafe
The original interval analysis code was created while study at Notre Dame under the direction from Prof. Mark Stadtherr and Prof. Joan Brennecke. The code was then modified into multithreading c++ mutex, with a main function under intText.cpp file to performance a bench mark computation through a list of examples from the paper, "Solving Polynomial Systems using a Branch and Prune Approach", by Hentenryck, Mcallester, and Kapur.

The difference in this repository was the use of multithreading c++ mutex to handle generated intervals during interval bisection. With the multithreading schemes, the newly generated intervals would be put onto a idle thread. Then the all threads with given intervals would performance interval anlaysis executions in parallel, while if no thread is in idle, then a newly generated interval would be pushed into stack for later execution. If a thread finished with conclusion (either no solution or single solution), then it would ask for pop from the stack to execute stored intervals.

During the execution of the multithreading schemes, the mutex would be applied to keep the thread safe issue.

Users are welcome to modify the main function to focus on their own interested problems. While on working with their own problems, users need to create new equationset c++ classes (or modifed an existing one such as EquationSet_Robot_kinematics.cpp and hpp files.

"multiReport.txt" and "multiReportSol.txt" are the output files for the build-in benchmark results. 

The whole project should be running under visual studio in windows.
