####
src/
####

All functions are implemented here.
For the directories, see the corresponding README file.

* array.c

   Function to initialise, destruct, load, and save multi-dimensional arrays with halo cells are implemented.

* common.c

   Utility functions and global parameters are defined.

* config.c

   Environment variables loader which is called when the solver is launched to acquire the runtime-parameters by the user is implemented.

* decide_dt.c

   A function to decide the time step of the simulator is implemented.

* linear_system.c

   Functions to initialise and destruct the tri-diagonal linear system solver are included.

* main.c

   Main function is here.

* save.c

   Functions to save flow field and parameters for restart are included.

* tdm.c

   Kernel functions to solve tri-diagonal matrices are implemented.
