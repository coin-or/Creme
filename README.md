# Crème

Crème is an implementation of the randomized thermal relaxation method to find a feasible solution of the Maximum Feasible Subsystem (MaxFS) problem.
The MaxFS problem consists, given a Linear Programming problem

A **x** &le; **b**,

generally infeasible, in finding a feasible subsystem containing a maximum number of inequalities.

The algorithm is described in the paper 

> E. Amaldi, P. Belotti, R. Hauser, Randomized Relaxation Methods for the Maximum Feasible Subsystem Problem. In M. Jünger, V. Kaibel (ed.), Integer Programming and Combinatorial Optimization, 11th International IPCO Conference, Berlin, Germany, June 8-10, 2005, pages 249-264, 2005.

## Download, installation and usage

Crème is found in the COIN-OR group at GitHub. It can be downloaded with Git. Run the command:

```
git clone https://github.com/coin-or/Creme.git
```

to get the source code of the trunk version (the only version available for now). 
Please read the `INSTALL` file for instructions on how to obtain them.

In general, the master version of Crème is subject to changes that might make it unstable or incorrect at times.
In order to be up-to-date with such changes, you may run the command `git pull` within the `Creme` directory.
A stable version of Crème will be available soon, and will only be changed for bug fixes.
A stable release will be available at <https://github.com/coin-or/Creme/releases>.
This release is not subject to change, and can also be downloaded as a tar ball.

To install Crème, plese refer to general installation [instructions](https://projects.coin-or.org/BuildTools/) for COIN-OR projects.

The impatient may want to issue the following commands:

```
git clone https://github.com/coin-or/Creme.git
cd Creme
mkdir build
cd build
../configure -C
make
make test
make install
```

The above commands place the executable file *creme* in the `Creme/build/bin/` directory, libraries in `Creme/build/lib/`, and include files in `Creme/build/include/`.
An alternative directory can be specified with the `--prefix` option of configure.
For instance, when replacing "`../configure -C`" above with

```
../configure -C --prefix=/usr/local
```

the executable will be installed in `/usr/local/bin/`, the libraries in `/usr/local/lib/`, and the include files in `/usr/local/include/`.
Crème is run as follows:

```
creme myinstance.bz2
```

where `myinstance.bz2` is a bzip2'ed file (`.bz2`) containing number of inequalities and of variables of the infeasible system, coefficients of matrix A, and elements of the right-hand side vector b.
A definition of the file format is available [here](Creme/data/format.txt).

You may also specify a set of options to tweak the performance of Crème. You can visualize them by simply running

```
creme
```

at the command line.


## Documentation

A user manual will be available soon, with explanations on most options available in Crème.
Doxygen documentation is also still unavailable, and when complete it will be generated by running

```
make doxydoc
```

from the same `build/` directory where you ran configure, make, and make
install. Documentation in both html and LaTeX format can be found in
the Doc/ subdirectory. Fire up your browser and take a look at
Doc/html/index.html for documentation of Crème.

## Resources and links

Crème is maintained by [Pietro Belotti](http://myweb.clemson.edu/~pbelott) (creme@list.coin-or.org).

Web page: <https://github.com/coin-or/Creme>

Crème does not have any dependency. If you want to be able to make changes to Crème, download 
[BuildTools](https://github.com/coin-or-tools/BuildTools).

External resources:  [COIN-OR](http://www.coin-or.org), 
[Eclipse Public License](http://www.opensource.org/licenses/eclipse-1.0).

## Improve this page, report a bug, contribute to Crème

As an open-source code, contributions to Crème are welcome. To submit a contribution, please follow the [COIN-OR guidelines](http://www.coin-or.org/contributions.html).

In order to report a bug, use the [issue system](https://github.com/coin-or/Creme/issues).
In order to ensure that your ticket is addressed in a timely fashion, please try to be as exhaustive as you can in the bug report, for instance by reporting what version of Crème you have downloaded and what operating system you are using, and again by attaching the model/data files with which the crash occurred.


## Contributors

 * Pietro Belotti (FICO)


## Project Links

 * Crème executables (will be available in the future)
 * [COIN-OR Initiative](http://www.coin-or.org/)
 * [Crème mailing list](http://list.coin-or.org/mailman/listinfo/Creme)
 * [coin-discuss mailing list](http://list.coin-or.org/mailman/listinfo/coin-discuss)
 * [Report a bug](https://github.com/coin-or/Creme/issues/new)
