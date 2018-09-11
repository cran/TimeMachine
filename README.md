TimeMachine
============

The TimeMachine is a computationally efficient R implementation of the
[simulation approach for stochastic trees][paper] by the same name proposed by
Jasra, De Iorio, and Chadeau-Hyam.

In general, the simulation of genealogical trees backwards in time (from
observations up to the most recent common ancestor) is hindered by the fact
that, while approaching the root of the tree, coalescent events become rarer,
with a corresponding increase in computation time.
The Time Machine tackles this issue by stopping the simulation of
the tree before reaching the MRCA and correcting for the induced bias.

Installation
------------
The latest stable version of the package is available from [CRAN][cran], and can
be easily installed using `install.packages`.


[paper]: http://dx.doi.org/10.1098/rspa.2010.0497
[cran]: https://cran.r-project.org/package=TimeMachine

