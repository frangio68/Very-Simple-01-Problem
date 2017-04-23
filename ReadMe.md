ReadMe for the distribution of VerySimple01Problem, a small C++ class for
solving the simplest problem in 0-1 variables. While this is hardly
interesting, the class also "efficiently" enumerates any given subset of all
the (2^n) solutions of the problem in objective function value order. This
has a few possible applications, and it might be useful as a "blueprint" to
someone in need to do the same (enumerate the solutions in objective function
value order) for other combinatorial problems with appropriate structure,
e.g. some path and matching problems. The ideas behind the approach are
fairly extensively described in the interface.

This code is provided "as is", without any explicit or implicit warranty
that it will properly behave or it will suit you needs. Although codes
reaching the distribution phase have usually been extensively tested, we
cannot guarantee that they are absolutely bug-free (who can?). Any use of
the codes is at you own risk: in no case we could be considered liable for
any damage or loss you would eventually suffer, either directly or indirectly,
for having used this code. More details about the non-warranty attached to
this code are available in the license description file.

The code also comes with a "good will only" support: feel free to contact us
for any comments/critics/bug report/request help you may have, we will be
happy to try to answer and help you. But we cannot spend much time solving
your problems, just the time to read a couple of e-mails and send you fast
suggestions if the problem is easily solvable. Apart from that, we can't
offer you any support.

This code is provided free of charge under the "GNU Lesser General Public
License".

Current version is: 1.01

Current date is: May 04, 2012

This release comes out with the following files:

-  doc/: PDF and HTML doxygen documentation

-  doxygen/: doxygen files to produce the documentation

-  VrySmplP/: contains the definition and implementation of the class
   VerySimple01Problem

-  Main/: contains an example of use of the VerySimple01Problem,
   which also works as a correctness tester, with a small makefile

We have used this class in

	http://pages.di.unipi.it/frangio/abstracts.html#MP05a


