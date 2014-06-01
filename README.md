RJigsaw
=======

Package reproduces toy program's plots as far as I've checked. 

Now works with RootCore. If you've compiled StandAlone, be sure to remove ./Root/*Cint* before compiling with RootCore.

The tool is fairly general - it should work for any number of steps in a decay chain (with a small assumption that at each vertex you have at most one leg that is going to continue to decay) and is separated by hemispheres, so you could in principle do asymmetric topologies.

This can be used StandAlone:

cd cmt/
gmake -f Makefile.Standalone

which will give you the library you can link as shown in the ./run directory. If you

source compileme.sh

You'll compile RJ_ttbar.C (my own version that uses the tool now) linking against the tool library, and you can see the output of 

./RJ_ttbar.out

in output.root. Alternatively (and the real reason I wanted to make such a package...), you can use the tool in a pyROOT script as half-demonstrated in RJ_ttbar.py.
