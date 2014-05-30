#g++ -m32 -Wall RJ_ttbar.C -I../RJigsaw/ -I./ -I$ROOTSYS/include -L$ROOTSYS/lib -L../StandAlone/ `root-config --cflags`  `root-config --glibs`   -lRJigsaw

g++ RJ_ttbar.C  -lRJigsaw  -I$PWD/../RJigsaw/ -I$ROOTSYS/include -L$ROOTSYS/lib -L$PWD/../StandAlone/ `root-config --cflags`  `root-config --glibs`  