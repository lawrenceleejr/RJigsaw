#g++ -m32 -Wall RJ_ttbar.C -I../RJigsaw/ -I./ -I$ROOTSYS/include -L$ROOTSYS/lib -L../StandAlone/ `root-config --cflags`  `root-config --glibs`   -lRJigsaw

g++ -o RJ_ttbar.out RJ_ttbar.C  -lRJigsaw  -I$PWD/../RJigsaw/ -I$ROOTSYS/include -L$ROOTSYS/lib -L$PWD/../StandAlone/ `root-config --cflags`  `root-config --glibs`  $PWD/../StandAlone/libRJigsaw.so

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/../StandAlone/