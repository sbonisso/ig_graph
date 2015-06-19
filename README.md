ig_graph
=============

IgGraph is for performing VDJ classification of antibody reads using the colored antibody graph (i.e., a colored de Bruijn graph).

###### Dependencies
* gcc with C++11 support [GCC 4.7 and later](https://gcc.gnu.org/projects/cxx0x.html)
* [SeqAn 2.0+](http://packages.seqan.de/)
* [TCLAP](http://tclap.sourceforge.net/)
* [cereal](https://github.com/USCiLab/cereal)
* [cpptest](http://cpptest.sourceforge.net/)

###### Compilation

Clone the repository using the --recursive option to initialize each submodule dependency:

```
git clone --recursive git@bitbucket.org:sbonisso/ig_graph
```

Compile by running make in the top-level directory. 

```
cd ig_graph
make
```

To create tests, additionally run make test. To create the unit_test executable, you must install cpptest independently (e.g., on Ubuntu: apt-get install libcpptest0).

```
make test
./unit_test
```

