ig_graph
=============

IgGraph is for performing VDJ classification of antibody reads using the colored antibody graph (i.e., a colored de Bruijn graph).

###### Dependencies
* [SeqAn 2.0+](http://packages.seqan.de/)
* gcc with C++11 support [GCC 4.7 and later](https://gcc.gnu.org/projects/cxx0x.html)
* [TCLAP](http://tclap.sourceforge.net/)
* [cereal](https://github.com/USCiLab/cereal)
* [cpptest](http://cpptest.sourceforge.net/)

###### Compilation

Compile by running make in the top-level directory. 

```
cd ig_graph
make
```

To create tests, additionally run make test.

```
make test
./unit_test
```
