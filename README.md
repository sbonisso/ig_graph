ig_graph
=============

IgGraph is for performing VDJ classification of antibody reads using the colored antibody graph (i.e., a colored de Bruijn graph).

###### Dependencies
* lemon graph library
* boost regex
* SeqAn
* gcc with C++11 support [GCC 4.7 and later](https://gcc.gnu.org/projects/cxx0x.html)

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
