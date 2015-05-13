ig_graph
=============

IgGraph is for performing VDJ classification of antibody reads using the colored antibody graph (i.e., a colored de Bruijn graph).

###### Dependencies
* lemon graph library
* boost regex
* SeqAn

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
