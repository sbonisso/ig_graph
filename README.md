ig_graph
=============

IgGraph performs VDJ classification of antibody reads using the colored antibody graph (i.e., a colored de Bruijn graph).

[http://sbonisso.github.io/ig_graph/](http://sbonisso.github.io/ig_graph/)

###### Dependencies
* gcc with C++11 support [GCC 4.7 and later](https://gcc.gnu.org/projects/cxx0x.html)
* [SeqAn 2.0+](http://packages.seqan.de/)
* [TCLAP](http://tclap.sourceforge.net/)
* [cereal](https://github.com/USCiLab/cereal)
* [catch](https://github.com/philsquared/Catch)

###### Compilation

Clone the repository using the --recursive option to initialize each submodule dependency:

```
git clone --recursive https://github.com/sbonisso/ig_graph.git
```

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

###### Use

The iggraph executable can be run directly, or for simplicity, you can use the two run_iggraph_{hsap,mmus}.sh scripts, from inside the ig_graph directory, to run for either human or mouse.

```
./run_iggraph_hsap.sh --read_file igseq_reads.fa --output_file igseq_vdj.tsv --num_thread 4
```

The options and flags can be seen when running the built-in help:

```
./iggraph -h
```
