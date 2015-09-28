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

IgGraph has been tested on many different versions of Ubuntu, Mac OS X, and should work on most unix-based OSs. If building from source the dependencies listed are contained as subtrees in the git repository, so there is no need to install them separately.

###### Compilation

Clone the repository, all dependencies are contained as subtrees, so nothing special needs to be done. Alternatively, you can also download the tar/zip and build from those.

```
git clone https://github.com/sbonisso/ig_graph.git
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
./run_iggraph_hsap.sh --read_file igseq_reads.fa \
		      --output_file igseq_vdj.tsv \
		      --num_thread 4
```

The options and flags can be seen when running the built-in help:

```
./iggraph -h
```
