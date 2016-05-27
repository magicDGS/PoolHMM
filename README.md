Pool-HMM
========

[Pool-HMM](https://qgsp.jouy.inra.fr/index.php?option=com_content&view=article&id=56&Itemid=63)
is a program to estimate allele frequencies and detecting selective sweeps in NGS pools described 
in [Boitard _et al._ 2012](http://mbe.oxfordjournals.org/content/29/9/2177.long) and implemented 
in [Boitard _et al._ 2013](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12063/abstract). See more information
in the [program README](https://raw.githubusercontent.com/magicDGS/PoolHMM/master/README.txt).

This repository contains the code downloaded from [here](https://forge-dga.jouy.inra.fr/projects/pool-hmm)
where I included few modifications for support parallelization in MacOS computers and more features. I used the portable
Queue implemented in [LEMON](https://github.com/vterron/lemon/tree/9ca6b4b1212228dbd4f69b88aaf88b12952d7d6f)
and add the dependencies in the source files in the original code.

I do not own any of the code in this repository although my name is included in the modified files.
Pool-HMM is under a [CeCILL-B license](http://www.cecill.info/) and the portable Queue under a 
[GNU v.3 license](http://www.gnu.org/licenses/).


## Extra features

* Parallelization support in MacOS
* Support for .pileup.gz input files