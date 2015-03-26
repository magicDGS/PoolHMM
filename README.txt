=========================
       pool-hmm 
=========================

This program aims at estimating allele frequencies and detecting selective sweeps using NGS data from a sample of pooled individuals from a single population.
It implements the HMM method of Boitard et al (2012). In this approach, each polymorphic site on the genome is assumed to have a hidden state,
which can take one of the three following values : "Neutral", "Intermediate" and "Selection". These three values are associated with different patterns of allele frequencies.
Using the observed NGS data, Pool-hmm predicts the most likely hidden state at each site.
Any window of sites with hidden state "Selection" is then considered as a sweep signal.

The method is based on the infinite sites model, i.e. at each genomic position only two alleles can exist, the ancestral allele and the derived allele. 
By default, the method assumes that the ancestral allele is unknown. A folded frequency spectrum is computed, and the allele count estimated with --estim is the minor allele count. Two alternative strategies can be chosen using the option --ancestral-allele (or -a). If argument 'reference' is used,
the method assumes that the ancestral allele is the reference allele provided in the pileup (third column). 
If argument 'provided' is used, the user needs to provide the ancestral allele by adding a 7th column to the pileup. 
Ancestral character states can be A, C, G, T or N (for unknown). In both cases ('reference' and 'provided'),
an unfolded allele frequency spectrum is computed. At all sites where the ancestral allele is unknown, the allele count estimated with --estim is the minor allele count.

Needs :

	- This program needs a version of Python older than 2.5 and strictly earlier than 3.
	- NumPy and SciPy have to be installed (you can download Numpy and SciPy from here : http://new.scipy.org/download.html)

Infile :

	- .pileup : the NGS data in pileup format (samtools), as produced by the "samtools mpileup" command without -u or -g option.
	See the test data sets for an example.

Usage :

	- Open a terminal and go into the directory where the Python code is stored.
	- Put the pileup file in this same directory.
	- Execute the command "python pool-hmm.py" followed by a list of input parameters.
	- The definition and usage of all input parameters can be obtained by executing the command "python pool-hmm.py -h"

Parameters :
	
	General parameters :

	-f (or --prefix) prefix, where prefix.pileup is the name of the pileup file to be analyzed.
        -n (or --n) value, where value is the number of chromosomes in the pool.
        -R (or --region) chr:start..end, where chr is the name of a chromosome and start and end are positions on the reference genome ; 
	this specifies the region to be analyzed in the pileup file ; by default all sites are used.
	-a (or --ancestral-allele) value, where value can be 'unknown', 'reference' or 'provided' ; if 'provided' is used, 
	the pileup file must have an additional 7th column providing the ancestral allele ; default is 'unknown'.
	-e (or --quality-encoding) name, where name is the format of quality scores, 'illumina' (default) or 'sanger'.
	-P (or --nb-process) <value>, number of processes used by the program.

	Choice of analysis :

	-o (or --only-spectrum) ; if this option is used, the allele frequency spectrum is estimated and stored in the .spectrum file, but there is no prediction of the hidden 	states.
	-p (or --pred) ; if this option is used, the program predicts the hidden states and generates the .emit, .pred, .stat and .post files.
	If the pileup file covers several chromosomes, these chromosomes must be analyzed independently (i.e. one command --pred for each chromosome, using the data filter 
	--region described below) 
	-S (or --estim) ; if this option is used, the program estimates the derived allele counts in the sample and generates the .estim file.

	Data filters :

	-c (or --min-coverage) value, where value is the minimum coverage required for a site to be used in the analyses ; default is 1. 
        -C (or --max-coverage) value, where value is the maximum coverage required for a site to be used in the analyses ; default is 10000.
        -q (or --min-quality) value, where value is the minimum phred score required for a nucleotide to be used in the analyses ; default is 0.

	Specific options allele frequency spectrum estimation :

	-r (or --ratio) value, where 1/value is the proportion of sites that are used to estimate the frequency spectrum ; default is 10.
        -t (or --theta) where value is a starting value for theta=4N*mu (N effective population size, mu per site mutation rate) for the EM algorithm estimating the frequency 		spectrum ; default is 0.005.
	
        Specific options for sweep prediction or allele frequency estimation :

	-k (or --k) value, where value is the per site transition probability between hidden states.
	-E (or --emit-file) ; if this option is used, the emission probabilities of the HMM are not computed, but taken from the file prefix.segemit 
	(or prefix_chr:start..end.emit if the --region option is used) that has been generated by a previous execution of Pool-hmm.
	-s (or --spectrum-file) prefix2, where prefix2.spectrum provides a frequency spectrum (typically, but not necessarily, computed from earlier runs of the program) ; 
	if this option is used, the allele frequency spectrum is not estimated but taken from this file.       

Outfiles :

	- .spectrum : the allele frequency spectrum estimated from the sample. Vector of size (n+1), where n is the number of chromosomes 
	in the sample. This file can serve as an input in subsequent analyses.
	- .segemit : emission probabilities of the HMM. Array of size p * 5, where p is the number of observed polymrphic sites in the input file.
	Each line represents one observed polymorphic site. The two first columns indicate the chromosome and position of this site.
	The three following columns represent the emission probabilities in hidden states "Neutral", "Intermediate" and "Selection".
	This file can serve as an input in subsequent analyses.
	- .pred : predicted hidden states. Array of size p * 2. Each line represents one observed polymorphic site. The first column indicate the position of this site.
	The second column provides the predicted hidden state (1="Neutral",2="Intermediate",3="Selection") at this site.
	- .stat : summary of the predicted sweep windows. The first line indicates the number of detected sweep windows.
	Then each line represents one sweep window. The first column indicates the first genomic position within the sweep. 
	The second column indicates the last genomic position within the sweep. The value in the last column is prob2=-log10(1-prob),
	where prob is the maximum of the posterior probability of hidden state "Selection" within the window. 
	This value can be used to sort sweep windows (the larger prob2, the more significant the sweep window).
	The transformation of prob into prob2 is used because prob is generally very close to 1.
	- .post : posterior probabilities of hidden state "Selection". Vector of size p. Each line represents one observed polymorphic site. 
	The first column indicates the position of this site. The second column provides the posterior probability of hidden state "Selection" at this site.
	- .estim : maximum a posteriori estimates of the number of derived alleles in the sample (from 0 to n). Array of size p * 5. 
	Each line represents one genomic site covered by reads. The two first columns indicate the chromosome and position of this site.
	For options 'ancestral' and 'provided', columns 3 and 4 indicate the ancestral and derived alleles at this site, and the last column provides the estimated derived allele count.
	For option 'unknown', columns 3 and 4 indicate the major and minor alleles at this site, and the last column provides the estimated minor allele count.
	This latter format is also used for option 'provided', but only at positions where the ancestral allele is unknown (N).

	Information about the progress of Pool-hmm is also printed to the standard output during the execution. 

Example :

	The file test_droso.pileup was obtained from a pool of 194 drosophila haplotypes (97 flies). To estimate the allele frequency spectrum from this data, the command is:
	>> python pool-hmm.py --prefix test_droso -n 194 --only-spectrum --theta 0.005
	To detect selective sweeps, the command is (with a possibly different value of -k) :
	>> python pool-hmm.py --prefix test_droso -n 194 --pred -k 0.001 --theta 0.005
	This command also estimates the allele frequency spectrum, because it is necessary for sweep detection.
	If this spectrum has already been estimated (using pool-hmm or any other method) and is stored in test_droso.spectrum, one can avoid computing it again and directly 		predict sweeps using :
	>> python pool-hmm.py --prefix test_droso -n 194 --pred -k 0.001 --spectrum-file test_droso
	If the result of this first prediction is not satisfactory, additional predictions can be tried with different values of --k. 
	In this case, starting from the emission probabilities of the first analysis generally saves a lot of time. The command to do so is :
	>> python pool-hmm.py --prefix test_droso -n 194 --pred -k 0.01 --emit-file
	Finally, an estimation of the allele frequency at each site can be obtained by :
	>> python pool-hmm.py --prefix test_droso -n 194 --estim
	or
	>> python pool-hmm.py --prefix test_droso -n 194 --estim --spectrum-file test_droso
	(the second command being only possible if test_droso.spectrum already exists)
	
	The files test_droo.spectrum, test_droso.segemit, test_droso.pred (for -k 0.001), test_droso.post (for -k 0.001), test_droso.stat (for -k 0.001) and test_droso.estim 
	resulting from these commands are provided. Due to the random selection of sites when computing the allele frequency spectrum, each new execution of the commands
	might lead to slightly different results.

Parallelization :

	In order to speed up the computations, the program can be run simultaneously on several processes. The number of processes 
	specified in the command line is not bounded, but the number of processes that can effectively work simultaneously is obviously determined 
	by the architecture of the computer where the program is run.
	Note also that if the result of the product (10000 * number of processes) is greater than the number of lines in the pileup file,
	the efficiency of the parallelization will be reduced. In practice, this means that parallelization is only useful for sufficiently 
	large files.

Tuning parameters :

	Parameter --r determines the number of genomic positions that are used when estimating the allele frequency spectrum. If L is the number of genomic positions covered by the NGS data,
	only L/r genomic positions are used (r>=1) for the estimation. Small values of r provide a more accurate estimation (i.e. with smaller variance), but imply a longer computation time.
	To reduce the computation time while keeping a very accurate estimation, one can follow the quick rule of thumb r <= p*L/100,
	where p is the smallest probability in the allele frequency spectrum. This probability is of course unknown before the analysis, but is generally of order 10e-4 or 10e-5.
	It can also be roughly evaluated by a first execution of pool-hmm based on a small number of sites (i.e. a large value of r).
	The above rule of thumb ensures that the standard deviation in the estimation of p (and consequently of other probabilities in the allele frequency spectrum)
	is at least ten times smaller then the estimation itself.

	Parameter --k determines the sensitivity of sweep detection. Larger values of k result in a larger number of sweeps detected, including both true sweeps signals and false positives.
	To control the rate of false positives, parameter k should be calibrated using neutral simulations, as described and discussed in Boitard et al (2009, 2012).
	An easier alternative is to simply increase (decrase) k if the number of sweeps detected is considered too small (large).
	Although this empirical approach does not allow to quantify the statistical evidence of the sweep signals, it can be very useful to point out interesting regions of the genome 
	in the perspective of further analyses.
	
===========================================
   filter-pileup-by-feature.py
===========================================

This script filters a pileup file for a genomic feature, which is provided in a .gtf file.


Usage :

	- Open a terminal and go into the directory where the Python code is stored.
	- Obtain a pileup-file and a .gtf file containing the annotation of the species of interest
	- Execute the command "python filter-pileup-by-feature.py" followed by a list of input parameters.
	- The definition and usage of all input parameters can be obtained by executing the command "python filter-pileup-by-feature.py -h"


Parameters :

	'--pileup' the pileup file which should be filtered
	'--gtf' a file containing the annotation for the given species in the gtf format 
                (may for example be obtained from the UCSC genome browser: http://genome.ucsc.edu/ section Tables)
        '--feature' the genomic feature that should be filtered (e.g: exon, intron, CDS), has to be present in the gtf file
        '--output' a pileup file filtered for the given feature


Example :

	The folder testdata contains an example of a pileup file and of a gtf file.
        Following an example of how the given pileup file may be filtered for exons:
        >> python filter-pileup-by-feature.py --pileup testdata/test1.pileup --gtf testdata/annotation.gtf --feature exon --output testdata/exon.pileup

        Similarly the pileup file may be filterd for CDS or introns:
        >> python filter-pileup-by-feature.py --pileup testdata/test1.pileup --gtf testdata/annotation.gtf --feature CDS --output testdata/cds.pileup
        >> python filter-pileup-by-feature.py --pileup testdata/test1.pileup --gtf testdata/annotation.gtf --feature intron --output testdata/intron.pileup
 

=======

Licence :

	Copyright Simon Boitard, David Robelin, Robert Kofler.

	e-mail : simon.boitard@toulouse.inra.fr, david.robelin@toulouse.inra.fr, robert.kofler@vetmeduni.ac.at.

	This software is governed by the CeCILL-B license under French law and
	abiding by the rules of distribution of free software.  You can  use, 
	modify and/ or redistribute the software under the terms of the CeCILL-B
	license as circulated by CEA, CNRS and INRIA at the following URL
	http://www.cecill.info. 

	As a counterpart to the access to the source code and rights to copy,
	modify and redistribute granted by the license, users are provided only
	with a limited warranty  and the software's author,  the holder of the
	economic rights,  and the successive licensors  have only  limited
	liability. 

	In this respect, the user's attention is drawn to the risks associated
	with loading,  using,  modifying and/or developing or reproducing the
	software by the user in light of its specific status of free software,
	that may mean  that it is complicated to manipulate,  and  that  also
	therefore means  that it is reserved for developers  and  experienced
	professionals having in-depth computer knowledge. Users are therefore
	encouraged to load and test the software's suitability as regards their
	requirements in conditions enabling the security of their systems and/or 
	data to be ensured and,  more generally, to use and operate it in the 
	same conditions as regards security. 

	The fact that you are presently reading this means that you have had
	knowledge of the CeCILL-B license and that you accept its terms.
