# Introduction
Microbial Index of Pathogenic bacteria (MIP) is an easy-to-use bioinformatic package that calculates the quantitative microbial pathogenic risk. MIP assesses the disease risk to human from a microbiome sample and reports the quantitative indices from multiple aspects using a reference pathogen database and an artificially curated pathogen-disease interaction network. MIP works as a plug-in tools of Parallel-META 3, and supports 16S rRNA amplicon sequences as input.

# Software Requirement and Dependency
```
MIP requires Parallel-META 3, please refer to 
https://github.com/qibebt-bioinfo/parallel-meta#installation-guide
for installation.
```

# Installation Guide
## MIP provides a fully automatic installer for easy installation.
**a. Download the package**
```
git clone https://github.com/XXXX.git	
```

**b. Install by installer**
```
cd mip
source install.sh
```

The package should take less than 1 minute to install on a computer with the specifications recommended above.

The example dataset could be found at “example” folder. Check the “example/Readme” for details about the demo run.

#  Basic Usage
With a input 16S rRNA amplicon sequence file, e.g. sample1.fasta:
**a. Profiling by Parallel-META using MIP database**
```
PM-parallel-meta -r sample1.fasta -D P -f F -o sample1.out

The “sample1.out” folder is the profiling result.
```

**b. Parse out the microbial index of pathogenic bacteria**
```
PM-parse-mip -i sample1.out/classification.txt -o sample1.mip
```

The output contains 5 files：
```
sample1.mip.summary.out: The overall MIP;
sample1.mip.taxa.out: The relative abundance of species contributed to the overall MIP (MIP LV1);
sample1.mip.OTU.Abd.out: The relative abundance of reference OTUs contributed to the overall MIP (MIP LV1);
sample1.mip.infection.out: Human diseases associated with the identified pathogenic bacteria (MIP LV2);
sample1.mip.site.out: Targeted human organs or body sites (MIP LV3).
```

#  Batch Processing
MIP also supports the batch input of profiling results by the following alternative two forms (compatible with Parallel-META 3):
**a. Sample list**
```
PM-parse-mip -l samples.list -o samples.mip
```

in which parameter “-l” assigns the file list of profiling results of multiple samples.
The format of a sample list:
```

Sample1	/home/data/sample1.out/classification.txt
Sample2	/home/data/sample2.out/classification.txt
...	
SampleN	/home/data/sampleN.out/classification.txt
```

**b. Abundance tables**
```
PM-parse-mip -T samples.OTU.Abd -o samples.mip
```

in which parameter “-T” assigns the profiling result of OTU table of multiple samples. The format of a OTU table:
```

	OTU_1	OTU_2	OTU_3	…	OTU_M
Sample1	100	200	0	…	50
Sample2	0	300	600	…	100
…	…	…	…	…	60
SampleN	50	80	0	…	200
```

# Example Dataset
Here we provide a demo dataset with profiling results of 20 microbiome. To run the demo, you can either run the script “Readme”:
```
cd example
sh Readme
```

or type the following command:
```
PM-parse-mip -l samples.list -o samples.mip
```