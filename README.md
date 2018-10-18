# Ginga

<p align="center">
  <img width="250" height="284" src="https://github.com/danielnzg85/2018-EPFL-IGEM/blob/master/Media/Ginga_logo.png">
</p>

Ginga is the first step of the CAPOEIRA to bring one step closer personalised oncology. Ginga is a pipeline base on state of the art genomic tools that are seemlesly integrated in a single streamline pipeline.
Ginga translates whole exome sequencing (WES) data from patients into a library of specific neoantigens with the goal to create a personalised therapy that its effective at targeting uniquely the cancer cells.

Ginga can output multiple metrics that are key to fully treat and diagnose cancer patients. Ginga outputs 3 main patient specific markers:

* A library of patient specific neoantigens --> These are used in CAPOEIRA's vaccine design to target specifically the cancer cells in the patient using the power of the immune system.
* A library of tumor short DNA sequences containing SNPs --> These short sequences can be found in the blood in the form of ctDNA and can be targeted using CAPOERA's detection system to monitor the effectivity of the vaccine. 
* A set of tumor specific chromosomal rearrangements (CR) --> Concentration of CRs found in the ctDNA of the patient can indicate the state of the cancer and if there is relapse.

Ginga's Pipeline:
![](https://github.com/danielnzg85/2018-EPFL-IGEM/blob/master/Media/Tree_Bioinfo_Dani_17-10_SF.png)

## Getting Started

Ginga was developed in Ubuntu 18.04 LTS and it is meant to be deployed in a Linux based system. It uses multiple establish packages and uses algotrithms code in python to .
What you will find in this repository are the python scripts use to streamline the process:

* neoExtract --> Converts the output of the annotated coded amino-acid sequences from Annovar and translate them into short peptides of the form of neoantigens.
* neoSearch --> Ranks the neoantigens according to their binding affinity of MHC-I and retrieves 

### Prerequisites

[Python]()

[libpng 1.6.35](http://www.libpng.org/pub/png/libpng.html)
```
sudo apt-get install libpng-dev
```
[zlib 1.2.11](https://zlib.net/)
```
sudo apt-get install zlib1g-dev
```

[libzip  1.5.1](https://libzip.org/)

```
sudo apt-get install libbz2-dev
```

[libzip  1.5.1](https://libzip.org/)

```
sudo apt-get install libbz2-dev
```


## Packages 

[Trimmomatic 0.38](http://www.usadellab.org/cms/?page=trimmomatic)
[BWA 0.7.17](https://sourceforge.net/projects/bio-bwa/files/)
[Samtools 1.9](http://samtools.sourceforge.net/)
[Picard Tools 2.18.14](https://github.com/broadinstitute/picard/releases/tag/2.18.14)
[GATK 4.0.9.0](https://software.broadinstitute.org/gatk/)
[Annovar 2018Apr16](http://annovar.openbioinformatics.org/en/latest/)
[NetMHC 4.0](http://www.cbs.dtu.dk/services/NetMHC/)




## Authors
* **2018 EPFL IGEM Team**
* **Daniel Nakhaee-Zadeh Gutierrez**


