# Reference-free discovery of millions of SNPs permits species and hybrid identification in Carya (Hickory)  
## Robert Literman, Brittany M. Ott, Jun Wen, LJ Grauke, Rachel Schwartz, Sara M. Handy
## Code Walkthrough

### 1) Study Data  
In this manuscript, we apply the SISRS (CITE) bioinformatics pipeline to generate millions of species-identifying SNPs for diploid *Carya* (pecan and hickory), using only low-coverage, Illumina short-read data (i.e. genome skims). Samples from this study can be broken down into three groups:  

1) **Species samples**: Genome skim data for 8 species of diploid *Carya*  

- *Carya aquatica* (bitter pecan): n = 3  
- *Carya cathayensis* (Chinese hickory): n = 2  
- *Carya cordiformis* (bitternut hickory): n = 5  
- *Carya illinoinensis* (pecan): n = 5  
- *Carya laciniosa* (shellbark hickory): n = 2  
- *Carya myristiciformis* (nutmeg hickory): n = 1  
- *Carya ovata* (shagbark hickory): n = 3  
- *Carya palmeri* (Mexican hickory): n = 2  

2) **Hybrid crosses**: Genome skim data for crosses between diploid *Carya* species. There were 5 examples of crosses between pecan (*C. illinoinensis*) and one other diploid *Carya*:  

- *C. illinoinensis* x *C. aquatica*: n = 4  
- *C. illinoinensis* x *C. cordiformis* (*Carya x brownii*): n = 3  
- *C. illinoinensis* x *C. laciniosa* (*Carya x nussbaumeri*): n = 1  
- *C. illinoinensis* x *C. myristiciformis*: n = 1  
- *C. illinoinensis* x *C. ovata* (‘Henke’s Hican’): n = 1  

There were also two crosses where one of the two parents were putative hybrids:  

- *C. x brownii* x *C. laciniosa*: n = 1  
- (*C. ovata* x *C. cordiformis*) x *C. illinoinensis*: n = 2  

3) **Companion Species**: One set of independent, higher-depth sequencing runs for each of the 8 diploid *Carya* species.  

- *Carya aquatica*: SRR6804841  
- *Carya cathayensis*: SRR6784938  
- *Carya cordiformis* SRR6804840  
- *Carya illinoinensis*: SRR6793970  
- *Carya laciniosa*: SRR6804855   
- *Carya myristiciformis*: SRR6804845  
- *Carya ovata*: SRR6804843  
- *Carya palmeri*: SRR6804847   
