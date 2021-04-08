# Reference-free discovery of millions of SNPs permits species and hybrid identification in Carya (Hickory)  
## Robert Literman, Brittany M. Ott, Jun Wen, LJ Grauke, Rachel Schwartz, Sara M. Handy
## Code Walkthrough

### 01) Study Data  
In this manuscript, we apply the SISRS (CITE) bioinformatics pipeline to generate millions of species-identifying single-nucleotide polymorphisms (**SNPs**) for diploid *Carya* (pecan and hickory), using only low-coverage, Illumina short-read data (i.e. genome skims). Samples from this study can be broken down into three groups:  

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
- *C. x laneyi* (*C. ovata* x *C. cordiformis*) x *C. illinoinensis*: n = 2  

3) **Companion Species**: One set of independent, higher-depth sequencing runs for each of the 8 diploid *Carya* species.  

- *Carya aquatica*: SRR6804841  
- *Carya cathayensis*: SRR6784938  
- *Carya cordiformis* SRR6804840  
- *Carya illinoinensis*: SRR6793970  
- *Carya laciniosa*: SRR6804855   
- *Carya myristiciformis*: SRR6804845  
- *Carya ovata*: SRR6804843  
- *Carya palmeri*: SRR6804847   

### 02) Read Preparation  
All reads were trimmed and processed identically.  

1) **Read merging**: Paired-end reads were merged using *bbmerge* from the *BBMap* suite (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/).  

```
bbmerge.sh in=<SAMPLE>_Raw_R1.fastq.gz in2=<SAMPLE>_Raw_R2.fastq.gz outa=<SAMPLE>_Adapters.fa outm=<SAMPLE>_Merged_Raw.fastq.gz outu=<SAMPLE>_Unmerged_Raw_1.fastq.gz outu2=<SAMPLE>_Unmerged_Raw_2.fastq.gz
```

2) **Adapter trimming**: Reads were adapter-trimmed using *bbduk*.  

```
bbduk.sh ref=<SAMPLE>_Adapters_BBMap.fa ktrim=r in=<SAMPLE>_Raw_R1.fastq.gz in2=<SAMPLE>_Raw_R2.fastq.gz out=<SAMPLE>_Adapter_Trim_1.fastq.gz out2=<SAMPLE>_Adapter_Trim_2.fastq.gz
```

3) **Chloroplast assembly**: Adapter-trimmed reads were used to generate chloroplast assemblies for each sample using *getOrganelle* (https://github.com/Kinggerm/GetOrganelle).  

```
get_organelle_from_reads.py -1 <SAMPLE>_Unmerged_Adapter_Trim_1.fastq.gz -2 <SAMPLE>_Unmerged_Adapter_Trim_2.fastq.gz -u <SAMPLE>_Merged_Adapter_Trim.fastq.gz -t 20 -o <DIR>/<SAMPLE> -F embplant_pt -R 50 -k 21,45,65,85,105
```

4) **Quality-trimming**: Reads were quality-trimmed using *bbduk*.  

```
bbduk.sh ref=<SAMPLE>_Adapters_BBMap.fa qtrim=w ktrim=r trimq=10 maq=15 minlength=50 in=<SAMPLE>_Adapter_Trim_1.fastq.gz in2=<SAMPLE>_Adapter_Trim_2.fastq.gz out=<SAMPLE>_Trim_1.fastq.gz out2=<SAMPLE>_Trim_2.fastq.gz
```

5) **Nuclear read separation**: Reads were mapped against the pooled chloroplast assemblies from Step 3 above using *bbmap*, and any reads that mapped were removed from the dataset leaving predominantly nuclear reads.  

```
bbmap.sh in=<SAMPLE>_Trim_1.fastq.gz in2=<SAMPLE>_Trim_2.fastq.gz ambiguous=all threads=20 outm=Chloroplast_Reads/<SAMPLE>_Chloroplast_Trim_1.fastq.gz outm2=Chloroplast_Reads/<SAMPLE>_Chloroplast_Trim_2.fastq.gz outu=Nuclear_Reads/<SAMPLE>_Nuclear_Trim_1.fastq.gz outu2=Nuclear_Reads/<SAMPLE>_Nuclear_Trim_2.fastq.gz
```

### 03) Folder Setup  
If you are trying to follow this Walkthrough as instructional, there is a particular folder structure that needs to be in place.  
- Base_Directory  
    - Reads  
        - TrimReads  
            - Species_1  
                - Sample_A_Trim_1.fastq.gz  
                - Sample_A_Trim_2.fastq.gz  
                - Sample_B_Trim.fastq.gz  
            - Species_2  
                - Sample_C_Trim_1.fastq.gz  
                - Sample_C_Trim_2.fastq.gz  
            - ETC.
        - scripts  
            - Read_Subsetter.py  

### 04) Composite Genome Assembly  
In order to isolate SNPs that can identify *Carya* species, we first need to isolate orthologous loci that were conserved among the *Carya* diploids. While a reference genome for *C. illinoinensis* has been published (https://academic.oup.com/gigascience/article/8/5/giz036/5484800), many clades lack a reference genome and in this study we wanted to provide steps for researchers that have WGS data, but no reasonable reference. SISRS generates orthologous sequence data through a 'composite genome' assembly process (i.e. a pan genome assembled using reads pooled across species). We only used our genome skim data to assembled the composite genome (i.e. no hybrid or companion data was used).  

1) **Read Subsetting**: SISRS composite genomes are assembled by first subsetting reads equally among taxa, and among samples therein. By default, the subsampling targets a final assembly depth of 10X depending on the approximate size of the clade genome. For this study, we used a genome size estimate of 750Mb.  
  
- The subsetting script can be found in [**scripts/SISRS/Read_Subsetter.py**](scripts/SISRS/Read_Subsetter.py), and was run as:  

```
python Read_Subsetter.py -g 750000000
```

2) **Genome Assembly**: SISRS uses *Ray* for genome assembly by default. The subset reads from Step 1 are used in assembly. In this study, both merged and unmerged reads were used, but all reads were treated as single-ended.  

```
mpirun -np 200 Ray -k 31 -p CarAqu_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarAqu_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarAqu_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarAqu_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarAqu_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarAqu_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarAqu_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarAqu_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarAqu_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarCat_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCat_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCat_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarCat_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCat_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCat_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_4_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_4_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_4_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_5_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_5_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_5_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_6_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_6_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_6_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_4_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_4_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_4_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_5_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_5_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_5_Nuclear_Merged_GenomeReads.fastq.gz -p CarLac_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarLac_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarLac_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarLac_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarLac_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarLac_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarMyr_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarMyr_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarMyr_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarOva_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarOva_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarOva_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarOva_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarOva_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarOva_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarOva_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarOva_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarOva_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarPalm_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarPalm_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarPalm_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarPalm_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarPalm_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarPalm_2_Nuclear_Merged_GenomeReads.fastq.gz -o <DIR>/Composite/Ray_Nuclear_Carya
```

3) **Composite Genome Processing**: The resulting contigs assembled by Ray are then prepped for mapping.  

```
mkdir <BASE_DIR>/Composite/Ray_Nuclear_Carya/Composite_Genome
cd <BASE_DIR>/Composite/Ray_Nuclear_Carya/Composite_Genome

# Add SISRS_ to contigs
rename.sh in=<BASE_DIR>/Composite/Ray_Nuclear_Carya/Contigs.fasta out=<BASE_DIR>/Composite/Ray_Nuclear_Carya/Composite_Genome/contigs.fa prefix=SISRS addprefix=t trd=t

# Index contigs
bowtie2-build contigs.fa contigs -p 20
samtools faidx contigs.fa

# Process Ray contig lengths
python <BASE_DIR>/scripts/Genome_SiteLengths.py <BASE_DIR>/Composite/Ray_Nuclear_Carya/Composite_Genome

# Generate BBMap stats
stats.sh in=contigs.fa &> BBmap_Stats
```

