# Reference-free discovery of millions of SNPs permits species and hybrid identification in Carya (Hickory)  
## Robert Literman, Brittany M. Ott, Jun Wen, LJ Grauke, Rachel Schwartz, Sara M. Handy
## Code Walkthrough

### 01) Study Data  
In this manuscript, we apply the SISRS (CITE) bioinformatics pipeline to generate millions of species-identifying single-nucleotide polymorphisms (**SNPs**) for diploid *Carya* (pecan and hickory), using only low-coverage, Illumina short-read data (i.e. genome skims). Samples from this study can be broken down into three groups:  

1) **Species samples**: Genome skim data for 8 species of diploid *Carya*  

    - *Carya aquatica* (CarAqu; bitter pecan): n = 3  
    - *Carya cathayensis* (CarCat; Chinese hickory): n = 2  
    - *Carya cordiformis* (CarCor; bitternut hickory): n = 5  
    - *Carya illinoinensis* (CarIll; pecan): n = 5  
    - *Carya laciniosa* (CarLac; shellbark hickory): n = 2  
    - *Carya myristiciformis* (CarMyr; nutmeg hickory): n = 1  
    - *Carya ovata* (CarOva; shagbark hickory): n = 3  
    - *Carya palmeri* (CarPalm; Mexican hickory): n = 2  

2) **Hybrid crosses**: Genome skim data for crosses between diploid *Carya* species. There were 5 examples of crosses between pecan (*C. illinoinensis*) and one other diploid *Carya*:  

    - *C. illinoinensis* x *C. aquatica* (xlc): n = 4  
    - *C. illinoinensis* x *C. cordiformis* (xbr; *Carya x brownii*): n = 3  
    - *C. illinoinensis* x *C. laciniosa* (xnuss; *Carya x nussbaumeri*): n = 1  
    - *C. illinoinensis* x *C. myristiciformis* (myrxill): n = 1  
    - *C. illinoinensis* x *C. ovata* (xio; ‘Henke’s Hican’): n = 1  

There were also two crosses where one of the two parents were putative hybrids:  

      - *C. x brownii* x *C. laciniosa* (xbrl): n = 1  
      - *C. x laneyi* (*C. ovata* x *C. cordiformis*) x *C. illinoinensis* (xila): n = 2  

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

**Results**:  
| **Study_ID**  | **Dataset**   | **Organism_Type** | **Sample_Type** | **Total_Trim_Bases** | **Nuclear_Trim_Bases** | **Percent_Nuclear** | **Nuclear_Coverage** |
|---------------|---------------|-------------------|-----------------|----------------------|------------------------|---------------------|----------------------|
| CarAqu        | Study         | Species           | Species         | 2,906,351,845        | 2,819,882,985          | 97%                 | 3.76                 |
| CarCat        | Study         | Species           | Species         | 3,328,074,932        | 3,148,712,883          | 95%                 | 4.20                 |
| CarCor        | Study         | Species           | Species         | 3,513,587,342        | 3,379,188,051          | 96%                 | 4.51                 |
| CarIll        | Study         | Species           | Species         | 4,140,494,431        | 4,061,837,481          | 98%                 | 5.42                 |
| CarLac        | Study         | Species           | Species         | 1,512,372,923        | 1,486,027,288          | 98%                 | 1.98                 |
| CarMyr        | Study         | Species           | Species         | 856,140,514          | 839,708,135            | 98%                 | 1.12                 |
| CarOva        | Study         | Species           | Species         | 2,427,825,202        | 2,370,364,475          | 98%                 | 3.16                 |
| CarPalm       | Study         | Species           | Species         | 1,099,075,171        | 1,084,591,557          | 99%                 | 1.45                 |
| CarAqu        | Companion     | Species           | Species         | 6,316,210,051        | 6,152,118,179          | 97%                 | 8.20                 |
| CarCat        | Companion     | Species           | Species         | 7,776,358,214        | 7,455,973,085          | 96%                 | 9.94                 |
| CarCor        | Companion     | Species           | Species         | 7,711,983,355        | 7,577,602,373          | 98%                 | 10.10                |
| CarIll        | Companion     | Species           | Species         | 5,979,862,916        | 5,864,294,842          | 98%                 | 7.82                 |
| CarLac        | Companion     | Species           | Species         | 4,714,724,473        | 4,636,123,901          | 98%                 | 6.18                 |
| CarMyr        | Companion     | Species           | Species         | 5,211,862,708        | 5,083,260,783          | 98%                 | 6.78                 |
| CarOva        | Companion     | Species           | Species         | 5,401,163,361        | 5,271,766,464          | 98%                 | 7.03                 |
| CarPalm       | Companion     | Species           | Species         | 4,913,799,132        | 4,783,538,578          | 97%                 | 6.38                 |
| CarAqu        | Pooled        | Species           | Species         | 9,222,561,896        | 8,972,001,164          | 97%                 | 11.96                |
| CarCat        | Pooled        | Species           | Species         | 11,104,433,146       | 10,604,685,968         | 95%                 | 14.14                |
| CarCor        | Pooled        | Species           | Species         | 11,225,570,697       | 10,956,790,424         | 98%                 | 14.61                |
| CarIll        | Pooled        | Species           | Species         | 10,120,357,347       | 9,926,132,323          | 98%                 | 13.23                |
| CarLac        | Pooled        | Species           | Species         | 6,227,097,396        | 6,122,151,189          | 98%                 | 8.16                 |
| CarMyr        | Pooled        | Species           | Species         | 6,068,003,222        | 5,922,968,918          | 98%                 | 7.90                 |
| CarOva        | Pooled        | Species           | Species         | 7,828,988,563        | 7,642,130,939          | 98%                 | 10.19                |
| CarPalm       | Pooled        | Species           | Species         | 6,012,874,303        | 5,868,130,135          | 98%                 | 7.82                 |
| CarAqu_1      | Study         | Species           | Specimen        | 1,240,836,872        | 1,216,062,794          | 98%                 | 1.62                 |
| CarAqu_2      | Study         | Species           | Specimen        | 643,386,892          | 629,849,427            | 98%                 | 0.84                 |
| CarAqu_3      | Study         | Species           | Specimen        | 1,022,128,081        | 973,970,764            | 95%                 | 1.30                 |
| CarCat_1      | Study         | Species           | Specimen        | 2,018,656,766        | 1,896,603,905          | 94%                 | 2.53                 |
| CarCat_2      | Study         | Species           | Specimen        | 1,309,418,166        | 1,252,108,978          | 96%                 | 1.67                 |
| CarCor_2      | Study         | Species           | Specimen        | 1,166,373,344        | 1,087,952,554          | 93%                 | 1.45                 |
| CarCor_3      | Study         | Species           | Specimen        | 817,541,788          | 801,973,640            | 98%                 | 1.07                 |
| CarCor_4      | Study         | Species           | Specimen        | 556,963,271          | 546,443,881            | 98%                 | 0.73                 |
| CarCor_5      | Study         | Species           | Specimen        | 534,591,918          | 521,583,259            | 98%                 | 0.70                 |
| CarCor_6      | Study         | Species           | Specimen        | 438,117,021          | 421,234,717            | 96%                 | 0.56                 |
| CarIll_1      | Study         | Species           | Specimen        | 699,631,159          | 670,188,589            | 96%                 | 0.89                 |
| CarIll_2      | Study         | Species           | Specimen        | 1,463,949,510        | 1,435,537,681          | 98%                 | 1.91                 |
| CarIll_3      | Study         | Species           | Specimen        | 895,171,897          | 884,412,347            | 99%                 | 1.18                 |
| CarIll_4      | Study         | Species           | Specimen        | 518,581,163          | 514,102,543            | 99%                 | 0.69                 |
| CarIll_5      | Study         | Species           | Specimen        | 563,160,702          | 557,596,321            | 99%                 | 0.74                 |
| CarLac_1      | Study         | Species           | Specimen        | 1,096,986,140        | 1,079,857,490          | 98%                 | 1.44                 |
| CarLac_3      | Study         | Species           | Specimen        | 415,386,783          | 406,169,798            | 98%                 | 0.54                 |
| CarMyr_1A     | Study         | Species           | Specimen        | 856,140,514          | 839,708,135            | 98%                 | 1.12                 |
| CarOva_1      | Study         | Species           | Specimen        | 809,251,226          | 788,968,773            | 97%                 | 1.05                 |
| CarOva_2      | Study         | Species           | Specimen        | 1,043,384,155        | 1,020,734,890          | 98%                 | 1.36                 |
| CarOva_3      | Study         | Species           | Specimen        | 575,189,821          | 560,660,812            | 97%                 | 0.75                 |
| CarPalm_1     | Study         | Species           | Specimen        | 433,018,606          | 423,807,716            | 98%                 | 0.57                 |
| CarPalm_2     | Study         | Species           | Specimen        | 666,056,565          | 660,783,841            | 99%                 | 0.88                 |
| myrxill_1     | Study         | Hybrid            | Specimen        | 796,277,556          | 783,476,660            | 98%                 | 1.04                 |
| xbr_1         | Study         | Hybrid            | Specimen        | 1,092,545,811        | 1,059,655,584          | 97%                 | 1.41                 |
| xbr_2         | Study         | Hybrid            | Specimen        | 671,977,744          | 662,873,028            | 99%                 | 0.88                 |
| xbr_3         | Study         | Hybrid            | Specimen        | 993,055,538          | 988,392,035            | 100%                | 1.32                 |
| xbrl_1        | Study         | Hybrid            | Specimen        | 931,978,537          | 925,654,760            | 99%                 | 1.23                 |
| xila_1        | Study         | Hybrid            | Specimen        | 956,297,607          | 937,470,941            | 98%                 | 1.25                 |
| xila_2        | Study         | Hybrid            | Specimen        | 659,509,130          | 654,559,563            | 99%                 | 0.87                 |
| xio_1         | Study         | Hybrid            | Specimen        | 777,233,207          | 770,030,411            | 99%                 | 1.03                 |
| xlc_1         | Study         | Hybrid            | Specimen        | 1,019,906,968        | 1,016,113,773          | 100%                | 1.35                 |
| xlc_2         | Study         | Hybrid            | Specimen        | 1,215,089,683        | 1,206,563,163          | 99%                 | 1.61                 |
| xlc_3         | Study         | Hybrid            | Specimen        | 759,828,575          | 751,298,527            | 99%                 | 1.00                 |
| xlc_4         | Study         | Hybrid            | Specimen        | 837,679,541          | 828,246,350            | 99%                 | 1.10                 |
| xnuss_1       | Study         | Hybrid            | Specimen        | 441,618,669          | 427,508,151            | 97%                 | 0.57                 |

### 03) Folder Setup  
If you are trying to follow this Walkthrough as instructional, there is a particular folder structure that needs to be in place. Plans are underway to create a more streamlined, user-friendly version of this pipeline.  
- BASE_DIR  
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

2) **Genome Assembly**: SISRS uses *Ray* for genome assembly by default. The subset reads from Step 1 are used in assembly. In this study, both merged and unmerged reads were used.

```
mpirun -np 200 Ray -k 31 -p CarAqu_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarAqu_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarAqu_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarAqu_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarAqu_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarAqu_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarAqu_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarAqu_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarAqu_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarCat_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCat_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCat_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarCat_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCat_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCat_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_4_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_4_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_4_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_5_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_5_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_5_Nuclear_Merged_GenomeReads.fastq.gz -p CarCor_6_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarCor_6_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarCor_6_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_4_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_4_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_4_Nuclear_Merged_GenomeReads.fastq.gz -p CarIll_5_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarIll_5_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarIll_5_Nuclear_Merged_GenomeReads.fastq.gz -p CarLac_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarLac_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarLac_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarLac_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarLac_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarLac_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarMyr_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarMyr_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarMyr_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarOva_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarOva_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarOva_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarOva_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarOva_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarOva_2_Nuclear_Merged_GenomeReads.fastq.gz -p CarOva_3_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarOva_3_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarOva_3_Nuclear_Merged_GenomeReads.fastq.gz -p CarPalm_1_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarPalm_1_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarPalm_1_Nuclear_Merged_GenomeReads.fastq.gz -p CarPalm_2_Nuclear_Unmerged_GenomeReads_1.fastq.gz CarPalm_2_Nuclear_Unmerged_GenomeReads_2.fastq.gz -s CarPalm_2_Nuclear_Merged_GenomeReads.fastq.gz -o <BASE_DIR>/Composite/Ray_Nuclear_Carya
```

3) **Composite Genome Processing**: The resulting contigs assembled by Ray are then prepped for mapping.  

- The Genome_SiteLengths.py script can be found in [**scripts/SISRS/Genome_SiteLengths.py**](scripts/SISRS/Genome_SiteLengths.py).

```
# Make a new folder for the composite genome
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

4) **Mapping SISRS orthologs onto *C. illinoinensis* reference genome**: While no reference genome or annotation data was used in generating our datasets, we also wanted to get a sense of what we put together. To that end, we mapped SISRS orthologs against the pecan reference genome from NCBI (GCA_011037805.1_ASM1103780v1).  

- The Genome_Mapper.py script can be found in [**scripts/Reference_Genome_Mapping/Genome_Mapper.py**](scripts/Reference_Genome_Mapping/Genome_Mapper.py).

```
# Map SISRS orthologs from C. illinoinensis (pooled) against reference genome

bowtie2 -p 20 -x <DIR>/CarIll_Ref -f -U <DIR>/CarIll.fa -S CarIll.sam

# Extract contigs that map uniquely

samtools view -Su -@ 20 -F 4 CarIll.sam | samtools sort -@ 20 - -o <DIR>/CarIll_Temp.bam

samtools view -@ 20 -H <DIR>/CarIll_Temp.bam > <DIR>/CarIll_Header.sam

samtools view -@ 20 <DIR>/CarIll_Temp.bam | grep -v "XS:" | cat <DIR>/CarIll_Header.sam - | samtools view -@ 20 -b - > <DIR/CarIll.bam

# Create genome mapping file

samtools view <DIR>/CarIll.bam | awk 'BEGIN {OFS = "\t"} { print $1, $3, $4, $2, $6}' > <DIR>/CarIll_MapData.tsv

cut -f1 <DIR>/CarIll_MapData.tsv | sort > <DIR>/Uniquely_Mapping_Contigs

# Create a list of sites from the SISRS orthologs that could be uniquely mapped to the reference genome

python <DIR>/Genome_Mapper.py <DIR>/CarIll_MapData.tsv
```
### Assembly Statistics

| **Contig_Count** | **N50**     | **L50** | **Uniquely_Mapped_Contigs** | **Multiply_Mapped_Contigs** | **Unmapped_Contigs** | **Percent_Uniquely_Mapped** | **Percent_Multiply_Mapped** | **Percent_Unmapped** | **Total_Bases** | **Mapped_Bases** | **Uniquely_Mapped_Bases** | **Percent_Mapped_Bases** | **Percent_Uniquely_Mapped_Bases** | **Percent_Reference_Genome** |
|------------------|-------------|---------|-----------------------------|-----------------------------|----------------------|-----------------------------|-----------------------------|----------------------|-----------------|------------------|---------------------------|--------------------------|-----------------------------------|------------------------------|
| 820,113          | 250,025     | 188     | 475,539                     | 81,970                      | 262,604              | 57.98%                      | 9.99%                       | 32.02%               | 169,107,690     | 100,495,556      | 93,513,132                | 59.43%                   | 55.30%                            | 14.36%                       |