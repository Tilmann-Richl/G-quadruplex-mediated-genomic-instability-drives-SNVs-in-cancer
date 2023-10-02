# G-quadruplex-mediated genomic instability drives SNVs in cancer
This repository contains the computational scripts used in the manuscript **G-quadruplex-mediated genomic instability drives SNVs in cancer** by Richl et al. 


## GDC Download
Prepare files fro donwload from https://portal.gdc.cancer.gov/repository by adding them to cart and downloading the **Manifest, metadata and biospecimen (json)** files. Then, download the individual patient files in parallel: 
```bash
#!/bin/bash
# Set variables
set -a
cpus=$(ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w)
# Define download function using the gdc-client
DOWNLOAD (){
    ID=$1
    ./gdc-client download $ID
    echo "Downloaded $ID"
    }
# Invoke function with GNU parallel
parallel --timeout 200% -j $cpus -a GDC_FILES "DOWNLOAD $(echo {})"
```
## Preparing genomic windows
Then, genomic windows are prepared for the subsequent analysis: 
```bash
# Make 5 kbps windows in hg38 exons
bedtools makewindows -b gene_exons_hg38.bed -w 5000 | awk '{print $1,$2,$3,NR}' OFS="\t" > exons_hg38_5kb_windows
```
## Processing GDC Data 
The data is pre-processed, formatted and joined into one file holding all cSNVs from all patients. Additionally, each cSNV is assigned to the respective patient it originates from. 


```bash
# Unzip all files
tar -xvf *.tar 
gunzip *
# Organize all *.maf files in a list
ls *.maf > list
# Count Deletions in patient and write to file
while IFS=" " read -r line
do 
    DEL=$(grep "DEL" $line | wc | awk '{print $1}')
    COUNT=$(wc $line | awk '{print $1}')
    RATIO=$(echo "$DEL / $COUNT" | bc -l)
    echo -e "$line\t$COUNT\t$DEL\t$RATIO" >> OUTPUT_FILE
done < list

# Write mutations with patient ID to one file
while IFS=" " read -r line
do 
    awk '{print $0,FILENAME;}' OFS="\t" $line | grep -v "#" | grep -v Hugo >> out
    echo $line
done < list

# Sort and rename

sort -k1,1 -k2,2n out > ALL_GDC_sorted.bed

```
Then, cSNVs are intersected with 5 kpbs genomic windows to assign each cSNV to the genomic window it is contained within. Likewise, G4s are intersected with genomic windows. 
```bash
# Intersect windows with GDC cSNvs and G4s
bedtools intersect -a ALL_GDC_sorted.bed -b exons_hg38_5kb_windows -wa -wb > exons_hg38_5kb_windows_cSNVs
bedtools intersect -a Na_PDS_G4s.hg38 -b exons_hg38_5kb_windows -wa -wb > exons_hg38_5kb_windows_G4s
# Count unique intersects per window for G4s and cSNVs
cut -f 8 exons_hg38_5kb_windows_G4s | sort | uniq -c > exons_hg38_5kb_windows_G4s.count
cut -f 19 exons_hg38_5kb_windows_SNPs | sort | uniq -c > exons_hg38_5kb_windows_SNPs.count

```
## Analyse the local density of G4s vs. randomized G4s 
```bash
# Generate 50 Mio random G4s with a length of 200 bp
bedtools random -l 200 -n 50000000 -g gene_exons_hg38.bed | sort -k1,1 -k2,2n > random_exon_G4s
# Find closest cSNVs to each random G4 and G4 and write distances to file
bedtools closest -a random_exon_G4s -b ALL_GDC_sorted.bed -D a > random_exon_SNPs_closest
bedtools closest -a G4s -b ALL_GDC_sorted.bed -D a > SNPs_closest
# Count unique occurences of each Distance
cut -f 22 random_SNPs_closest | sort | uniq -c | sort -r > random_exon_SNPs_closest.count
cut -f 22 SNPs_closest | sort | uniq -c | sort -r > SNPs_closest.count
```
## Analyze the local density of G4s and cSNVs in CDK4 and MYC genes
```bash
# Extract genomic positions of coding exons from the fasta file
grep ">" MYC.fa | awk '{print $2}' | sed 's/range=//g' | sed 's/:/\t/g' | sed 's/-/\t/g' | awk '{print $0, "ex" NR}' OFS="\t" > MYC_cds.bed
# Intersect coding exons with G4s 
bedtools intersect -b Na_PDS_G4s.hg38 -a MYC_cds.bed -c > MYC_cds_G4s.bed
# Intersect with cSNVs
bedtools intersect -a MYC_cds_G4s.bed -b ALL_GDC.sorted.bed -c > MYC_G4s_cSNVs.tsv
```
