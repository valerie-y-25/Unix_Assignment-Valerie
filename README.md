# UNIX Assignment

## Data Inspection

### Attributes of `fang_et_al_genotypes`

```
$ wc fang_et_al_genotypes.txt
$ awk -F "\t" '{print NF; exit}' fang_et_al_genotypes.txt
$ file fang_et_al_genotypes.txt
$ cut -f 3 fang_et_al_genotypes.txt | sort | uniq -c
```

By inspecting this file I learned that:

* This file has 2,783 lines, 2,744,038 words and 11,051,939 bytes (11M).
* 2782 sample rows and 1 header row. 
* Names of headers include: `Sample_ID`, `JG_OTU`, `Group`, and then genotype samples. 
* `awk` provides column amount. There are 986 : 983 genotype columns and 3 metadata columns
* ASCII text, with very long lines
* last line: I am cutting column 3 from every row since that is the `Group`, then piping it to sort to arrange alphabetically, then piping to `uniq -c` to count how many per group. I saw that there was 290 ZMMIL, 1256 ZMMLR, 27 ZMMMR (1573 total maize) and 900 ZMPBA, 41 ZMPIL,34 ZMPJA (975 total teosinte).


### Attributes of `snp_position.txt`

```
$ wc snp_position.txt
$ awk -F "\t" '{print NF; exit}' snp_position.txt
$ file snp_position.txt
```

By inspecting this file I learned that:

* This file has 984 lines, 13,198 words and 82,763 bytes (81K).
* 983 rows of SNP entries and 1 header row.
* 1st column is SNP ID, 3rd column is chromosome location, 4th column is nucleotide location.
* `awk` shows that there are 15 columns
* ASCII text


## Data Processing

###

### Maize Data

EXTRACT
```
awk -F $'\t' 'NR==1 || $3=="ZMMIL" || $3=="ZMMLR" || $3=="ZMMMR"' original/fang_et_al_genotypes.txt > altered/maize/maize_genotypes.txt
```
* `awk` allows us to look at each row and include it if it's the right conditions.
* -F $'\t' because our files are tab-delimited. -F field delimiter and tells to split lines into columns using tab (\t)
* `'NR==1 || $3=="ZMMIL" || $3=="ZMMLR" || $3=="ZMMMR"` : This is the condition saying I want the maize groups and I use `OR aka ||`because there are three different groups, but all are maize and I want all three. `NR==1` keeps our first header line.
* Reads from the file in original and outputs into altered, maize folder, named `maize_genotypes.txt`
* This gives us groups we want for maize (Group = ZMMIL, ZMMLR, and ZMMMR) and for teosinte (Group = ZMPBA, ZMPIL, and ZMPJA).

TRANSPOSE
```
$ awk -f ../../original/transpose.awk maize_genotypes.txt > transposed_maize_genotypes.txt
```
We were given this! Transposing genotype data so the columns become rows.

SORT AND CUT
```
$ cut -f 1 transposed_maize_genotypes.txt | head
$ tail -n +2 original/snp_position.txt > altered/snp_position_noheader.txt
$ tail -n +4 transposed_maize_genotypes.txt > maize_genotypes_noheader.txt

$ sort -k1,1 maize/maize_genotypes_noheader.txt > maize/maize_genotypes_sorted.txt
$ sort -k1,1 snp_position_noheader.txt > snp_position_sorted.txt

$ cut -f 1,3,4 snp_position_sorted.txt > snp_position_3col_sorted.txt
```
I want to get rid of first header line from snp_position.txt and redirecting it to a new file it incase I mess up and need to go back. (`snp_position_noheader.txt`).<br>
Next, `cut` line shows me that in column 1 and that Sample_ID, JG_OTU, Group as first 3 metadata rows. I get rid of first 3 lines from transposed_maize_genotypes.txt and do the same thing but with teosinte ( `maize_genotypes_noheader.txt` and `teosinte_genotypes_noheader.txt`). In hindsight, noheader was not the best file name as I am removing metadata.<br>
Then, I sort them: `teosinte_genotypes_noheader.txt`, `maize_genotypes_noheader.txt`, and `snp_position_noheader.txt` into `maize_genotypes_sorted.txt`, `snp_position_sorted.txt`, and `teosinte_genotypes_sorted.txt`.<br> 
`sort` on field 1 for all files allows `join` to work properly!<br>
Lastly, I cut out the only columns I needed from `snp_position_sorted.txt`. Assignment indicates that the first column should be "SNP_ID", the second column should be "Chromosome", and the third column should be "Position" . I only need column 1, 3, and 4.

JOIN
```
$ join -1 1 -2 1 -t $'\t' snp_position_3col_sorted.txt maize/maize_genotypes_sorted.txt > maize/maize_joined.txt
```
`-t $'\t'` : This creates tab-delimited output and prevents any issues with white spaces.<br>
`join -1 1 -2 1` means to use column 1 of file 1 and column 1 of file 2 as the key joining column.<br>
This creates files with SNP_ID Chromosome position, and genotypes tacked on.<br>

GENERATING INCREASING AND DECREASING FILES
```
$ awk -F $'\t'  '$2==1' maize/maize_joined.txt | sort -t $'\t'   -k3,3n | sed 's/?\/?/?/g' > ../final/maize/maize_chr1_increase.txt
$ awk -F $'\t'  '$2==1' maize/maize_joined.txt | sort -t $'\t'   -k3,3nr | sed 's/?\/?/-/g' > ../final/maize/maize_chr1_decrease.txt
```
The first line is for chromosomes sorted by increasing position and missing data encoded by ?.
* `$2==1` means column 2 (chromosome) must equal 1 because we are looking at chromosome 1.
* `sort` orders column 3 (position) numerically
* `sed` finds ?/? and replaces it with ?, we have to use \/ because / is a default character must be escaped by \.
The second line is for chromosomes sorted by decreasing position and missing data encoded by -.
* Explanation is similar, except `3nr`, with `r` added on for reverse. 
* `sed` still finds ?/?, but replaces with -
I repeated these for chromosomes 1 through 10 to create 10 increasing and 10 decreasing files.

GENERATING MULTIPLE AND UNKNOWN FILES
```
$ awk -F $'\t'  '$3=="unknown"' maize/maize_joined.txt > ../final/maize/maize_unknown.txt
$ awk -F $'\t'  '$3=="multiple"' maize/maize_joined.txt > ../final/maize/maize_multiple.txt
```
Here, since I’m looking for all the SNPs with unknown (or multiple position), that means I’m looking at column 3 which is the position. `$3=="unknown"` means to find in column 3, those positions that equal “unknown” (or “multiple”). Both are directed to their respective files.

In total, there are 44 files in the `final` folder. 22 maize files in the folder `maize` and 22 teosinte files in the folder `teosinte`.


### Teosinte Data

I did basically the same thing for the Teosinte Data processing so the explanations should apply to all of the teosinte code.

```
$ awk -F $'\t' 'NR==1 || $3=="ZMPBA" || $3=="ZMPIL" || $3=="ZMPJA"' original/fang_et_al_genotypes.txt > altered/teosinte/teosinte_genotypes.txt

$ awk -f ../../original/transpose.awk teosinte_genotypes.txt > transposed_teosinte_genotypes.txt

$ tail -n +4 transposed_teosinte_genotypes.txt > teosinte_genotypes_noheader.txt

$ sort -k1,1 teosinte/teosinte_genotypes_noheader.txt > teosinte/teosinte_genotypes_sorted.txt

$ join -1 1 -2 1 -t $'\t' snp_position_3col_sorted.txt teosinte/teosinte_genotypes_sorted.txt > teosinte/teosinte_joined.txt

$ awk -F $'\t'  '$2==1' teosinte/teosinte_joined.txt | sort -t $'\t'   -k3,3n | sed 's/?\/?/?/g' > ../final/teosinte/teosinte_chr1_increase.txt
$ awk -F $'\t'  '$2==1' teosinte/teosinte_joined.txt | sort -t $'\t'   -k3,3nr | sed 's/?\/?/-/g' > ../final/teosinte/teosinte_chr1_decrease.txt

awk -F $'\t'  '$3=="unknown"' teosinte/teosinte_joined.txt > ../final/teosinte/teosinte_unknown.txt
awk -F $'\t'  '$3=="multiple"' teosinte/teosinte_joined.txt > ../final/teosinte/teosinte_multiple.txt
```

