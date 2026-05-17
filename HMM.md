## PIPELINE: DEVELOPMENT OF A HIDDEN MARKOV MODEL PROFILE FOR ANNOTATION OF KUNITZ-TYPE PROTEINS DOMAIN 

To build a Hidden Markov Model (HMM) profile for the Kunitz-type protease inhibitor domain, was used the advanced query builder in the Protein Data Bank, we filtered entries according to the following criteria:
- Sequence length between 40 and 80 amino acids
- PFAM identifier: PF00014
- Structure resolution ≤ 3 Å
The filtered search returned 135 entries. Before downloading the results, we generated a custom tabular report (available in this repository) containing the Entry ID, Sequence, and Chain identifier.
The search specifically targeted PF00014, but many structures contained multiple chains. Thus, in homodimeric entries, identical chains (e.g., Chains A and B) were reduced to a single representative to avoid sequence redundancy. In heterodimeric complexes containing auxiliary fragments (e.g., Chains Q or J), only the primary Kunitz monodomain was retained based on sequence length, in order to produce a non-redundant, high-resolution dataset suitable for HMM construction.

---
### ASSESSMENT OF THE TRAINING SET 
 After inspecting the tabular report and selecting only the relevant fields from the downloaded file: Entry ID, chain (auth asym ID), and sequence; quotation marks (") were removed to standardize the dataset prior to further processing.

```
awk -F "," '{if ($1!= "") print $1,$3,$2}' rcsb_pdb_custom_report_20260415204951.csv | tr -d '"' 
```

The filtered table was converted into FASTA format, exported as a .txt file, and further refined by retaining only sequences between 40 and 80 residues in length:
```
awk -F "," '{if ($1 != "" && length($2)<80 && length($2)>40) print $1,$3,$2}' rcsb_pdb_custom_report_20260415204951.csv | tr -d '"' | awk '{print ">"$1$3;print $2}'> pdb_seqs.txt
```

After filtering were retained a total of 102 sequences.
```
grep ">" pdb_seqs.txt | wc
```
102     102     714

---

### Removing redundancy with MMseqs2
Clustering was performed by aploading the FASTA file to the MMseqs2 web server, selecting the following parameters: 
- Sequence identity threshold: 0.9
- Coverage threshold: 0.9
- Input sequences: 102 -> Reduced representative set: 19

Paste the clustering output obtained from MMseqs2 into a file.txt, by using _vi_ terminal text editor. One representative sequence was extracted for each cluster and  corresponding IDs were stored in a new file after removing associated chains:
```
grep -A 1 ^Clus pdb_seqs.clust | grep ">" | tr -d ">" > pdb_id_chain.rep
```
Create a new file whithout chains identifiers
```
cut -c 1-4 pdb_id_chain.rep > pdb_id.rep
```
---
### PERFORMING THE STRUCTURAL ALIGNMENT WITH MTMALIGN 
To perform the structural alignment the representative PDB structures were first downloaded from PDB and A dedicated directory was created to store the extracted chain structures ```mkdir pdb_chains```: 

```
for i in $(cat pdb_id.rep); do wget "http://files.rcsb.org/view/$i.pdb"; done
```

Since the list of representative structures contained both the PDB identifier and the chain identifier in a single string, these values were separated into two columns using:
```
awk '{print substr($1,1,4) "\t" substr($1,5)}' pdb_id_chain.rep > pdb_chains.rep
```

To isolate the chains of interest, was used a Python script: get_chain.py 
Since the structural coordinates in PDB files are stored in lines beginning with the ATOM record, and the chain identifier is located at column 22, the script selectively extracted the desired chain from each structure.

```
while IFS=$'\t' read -r pdb chain; do python get_chain.py $pdb".pdb" $chain > "$pdb$chain.pdb"; done < pdb_chains.rep
```
In this loop, the first column was assigned to the variable pdb and the second to chain. The resulting files contained only the structural coordinates of the selected Kunitz-domain chains.
Finally, all extracted chain files were listed into an input file for mTM-align:

```
ls pdb_chains/* > my_input.txt
```

After an initial structural alignment, two structures (5YW1A and 1D0DA) were excluded from the dataset because they introduced inconsistencies in the alignment and a clened input file was generated: 

```
sed -i '' '/5YW1A/d' my_clean_input.txt
sed -i '' '/1D0DA/d' my_clean_input.txt
```

After downloading, compiling, and initializing the program, the multiple structural alignment was performed by running mTM-align.

```
cd mTM-align/src
make
./mTM-align/src/mTM-align -i my_clean_input.txt > my_clean_output.txt
```

The alignment section was then extracted from the output file and converted into FASTA format:

```
tail -n +26 my_clean_output.txt | head -n 19 | awk '{print ">"$1;print $2}' > kunitz_clean.ali
```
---
### Build the HMM profile
The resulting multiple structural alignment was subsequently used to build a profile Hidden Markov Model (HMM) with HMMER:
```
hmmbuild ./kunitz_clean.hmm ./kunitz_clean.ali
```
The profile HMM logo was visualized by uploading the generated HMM file to Skylign, which provides a graphical representation of residue conservation and positional variability within the alignment.

---
### VALIDATION
#### TEST SET GENERATION
The integrity of both FASTA datasets was first verified by counting the number of protein entries. The negative dataset contained 574,229 sequences, while the positive dataset contained 398 sequences. Both compressed FASTA files were then extracted to generate the datasets: 

```
zcat -f ./uniprotkb_reviewed_true_NOT_xref_pfam_P_2026_05_06.fasta | grep ">" | wc -l

zcat -f ./uniprotkb_reviewed_true_NOT_xref_pfam_P_2026_05_06.fasta > ./negative_kunitz.fasta
```

```
zcat -f ./uniprotkb_reviewed_true_AND_xref_pfam_P_2026_05_06.fasta | grep ">" | wc -l

zcat -f ./uniprotkb_reviewed_true_AND_xref_pfam_P_2026_05_06.fasta > ./positive_kunitz.fasta
```
---
### BLAST filtering 
To avoid bias during model evaluation, sequences highly similar to those used for HMM training were removed from the positive test set. A BLAST protein database was created from the training seed alignments, and the positive dataset was searched against it using blastp: 
```
makeblastdb -in kunitz_clean.ali -dbtype prot -out training_seeds_db_clean

blastp -query positive_kunitz.fasta -db training_seeds_db_clean -out blast_results_clean.bl8 -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore"
```
Sequences showing ≥95% identity and an E-value ≤ 1e-5 were considered too similar to the training data and were excluded. 

```
awk '($3 >= 95) && ($7 <= 1e-5) {print $1}' \
blast_results_clean.bl8 | cut -d'|' -f2 | sort -u > to_remove_clean.list

grep ">" positive_kunitz.fasta | cut -d'|' -f2 | sort > tot_ids.txt

comm -23 tot_ids_clean.txt <(sort to_remove_clean.list) > IDS_clean.txt
```

The remaining non-redundant positive sequences were collected into a filtered FASTA file for downstream validation.
```
awk -F'|' 'NR==FNR{a[$1];next} /^>/{p=($2 in a)} p' IDS_clean.txt positive_kunitz.fasta > filtered_clean.fasta
```
---
### MODEL TESTING: Optimization and Assessment 
The HMM model was evaluated by running hmmsearch against both the negative and positive datasets. 
The --max option to bypass heuristic filters. To maximize sensitivity and ensure that all sequences were fully evaluated. 
The parameter -Z 1000 sets the effective size of the sequence database to 1000.
This affects the calculation of the E-value, which estimates the number of expected false positives. An E-value equal to 1 indicates that approximately one false positive hit is expected among 1000 searches, corresponding to an estimated false discovery rate (FDR) of about 0.1%.
The negative dataset was used to identify false positives, while the filtered positive dataset was used to evaluate true positive detection.

```
hmmsearch --max --noali --tblout negative_kunitz.search -Z 1000 kunitz.hmm negative_kunitz.fasta
```

```
hmmsearch --max --noali --tblout positive_kunitz_clean.search -Z 1000 kunitz_clean.hmm filtered_clean.fasta
```

To retain only the information of interest from the HMMER output files, the results were parsed to extract three fields: the sequence identifier, the E-value, and the real class label, for both datasets

- Negative dataset:
```
grep -v '^#' negative_kunitz.search | awk '{print $1"\t"$8"\t0"}' > negative_kunitz.match
``` 

- Positive dataset:
    
```
grep -v '^#' positive_kunitz_clean.search | awk '{print $1"\t"$8"\t1"}' > positive_kunitz_clean.match
```
Were appedend _0_ and _1_ as label to distinguish between members belonging to negative and positive members of dataset, respectively.

---
### Reintroducing Missing negative sequences 
To calculate the real performance of the model, it is necessary to also include the negative sequences that produced no match during the HMMER search. 

The identifiers of all negative sequences were extracted from the FASTA file using the following command to creat a list of identifiers:
```
grep ">" negative_kunitz.fasta | awk '{print $1}' | tr -d ">" | sort > negative_kunitz.ids
``` 
Then the identifiers corresponding to matched sequences were extracted:
```
awk '{print $1}' negative_kunitz.match | sort > negative_kunitz_match.ids 
``` 
To identify sequences present in the full dataset but absent from the HMMER matches was used:
```
comm -23 negative_kunitz.ids negative_kunitz_match.ids | less | awk '{print $1"\t100\t0"}' > negative_kunitz.nomatch
``` 
- comm -23 returns only the identifiers present in negative_kunitz.ids but absent from negative_kunitz_match.ids, corresponding to negative sequences with no HMMER hit.

- E-value of 100 was assigned for sequences that produced no match in the HMMER search associated, so they could still be included in the perfomance evaluation.  

The matched and non-matched negatives were then merged into a single file, which contains the complete set of negative sequences (unreviewd) from SwissProt
```
cat negative_kunitz.match negative_kunitz.nomatch > negative_kunitz_tot.match
```
---
# K-Fold Cross-Validation 
To evaluate the robustness of the model, the positive and negative datasets were divided into two subsets for cross-validation.
Before splitting, the negative dataset was randomly shuffled:

```
sort -R negative_kunitz_tot.match
```

The datasets were splitted into two aproximately equal subsets:
- Positives:
    ```
    head -n 191  positive_kunitz.match > kunitz_set_1.txt
    tail -n 191 positive_kunitz.match > kunitz_set_2.txt
    ```
- Negatives:
    ```
    head -n 287115 negative_kunitz_tot.match >> kunitz_set_1.txt
    tail -n 287114 negative_kunitz_tot.match >> kunitz_set_2.txt
    ```

Then was verified that the division was performed correctly:
```
wc kunitz_set_1.txt kunitz_set_2.txt
``` 
287306  861918 7773285 kunitz_set_1.txt
287305  861915 7790552 kunitz_set_2.txt
574611 1723833 15563837 total

---
### Threshold optimization and Performance evaluation 
To evaluate model performance across multiple E-value thresholds (10ˆ-1 - 10ˆ-15), the script performance.py was executed iteratively using a for loop.
- Optimization Set 1
```
for i in `seq 1 15`; do python ./performance.py kunitz_set_1.txt 1e-$i; done > kunitz_set_1.results
``` 


- Optimization Set 2
```
for i in `seq 1 15`; do python ./performance.py kunitz_set_2.txt 1e-$i; done > kunitz_set_2.results
```
Each execution calculated the classification performance metrics at a specific E-value cutoff.


To identify the optimal threshold, the results were sorted according to the Matthews Correlation Coefficient (MCC), which provides a balanced measure of binary classification quality even for highly imbalance datasets. 

```
sort -k6,6gr kunitz_set_2.results
sort -k6,6gr kunitz_set_1.results
```

The threshold corresponding to the highest MCC value was selected as the optimal cutoff for the HMM profile