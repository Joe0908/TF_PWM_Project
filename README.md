Pipeline for TF PWM Annotation

Overview:
This pipeline processes a FASTA file containing human transcription factor (TF) sequences and annotates them with Position Weight Matrix (PWM) information using similarity searches and external tools. The goal is to assign the most likely PWM to each input TF using multiple methods: sequence identity, homology, ModCRE, and AlphaFold.
 
1. Input
Input file: A FASTA file of TF sequences.

2. MMseqs2 Similarity Search
MMseqs2 is used to compare the input sequences against a merged reference database of JASPAR and CisBP sequences.
Script used: run_mmseqs.sh
This will:
Identify identical sequences → Column 2: Identical_PWM
Identify sequences with >70% similarity → Column 3: Homologous_PWM
Identify sequences with 50–70% similarity → Column 4: Relatively_Homologous_PWM

3. Chart Generation
Script used: generate_tf_pwm_chart.py

4. Add TF Family Information
Using the file TF_accession_family.csv, which maps UniProt IDs to TF families, we fill in the first column: TF_family.
Script used: update_chart_with_TF_family.py
Output: TF_PWM_chart_2.tsv


5. Run ModCRE on Remaining Sequences
For sequences without assigned PWM from identity/homology matches, we run the ModCRE model.
The PWM name assigned by ModCRE is recorded in Column 5: ModCRE.
Script used: modcre.py
→ Output: TF_PWM_chart_3.tsv


6. Run AlphaFold on Remaining Sequences
For sequences not resolved by ModCRE, AlphaFold is used. Predicted models are matched to PWMs.
The PWM name is recorded in Column 6: AlphaFold.
Script used: alphafold.py
→ Output: TF_PWM_chart_4.tsv

7. Filter Invalid Results
Two text files (wrong_model_modcre.txt and wrong_model_af3.txt) contain sequences deemed non-TF by ModCRE or AlphaFold. These sequences should be removed.
Script used: final_table.py
Final output: TF_PWM_chart_final.tsv
