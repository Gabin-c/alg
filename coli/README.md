# Données test *E. coli* 



Ce répertoire contient 

- `ecoli_sample.fasta`: un fichier fasta contenant les 150 000 premières nucléotides de *E. coli*. 

- `ecoli_mutated_reads_1000.fasta`: un fichier de reads au format fasta, où chaque read est sur une ligne.

  - Proviennent d'une version mutée de `ecoli_sample.fasta`, distante de 1000 substitutions, 2 insertions et 2 délétions. 
  - Couverture 20x
  - taux d'erreurs de séquençage 1%

-  `ecoli_mutated_truth.vcf`; un fichier simplifié de variants, au format identique (hormis les commentaires) au format demandé pour le projet (notamment 0-based). 

- `validation.py`: un script permettant de calculer la qualité des résultats d'un fichier produit par votre outil. 

  - Exemple: 

    ```bash
    python validation.py ecoli_mutated_truth.vcf my_predictions.vcf
    Nb Variants 1004 in truth
    Nb TP 535 in truth & in pred
    Nb FP 0 not in truth & in pred
    Nb FN 469 in truth & not in pred
    Precision 100.0 %
    Recall 53.5 %
    ```

    
