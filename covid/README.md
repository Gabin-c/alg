## Projet ALG 2020 -- Mapping et SNP calling de données de séquençage COVID



## Ce dossier contient :

* `SRX9435498_subset10000.fasta` : un fichier de 10000 reads d'un vrai séquençage d'un échantillon d'un patient atteint de la covid-19, id = SRX9435498.

* `MN908947.fasta` : le génome de référence du virus sars-cov2, taille = 29903 bp.

  

## D'où viennent ces données ?

* Données de séquençage [SRX9435498](https://www.ncbi.nlm.nih.gov/sra/?term=SRX9500342) : 
  * Origine : fichier soumis par le  New Mexico Department of Health Scientific Laboratory le 4 novembre 2020
  * Echantillon : prélèvement sur un patient atteint du COVID au Nouveau mexique (USA), le 24/08/2020.
  * Technologie : séquençage Illumina sur une librairie tiled-amplicon (amplification spécifique du génome viral sars-cov2)
  * Lien download : https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?exp=SRX9435498&cmd=search&m=downloads&s=seq
  * Download fasta : 163 Mo, 556,596 paires de reads de tailles variables (max 150 pb)
  * Sous-échantillonnage de 10000 reads de taille au moins 100 bp sans autre caractère que A,C,G, et T en évitant les 50000 premiers. Fichier `SRX9435498_subset10000.fasta`

* Génome de référence : sars-cov2, première version du génome de référence (janv 2020)

  * Publication : Wu et al , Feb 2020 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7094943/
* ID : MN908947 :https://www.ncbi.nlm.nih.gov/nuccore/MN908947
  * Fichier `MN908947.fasta`



Note : autres sites pour télécharger des données de séquençage COVID-19 :

* https://www.ncbi.nlm.nih.gov/sars-cov-2/ (puis cliquer sur "View in SRA")

* https://www.covid19dataportal.org/sequences?db=sra-experiment-covid19#search-content

  Sur ce site de l'EBI, on peut filtrer par l'origine géographique. Il y a quelques jeux de données français, mais ils étaient un peu trop gros pour cette expérience (plusieurs Go chacun).

