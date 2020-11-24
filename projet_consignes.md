# Mappeur-SNPcaller de reads sur génome de référence.

<h6>Substitutions autorisées (pas d’insertion ni de délétion).</h6>

<b><u>Objectif :</u></b>	
identifier les SNPs contenus dans les reads et les fournir dans un fichier vcf.

<b><u> Consigne :</u></b>
l’alignement d’un read est testé sur les régions du génome qui partagent
au moins un k-mer avec ce read

Le génome de référence est indexé sous forme d’un FM-index

<b><u> 2 programmes :</u></b><br>
    <ul><i>index.py :</i>
		<ul>Indexer un génome de référence fourni en fasta.</ul>
		<ul>Créer un FM-index : transformée de BW + structure de données 
		additionnelles pour le pattern matching</ul>
		<ul>FMI stocker dans <i>dumped_index.dp</i>.</ul>
		<ul>Utiliser <i>tools_karkkainen_sanders.py</i> pour le tableau de suffixes</ul></ul>
<h6>Prend en entrée les arguments suivants : </h6>

    python index.py --ref [genome_file.fa] --out [dumped_index.dp]
<ul><i>map.py :</i>
	<ul>Mapper un ensemble de reads sur le génome de ref</ul></ul>

<h6>Prend en entrée les arguments suivants : </h6>
	
    python map.py --ref [genome_file.fa] --index [dumped_index.dp] --reads [reads.fa] -k [k_value] --max_hamming [h_value] --min_abundance [m_value] --out snps.vcf

<b><u> Fichier VCF de sortie :</u></b>
<br>Contient 4 informations : 
<ul>
<li>La position du génome muté</li>
<li>L'allèle de référence</li>
<li>L'allèle alternatif</li>
<li>Le nombre de reads mappés avec cet allèle alternatif à cette position</li></ul>
Exepmle de fichier vcf :

    POS REF ALT ABUNDANCE 
    