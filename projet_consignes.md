<h1>Mappeur-SNPcaller de reads sur génome de référence. </h1>

<h6>Substitutions autorisées (pas d’insertion ni de délétion).</h6>

<h3><b><u>Objectif :</u></b></h3>
identifier les SNPs contenus dans les reads et les fournir dans un fichier vcf.

<h3><b><u> Consigne :</u></b></h3>
l’alignement d’un read est testé sur les régions du génome qui partagent
au moins un k-mer avec ce read

Le génome de référence est indexé sous forme d’un FM-index
#
<h3><b><u> 2 programmes :</u></b><br></h3>
    <ul><b><i>index.py :</b></i>
    <ul>
		<li>Indexer un génome de référence fourni en fasta.</li>
		<li>Créer un FM-index : transformée de BW + structure de données 
		additionnelles pour le pattern matching</li>
		<li>FMI stocker dans <i>dumped_index.dp</i>.</li>
		<li>Utiliser <i>tools_karkkainen_sanders.py</i> pour le tableau de suffixes</li></ul></ul>
<h6>Prend en entrée les arguments suivants : </h6>

    python index.py --ref [genome_file.fa] --out [dumped_index.dp]
<ul><b><i>map.py :</i></b>
	<ul><li>Mapper un ensemble de reads sur le génome de ref</li></ul></ul>

<h6>Prend en entrée les arguments suivants : </h6>
	
    python map.py --ref [genome_file.fa] --index [dumped_index.dp] --reads [reads.fa] -k [k_value] --max_hamming [h_value] --min_abundance [m_value] --out snps.vcf

#
<h3><b><u> Fichier VCF de sortie :</u></b></h3>
<br>Contient 4 informations : 
<ul>
<li>La position du génome muté</li>
<li>L'allèle de référence</li>
<li>L'allèle alternatif</li>
<li>Le nombre de reads mappés avec cet allèle alternatif à cette position</li></ul>
<h6>Exepmle de fichier vcf :</h6>

    POS REF ALT ABUNDANCE 
    