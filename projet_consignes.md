<h1>Mappeur-SNPcaller de reads sur génome de référence. </h1>

<h6>Implémenter un mappeur-SNPcaller de reads sur génome de référence. Les substitutions sont autorisées (pas d’insertion ni de délétion).</h6>

<h3><b><u>Objectif :</u></b></h3>
Identifier les SNPs contenus dans les reads et les fournir dans un fichier vcf.

<h3><b><u> Consigne :</u></b></h3>
L’alignement d’un read est testé sur les régions du génome qui partagent
au moins un <i>k-mer</i> avec ce read

Le génome de référence est indexé sous forme d’un FM-index
#
<h3><b><u> 2 programmes :</u></b><br></h3>
    <ul><b><i>index.py :</i></b>
    <ul>
		<li>Indexer un génome de référence fourni en fasta.</li>
		<li>Créer un FM-index : transformée de BW + structure de données 
		additionnelles pour le pattern matching</li>
		<li>FMI stocker dans <i>dumped_index.dp</i>.</li>
		<li>Utiliser <i>tools_karkkainen_sanders.py</i> pour le tableau de suffixes</li></ul></ul>
<h6>Prend en entrée les arguments suivants : </h6>

    python index.py --ref [genome_file.fa] --out [dumped_index.dp]
<ul><b><i>map.py :</i></b>
	<ul><li>Mapper un ensemble de reads (<i>reads.fa</i> sur un génome de référence (<i>genome_file.fa)</i></li></ul></ul>

<h6>Prend en entrée les arguments suivants : </h6>
	
    python map.py --ref [genome_file.fa] --index [dumped_index.dp] --reads [reads.fa] -k [k_value] --max_hamming [h_value] --min_abundance [m_value] --out snps.vcf

Pour chaque read : 
<lu>
<li>Les potentielles positions sont détectées vie l'ancrage de k-mers</li>
<li>Attention à ne pas tester plusieurs fois la même position de mapping</li>
<li>Si la mapping échoue (trop de substitution), rien n'est stocké</li>
<li>Si plusieurs positions possibles :
<br>On garde la position montrant moins de substitutions
<br>Si plusieurs dans ce cas : position la plus à gauche du génome</li>

#
<h3><b><u> Fichier VCF de sortie :</u></b></h3>
<br>On a 4 informations à notre disposition : 
<ul>
<li>La position du génome muté</li>
<li>L'allèle de référence</li>
<li>L'allèle alternatif</li>
<li>Le nombre de reads mappés avec cet allèle alternatif à cette position</li></ul>
<h6>Exepmle :</h6>

    POS REF ALT ABUNDANCE 
    
Indiquer seulement les variants dont l'abondance est supérieure à <i>m_value</i> définie par l'option <i>-min_abundance</i>.
<br>Les lignes seront ordonnées par positions sur le génome. La position est <i>0-based</i>, c'est-à-dire que la première lettre du génome est considérée comme la position 0.

<br>Le fichier débutera avec les lignes suivantes:

 
    REF: nom de fichier contentant le génome de référence
    READS: nom de fichier contentant les reads
    K: Valeur de k utilisée pour effectuer l'ancrage
    MAX_SUBST: la variable h_value
    MIN_ABUNDANCE: la variable m_value