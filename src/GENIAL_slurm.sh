#!/bin/bash
#SBATCH -p Research		   # utilisation de la queue recherche [OBLIGATOIRE]
#SBATCH -o %x.%N.%j.out            # fichier où sera écrit la sortie standart [OBLIGATOIRE]
#SBATCH -e %x.%N.%j.err            # fichier où sera écrit la sortie d'erreur [OBLIGATOIRE]
#SBATCH --cpus-per-task=10	   # nombre de coeur/threads réservés (1 à 40) ici 40 [A MODiFIER]
#SBATCH --job-name=genial_slurm		   # nom du job [A MODIFIER]

. "/global/conda/etc/profile.d/conda.sh"	# sourcer conda [OBLIGATOIRE]

conda activate base		# activer conda [OBLIGATOIRE]

# Pour les lignes suivantes, supprimer le 1er caractère '#' en début de ligne pour utiliser l'outil désiré
#conda activate artwork		# supprimer '#' pour utiliser iVARCall2 (oui oui ce n'est pas une erreur)
#conda activate fastosh		# supprimer '#' pour utiliser FasTosh
#conda activate naura		# supprimer '#' pour utiliser NAuRA
#conda activate genial		# supprimer '#' pour utiliser GENIAL
#conda activate scoary		# supprimer '#' pour utiliser roary/scoary

conda activate genial_db

# Ci-dessous écrire la commande que l'on souhaite lancer
#python GENIAL \
#-f ../../input_slurm_small.tsv \
#-defaultdb vfdb \
#-T 10 \
#-w ../../ \
#-r bacillus_mob_recon \
#-minid 90 \
#-mincov 50 \		
#--mob_recon


#spades.py \
#--plasmid \
#--careful \
#-1 /global/bio/data/GAMeR_DB/BACILLUS/18SBCL215B/18SBCL215B_R1.fastq.gz \
#-2 /global/bio/data/GAMeR_DB/BACILLUS/18SBCL215B/18SBCL215B_R2.fastq.gz \
#-o ../../test_plasmids_spades_18SBCL215B_bac \
#-t 10 

#mob_recon \
#-i ../../assembly_graph.cycs.fasta \
#-o ../../test_mobrecon_spades_18SBCL215B_bac2/ \
#-n 10

python GENIALslurm \
-nas /global/bio/ \
-T 10 \
-w ../../



