# Author
Nicolas Radomski

# Usage
The Python script PairedEndAssemblyConda.py aims at performing de novo assembly of bacterial genomes from a Python module genomic.py and a Conda environment PairedEndAssembly.
This workflow run BBnorn (step 1_normalization), Trimmomatic (step 2_trimming), Spades (step 3_assembly), prokka (4_annotation) and Quast (5_quality), successively.
The main script PairedEndAssemblyConda.py and module genomic.py (version 20200923) were prepared and tested with Python 3.6.2 (Septembre 2020).
The module genomic.py has to be with the present main script PairedEndAssemblyConda.py to lunch it.
The Conda environment PairedEndAssembly has to be prepared as presented below.

# Build the Condan Environment PairedEndAssembly
conda activate \
which conda \
conda create -n PairedEndAssembly \
conda activate PairedEndAssembly \
conda install -c conda-forge -c bioconda prokka \
conda install python=3.6.2=h02fb82a_12 \
conda install -c conda-forge biopython \
conda update -c conda-forge biopython \
conda install -c bioconda bbmap \
conda install -c bioconda spades \
conda install -c bioconda quast \
conda install -c bioconda trimmomatic

# Run a single the bash command to test the Python script PairedEndAssemblyConda.py with a single paired-end reads
#!/bin/bash \
#SBATCH -p Research \
#SBATCH -o %x.%N.%j.out \
#SBATCH -e %x.%N.%j.err \
#SBATCH --cpus-per-task=48 \
#SBATCH --job-name=test-20200922 \
source /global/conda/bin/activate;conda activate PairedEndAssembly; \
python /global/bio/projets/GAMeR/Nicolas-Radomski/Python/PairedEndAssemblyConda.py \
 -t 48 \
 -n 100 \
 -l 500 \
 -R1 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/ERR3997398_R1.fastq.gz \
 -R2 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/ERR3997398_R2.fastq.gz \
 -adap /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/NexteraPE-PE.fa \
 -prot /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109_prot.fasta \
 -nucl /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.fasta \
 -coor /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.gff3 \
 -norm /global/conda/envs/PairedEndAssembly/bin/bbnorm.sh \
 -trim /global/conda/envs/PairedEndAssembly/bin/trimmomatic \
 -denovo /global/conda/envs/PairedEndAssembly/bin/spades.py \
 -annot /global/conda/envs/PairedEndAssembly/bin/prokka \
 -qual /global/conda/envs/PairedEndAssembly/bin/quast.py

# Run mutiple Bash commands to test the Python script PairedEndAssemblyConda.py with multiple paired-end reads
## remove the file comands.lst
rm commands.lst
## creat a file list_of_IDs.lst including a list of ID sample to process
## creat a file commands including a list of Bash commands
for l in `cat list_of_IDs.lst`; do \
	echo "source /global/conda/bin/activate;conda activate PairedEndAssembly;python /global/bio/projets/GAMeR/Nicolas-Radomski/Python/PairedEndAssemblyConda.py -t 48 -n 100 -l 500 -R1 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/$l _R1.fastq.gz -R2 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/$l _R2.fastq.gz -adap /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/NexteraPE-PE.fa -prot /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109_prot.fasta -nucl /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.fasta -coor /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.gff3 -norm /global/conda/envs/PairedEndAssembly/bin/bbnorm.sh -trim /global/conda/envs/PairedEndAssembly/bin/trimmomatic -denovo /global/conda/envs/PairedEndAssembly/bin/spades.py -annot /global/conda/envs/PairedEndAssembly/bin/prokka -qual /global/conda/envs/PairedEndAssembly/bin/quast.py";
done >> commands.lst
## remove space before _R1.fastq.gz and _R2.fastq.gz in the file commands.lst
sed -i "s@ _R1.fastq.gz@_R1.fastq.gz@" commands.lst \
sed -i "s@ _R2.fastq.gz@_R2.fastq.gz@" commands.lst \
## lunch Bash commands of the file commands.lst with sarray
sarray -p Research --cpus-per-task=48 -e %x.%N.%j.err -o %x.%N.%j.out --job-name=test-20200925 commands.lst