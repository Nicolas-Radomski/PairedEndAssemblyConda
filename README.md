# Usage
The main Python script PairedEndAssemblyConda.py aims at performing de novo assembly of bacterial genomes from a Python module genomic.py and a Conda environment PairedEndAssembly.
- This workflow run BBnorn (step 1_normalization), Trimmomatic (step 2_trimming), Spades (step 3_assembly), prokka (4_annotation) and Quast (5_quality), successively.
- The main script PairedEndAssemblyConda.py and module genomic.py (version 20200923) were prepared and tested with Python and dependancies below (Septembre 2020).
- The module genomic.py has to be with the present main script PairedEndAssemblyConda.py to lunch it properly.
- The Conda environment PairedEndAssembly has to be prepared as presented below.
- The user can use his own dependancies in his own bin.
- The paired-end reads must be encoded ID_R1.fastq.gz and ID_R2.fastq.gz with ID meaning sample identifier.
# Dependencies
The main script PairedEndAssemblyConda.py and module genomic.py (version 20200923) were prepared and tested with Conda packages below. \
Name                    Version                   Build  Channel \
prokka                    1.14.6                  pl526_0    bioconda \
python                    3.7.8           h6f2ec95_1_cpython    conda-forge \
biopython                 1.78             py37h8f50634_0    conda-forge \
bbmap                     38.84                h516909a_0    bioconda \
spades                    3.13.1                        0    bioconda \
quast                     5.0.2           py37pl526hb5aa323_2    bioconda \
trimmomatic               0.39                          1    bioconda
# Build the Conda Environment PairedEndAssembly
## from available targeted Conda packages
```
conda activate
conda create -n PairedEndAssembly
conda activate PairedEndAssembly
conda search prokka
conda install -c bioconda prokka=1.14.6=pl526_0
conda search python
conda install -c conda-forge python=3.7.8=h6f2ec95_1_cpython
conda search biopython
conda install -c conda-forge biopython=1.78=py37h8f50634_0
conda search bbmap
conda install -c bioconda bbmap=38.84=h516909a_0
conda search spades
conda install -c bioconda spades=3.13.1=0
conda search quast
conda install -c bioconda quast=5.0.2=py37pl526hb5aa323_2
conda search trimmomatic
conda install -c bioconda trimmomatic=0.39=1
```
## from available updated Conda packages
```
conda activate
conda create -n PairedEndAssembly
conda activate PairedEndAssembly
conda install -c bioconda prokka
conda update -c bioconda prokka
conda install -c conda-forge python
conda update -c conda-forge python
conda install -c conda-forge biopython
conda update -c conda-forge biopython
conda install -c bioconda bbmap
conda update -c bioconda bbmap
conda install -c bioconda spades
conda update -c bioconda spades
conda install -c bioconda quast
conda update -c bioconda quast
conda install -c bioconda trimmomatic
conda update -c bioconda trimmomatic
```
# Test the Python script PairedEndAssemblyConda.py with a single set of paired-end reads
## Prepare a single command in a Bash script (sbatch_PairedEndAssemblyConda.sh)
```
#!/bin/bash
#SBATCH -p Research
#SBATCH -o %x.%N.%j.out
#SBATCH -e %x.%N.%j.err
#SBATCH --cpus-per-task=48
#SBATCH --job-name=test-20200922
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
```
## Run the Bash script sbatch_PairedEndAssemblyConda.sh
```
sbatch -p Research --cpus-per-task=48 -e %x.%N.%j.err -o %x.%N.%j.out sbatch_PairedEndAssemblyConda.sh
```
# Run mutiple Bash commands to test the Python script PairedEndAssemblyConda.py with multiple sets of paired-end reads
## remove the file comands.lst
```
rm commands.lst
```
## creat a file list_of_IDs.lst including a list of ID sample to process
```
gedit list_of_IDs.lst
```
## creat a file commands.lst including a list of Bash commands
```
for l in `cat list_of_IDs.lst`; do 
	echo "source /global/conda/bin/activate;conda activate PairedEndAssembly;python /global/bio/projets/GAMeR/Nicolas-Radomski/Python/PairedEndAssemblyConda.py -t 48 -n 100 -l 500 -R1 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/$l _R1.fastq.gz -R2 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/$l _R2.fastq.gz -adap /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/NexteraPE-PE.fa -prot /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109_prot.fasta -nucl /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.fasta -coor /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.gff3 -norm /global/conda/envs/PairedEndAssembly/bin/bbnorm.sh -trim /global/conda/envs/PairedEndAssembly/bin/trimmomatic -denovo /global/conda/envs/PairedEndAssembly/bin/spades.py -annot /global/conda/envs/PairedEndAssembly/bin/prokka -qual /global/conda/envs/PairedEndAssembly/bin/quast.py";
done >> commands.lst
```
## remove space before _R1.fastq.gz and _R2.fastq.gz in the file commands.lst
```
sed -i "s@ _R1.fastq.gz@_R1.fastq.gz@" commands.lst
sed -i "s@ _R2.fastq.gz@_R2.fastq.gz@" commands.lst
```
## lunch Bash commands of the file commands.lst with sarray
```
sarray -p Research --cpus-per-task=48 -e %x.%N.%j.err -o %x.%N.%j.out --job-name=test-20200925 commands.lst
```
# Author
Nicolas Radomski