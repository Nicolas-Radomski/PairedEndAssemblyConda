# Usage
The main Python script PairedEndAssemblyConda.py aims at performing de novo assembly of bacterial genomes from a Python module genomic.py and a Conda environment PairedEndAssembly.
- This workflow run BBnorn (step 1_normalization), Trimmomatic (step 2_trimming), Spades (step 3_assembly), Prokka (4_annotation) and Quast (5_quality), successively.
- The main script PairedEndAssemblyConda.py and module genomic.py (version 20200923, Septembre 2020) were prepared and tested with Python and dependencies below.
- The module genomic.py has to be with the present main script PairedEndAssemblyConda.py to launch it properly.
- The Conda environment PairedEndAssembly has to be prepared as presented below.
- The user can setup his own dependencies in his own bin.
- The paired-end reads must be named ID_R1.fastq.gz and ID_R2.fastq.gz for forward and reverse reads, respectively (ID means sample identifier).
- The IDs have to include a maximum of 16 alphanumeric characters (AZ, az , 09) including potential underscores (_).
- The accents (‘, ¨, ^), space ( ), hyphen (-), and special characters (/, », (, }, =, +, @) are not accepted in the IDs.
- The quality scores paired-end reads must be encoded with Phred33.
# Dependencies
The main script PairedEndAssemblyConda.py and module genomic.py (version 20200923) were prepared and tested with Conda packages below (Name/Version/Build/Channel).
- prokka/1.14.6/pl526_0/bioconda
- python/3.7.8/h6f2ec95_1_cpython/conda-forge
- biopython/1.78/py37h8f50634_0/conda-forge
- bbmap/38.84/h516909a_0/bioconda
- spades/3.13.1/0/bioconda
- quast/5.0.2/py37pl526hb5aa323_2/bioconda
- trimmomatic/0.39/1/bioconda
# Building of the Conda Environment PairedEndAssembly
## 1/ From available targeted Conda packages
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
## 2/ From available updated Conda packages
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
# Launching of the script PairedEndAssemblyConda.py
## 1/ With a single set of paired-end reads
### 1.1/ prepare a single command in a Bash script (bash_PairedEndAssemblyConda.sh)
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
### 1.2/ run the Bash script bash_PairedEndAssemblyConda.sh with sbatch
```
sbatch bash_PairedEndAssemblyConda.sh
```
## 2/ With multiple sets of paired-end reads
### 2.1/ creat a file list_of_IDs.lst including a list of ID samples to process (one ID per line with \n)
```
gedit list_of_IDs.lst
```
### 2.2/ creat a file commands.lst including a list of Bash commands
```
rm commands.lst
for l in `cat list_of_IDs.lst`; do 
	echo "source /global/conda/bin/activate;conda activate PairedEndAssembly;python /global/bio/projets/GAMeR/Nicolas-Radomski/Python/PairedEndAssemblyConda.py -t 48 -n 100 -l 500 -R1 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/$l _R1.fastq.gz -R2 /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/$l _R2.fastq.gz -adap /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/NexteraPE-PE.fa -prot /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109_prot.fasta -nucl /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.fasta -coor /global/bio/projets/GAMeR/Nicolas-Radomski/Python/data/Enteritidis_P125109/Enteritidis_P125109.gff3 -norm /global/conda/envs/PairedEndAssembly/bin/bbnorm.sh -trim /global/conda/envs/PairedEndAssembly/bin/trimmomatic -denovo /global/conda/envs/PairedEndAssembly/bin/spades.py -annot /global/conda/envs/PairedEndAssembly/bin/prokka -qual /global/conda/envs/PairedEndAssembly/bin/quast.py";
done >> commands.lst
sed -i "s@ _R1.fastq.gz@_R1.fastq.gz@" commands.lst
sed -i "s@ _R2.fastq.gz@_R2.fastq.gz@" commands.lst
```
### 2.3/ run the Bash commands of the file commands.lst with sarray
```
sarray -p Research --cpus-per-task=48 -e %x.%N.%j.err -o %x.%N.%j.out --job-name=test-20200925 commands.lst
```
# Illustration
![Workflow](https://github.com/Nicolas-Radomski/PairedEndAssemblyConda/blob/master/illustration.png)
# References
- First version (i.e. SPAdes assembly): Fritsch L., A. Felten, F. Palma, J.F Mariet, N. Radomski, M.Y. Mistou, J.C. Augustin and L. Guillier. Insights from genome-wide approaches to identify variants associated to phenotypes at pan-genome scale: Application to L. monocytogenes ability to grow in cold conditions. 2018, International Journal of Food Microbiology, 291(16): 181-188, doi.org/10.1016/j.ijfoodmicro.2018.11.028
- Second version (i.e. ARTwork): Palma F., Brauge T., Radomski N., Mallet L., Felten A., Mistou M.Y., Brisabois A., Guillier L. and G. Midelet-Bourdin. Dynamics of mobile genetic elements of Listeria monocytogenes persisting in ready-to-eat seafood processing plants in France. 2020, BMC Genomics, 21(1): 130, doi: 10.1186/s12864-020-6544-x
# Acknowledgment
My old colleagues Arnaud Felten and Ludovic Mallet with whom I learned a lot about Python
# Author
Nicolas Radomski
