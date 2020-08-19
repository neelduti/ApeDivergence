#Run in terminal
ProFl=$(pwd)
################### 1. data download
mkdir $ProFl/data $ProFl/data/cds $ProFl/data/gff $ProFl/doc
cd $ProFl/data/cds
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-98/fasta/*/cds/*.cds.all.fa.gz .
cd $ProFl/data/gff
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-98/gff3/*/*.98.gff3.gz .
cd ..
gzip -d ./*/Homos_sapiens*
gzip -d ./*/Pan_troglodytes*
gzip -d ./*/Gorilla_gorilla*
gzip -d ./*/Nomascus_leucogenys*
#delete all the zipped files once unzipping is complete
rm ./*/*.gz

################### 2. Primary transcript selection 
#if possible run these as four separate jobs on a HPC
python3 $ProFl/code/S2_PrimaryTranscriptSelection.py -s Homos_sapiens
python3 $ProFl/code/S2_PrimaryTranscriptSelection.py -s Pan_troglodytes
python3 $ProFl/code/S2_PrimaryTranscriptSelection.py -s Gorilla_gorilla
python3 $ProFl/code/S2_PrimaryTranscriptSelection.py -s Nomascus_leucogenys

################### 3. Pairwise all vs all BLASTp
python3 $ProFl/code/S3_PairwiseBlastP.py -s Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys



################### 4 Pairwise synteny
mkdir $ProFl/data/synteny $ProFl/doc/synteny
#if possible run this job on a HPC
python3 $ProFl/code/S4_Synteny.py -s Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys


################### 5 MSA
mkdir $ProFl/data/msa $ProFl/doc/msa
python3 $ProFl/code/S5_MSA.py -s Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Nomascus_leucogenys -f Homo_sapiens


