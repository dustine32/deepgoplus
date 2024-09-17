## Train model on GO release from GOWithHierarchy files
For the purpose of training a model with a GO release of annotations mapped to PANTHER long IDs. The standard procedure for DeepGOPlus is to start with a SwissProt file in TXT format. The protein and annotation data in this file is parsed into a pandas dataframe and exported into a pickle to be loaded as training data. I decided to insert our data into the DeepGOPlus process at this point by converting our usual GOWithHierarchy annotation files into SwissProt TXT format to be used as input into `uni2pandas.py`.
Step 5 requires the `diamond` executable, which can be downloaded from https://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz (or latest release).
### Step 1: Convert GOWithHierarchy files to SwissProt TXT format
Inputs
* `GOWithHierarchy-ALL-14.1.dat` - Concatenated GOWithHierarchy files for `ALL` three GO aspects. At the time generated we used PANTHER14.1 gene IDs, so an ID forward mapping file to PANTHER15.0 was also supplied.
* `all.fasta` - FASTA sequences for all proteins in the PANTHER15.0 database for inclusion in the TXT output.
* `organism.dat` - For obtaining taxon ID from the PANTHER ID's species code (e.g., `HUMAN`, `CAEEL`, `SCHPO`).
* `--id_forward_mapping` - librarySeqMap used to (if supplied) update the annotated protein long ID to its forwarded 15.0 ID
```bash
python3 go_hierarchy2uniprot.py GOWithHierarchy-ALL-14.1.dat PANTHER15.0/library_building/all.fasta organism.dat --id_forward_mapping PANTHER15.0/library_building/nodeForwardTracking/librarySeqMap_14.1to15 > go_hierarchy_20191021_uniprot.txt
```
### Step 2: Convert SwissProt TXT format to Pandas DataFrame pkl file
First need to gzip the `_uniprot.txt` file.
```bash
gzip -c go_hierarchy_20191021_uniprot.txt > go_hierarchy_20191021_uniprot.txt.gz
python3 uni2pandas.py --swissprot-file go_hierarchy_20191021_uniprot.txt.gz --go-file go.obo --out-file go_hierarchy_20191021_uniprot.pkl
```
### Step 3: Split the data into training and testing sets
```bash
python3 deepgoplus_data.py --go-file go.obo --data-file go_hierarchy_20191021_uniprot.pkl --out-terms-file data_dir/terms.pkl --train-data-file data_dir/train_data.pkl --test-data-file data_dir/test_data.pkl
```
### Step 4: Train the model
Requires a GPU.
```bash
python3 deepgoplus.py --go-file go.obo --terms-file data_dir/terms.pkl --train-data-file data_dir/train_data.pkl --test-data-file data_dir/test_data.pkl --model-file data_dir/model.h5 --out-file data_dir/predictions.pkl --logger-file data_dir/training.csv --device gpu:0
```
Running this with slurm on HPC, I had to request a GPU node on discovery with the following headers:
```bash
#SBATCH --partition=gpu
#SBATCH --gpus=p100:1
```
Also need to ensure the correct version of CUDA and cuDNN are loaded:
```slurm
module load gcc/8.3.0 intel/18.0.4 cuda/10.1.243 cudnn/7.6.5.32-10.1
```
An example slurm script is provided in `deepgoplus_train_gpu.slurm`.
### Step 5: Run diamond
Extract FASTA file for the sequences in our UniProt TXT file:
```bash
python3 diamond_data.py --data-file go_hierarchy_20191021_uniprot.pkl --out-file go_hierarchy_20191021_uniprot.fasta
```
Make diamond DB from protein FASTA
```bash
./diamond makedb --in go_hierarchy_20191021_uniprot.fasta --db data_dir/train_data.dmnd
```
Run diamond on PANTHER human sequences (homo_sapiens.fasta, the sequences we want to use DeepGOPlus on)
```bash
./diamond blastp -d data_dir/train_data.dmnd --more-sensitive -q homo_sapiens.fasta --outfmt 6 qseqid sseqid bitscore > data_dir/test_diamond.res
```
### Step 6: Predict GO annotations
Compress files as required by predict.py
```bash
gzip -c data_dir/test_diamond.res > data_dir/test_diamond.res.gz
gzip -c homo_sapiens.fasta > data_dir/homo_sapiens_15.fasta.gz
```
Predict annotations for PANTHER human sequences using trained model and diamond results.
```bash
python3 predict.py --in-file data_dir/homo_sapiens_15.fasta.gz --out-file results_20191021.tsv.gz --go-file data_dir/go.obo --model-file data_dir/model.h5 --terms-file data_dir/terms.pkl --annotations-file go_hierarchy_20191021_uniprot.pkl --diamond-file data_dir/test_diamond.res.gz
```