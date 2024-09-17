#!/bin/bash
#SBATCH --time=11:00:00
#SBATCH --mem=24gb
#SBATCH --partition=gpu
#SBATCH --gpus=p100:1
#SBATCH --mail-type=END,FAIL

module load gcc/8.3.0 intel/18.0.4 cuda/10.1.243 cudnn/7.6.5.32-10.1
. env/bin/activate

python3 deepgoplus.py --go-file go.obo --terms-file data_dir/terms.pkl --train-data-file data_dir/train_data.pkl --test-data-file data_dir/test_data.pkl --model-file data_dir/model.h5 --out-file data_dir/predictions.pkl --logger-file data_dir/training.csv --device gpu:0
