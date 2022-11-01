import argparse
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument('-protein', '--protein', required=True, help='Name of DD project')
parser.add_argument('-file_path', '--file_path', required=True, help='PAth to project folder without project folder name')
parser.add_argument('-n_it', '--n_it', required=True, help='Number of current iteration')
parser.add_argument('-mdd', '--morgan_directory', required=True, help='PAth to Morgan fingerprint directory')

# adding parameter for where to save all the data to:
parser.add_argument('-save', '--save_path', required=False, default=None)

io_args = parser.parse_args()

protein = io_args.protein
n_it = int(io_args.n_it)
mdd = io_args.morgan_directory

#filepath.
dp=io_args.file_path

DATA_PATH = dp+'/'+protein  # Now == file_path/protein
SAVE_PATH = io_args.save_path
# if no save path is provided we just save it in the same location as the data
if SAVE_PATH is None: SAVE_PATH = DATA_PATH

add = mdd

try:
    os.mkdir(SAVE_PATH + '/iteration_' + str(n_it) + '/simple_job_predictions')
except OSError:
    pass

for f in glob.glob(SAVE_PATH + '/iteration_' + str(n_it) + '/simple_job_predictions/*'):
    os.remove(f)

# TODO Add to GUI create a project
ct = 0.9
time = '0-10:30'

# temp = []
part_files = []

for i, f in enumerate(glob.glob(add + '/*.txt')):
    part_files.append(f)

ct = 1
for f in part_files:
    with open(SAVE_PATH + '/iteration_' + str(n_it) + '/simple_job_predictions/simple_job_' + str(ct) + '.sh',
              'w') as ref:
        ref.write('#!/bin/bash\n')
        ref.write('#SBATCH --account=def-jtus\n#SBATCH --nodes=1\n')
        ref.write('#SBATCH --gres=gpu:a100:1     #####  AWARE: MODIFY NAME GPU BASED ON THE USED CLUSTER #####\n')
        ref.write('#SBATCH --cpus-per-task=5\n#SBATCH --mem-per-cpu=2G\n#SBATCH --time 05:00:00\n#SBATCH --ntasks=2\n')
        ref.write('#SBATCH --mail-user=fab.alt@protonmail.com\n#SBATCH --mail-type=FAIL\n')
        ref.write('source ~/.bashrc\n')
        ref.write('module load cuda cudnn\nsource ~/envDeepDock/bin/activate\n')
        cwd = os.getcwd()
        ref.write('cd {}\n'.format(cwd))
        ref.write('python -u $DEEPDOCKNA/phase_2-3/Prediction_morgan_1024.py' + ' ' + '-fn' + ' ' + f.split('/')[
            -1] + ' ' + '-protein' + ' ' + protein + ' ' + '-it' + ' ' + str(n_it) + ' ' + '-mdd' + ' ' + str(
            mdd) + ' ' + '-file_path' + ' ' + dp + '\n')
        ref.write("\n echo complete")

    ct += 1
