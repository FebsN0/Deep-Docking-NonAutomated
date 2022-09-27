import argparse
import builtins as __builtin__
import glob
import gzip
import os
from contextlib import closing
from multiprocessing import Pool


parser = argparse.ArgumentParser()
parser.add_argument('-if', '--is_final', required=True, help='True/False for is this the final iteration?')
parser.add_argument('-n_it', '--iteration_no', required=True, help='Number of current iteration')
parser.add_argument('-protein', '--protein', required=True, help='Name of the DD project')
parser.add_argument('-file_path', '--file_path', required=True, help='Path to the project directory, excluding project directory name')
parser.add_argument('-t_pos', '--tot_process', required=True, help='Number of CPUs to use for multiprocessing')
parser.add_argument('-score', '--score_keyword', required=True, help='Score keyword. Name of the field storing the docking score in the SDF files of docking results')
parser.add_argument('-zincid', '--zincid_keyword', required=True, help='zincid keyword. Name of the field storing the zincid in the SDF files of docking results. This is not so obvious that zincid is at the first line of each molecule section in the sdf file')
parser.add_argument('-site', '--siteSelected', required=True, help='select the binding site used in this cluster [1,2,3,4]: ')
parser.add_argument('-set', '--datasetSelected', required=True, help='choose the dataset to process [train,test,valid]: ')

io_args = parser.parse_args()
is_final = io_args.is_final
n_it = int(io_args.iteration_no)
protein = io_args.protein
file_path = io_args.file_path
tot_process = int(io_args.tot_process)
key_word_score = str(io_args.score_keyword)
key_word_zincid = str(io_args.zincid_keyword)
siteSel = str(io_args.siteSelected)
setSel = io_args.datasetSelected

if is_final == 'False' or is_final == 'false':
    is_final = False
elif is_final == 'True' or is_final == 'true':
    is_final = True
else:
    raise TypeError('-if parameter must be a boolean (true/false)')
# mol_key = 'ZINC'
print("Keyword Score: ", key_word_score)
print("Keyword zincid: ", key_word_zincid)
print("Selected site: ", siteSel)
print("Selected dataset: ", setSel)

def get_scores(ref):
    scores = []
    for line in ref:  # Looping through the molecules
# it takes the first line.. but in many cases ZINCID is not in the first line
#        zinc_id = line.rstrip()
        line = ref.readline()
        # '$$$' signifies end of molecule info
        while line != '' and line[:4] != '$$$$':  # Looping through its information and saving scores

            tmp = line.rstrip().split('<')[-1]
#modified part: search for ZINCID
            if key_word_zincid == tmp[:-1]:
                zinc_id = ref.readline().rstrip()
            if key_word_score == tmp[:-1]:
                tmpp = float(ref.readline().rstrip())
                if tmpp > 100 or tmpp < -100:
                    print(zinc_id, tmpp)
                else:
                    scores.append([zinc_id, tmpp])

            line = ref.readline()
    return scores


def extract_glide_score(filen):
    scores = []
    try:
        # Opening the GNU compressed file
        with gzip.open(filen, 'rt') as ref:
            scores = get_scores(ref)

    except Exception as e:
        print('Exception occured in Extract_labels.py: ', e)
        # file is already decompressed
        with open(filen, 'r') as ref:
            scores = get_scores(ref)

    if 'test' in setSel:
        new_name = 'testing'
    elif 'valid' in setSel:
        new_name = 'validation'
    elif 'train' in setSel:
        new_name = 'training'
    else:
        print("Could not generate further training set")
        exit()

    with open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/' + new_name + '_' + 'labels.txt', 'a') as ref:
        for z_id, gc in scores:
            ref.write(str(gc) + ',' + z_id + '\n')


if __name__ == '__main__':
    files = []
    iter_path = file_path + '/' + protein + '/iteration_' + str(n_it)

    # Checking to see if the labels have already been extracted:
    sets = ["training", "testing", "validation"]
    files_labels = glob.glob(iter_path + "/*_labels.txt")
    foundAll = True
    for s in sets:
        found = False
        print(s)
        for f in files_labels:
            set_name = f.split('/')[-1].split("_labels.txt")[0]
            if set_name == s:
                found = True
                print('Found')
                break
        if not found:
            foundAll = False
            print('Not Found')
            break
    if foundAll:
        print('Labels have already been extracted...')
        print('Remove "_labels.text" files in \"' + iter_path + '\" to re-extract')
        exit(0)

    # Checking to see if this is the final iteration to use the right folder
    if is_final:
        path = file_path + '/' + protein + '/after_iteration/docking/site_'+ siteSel + '/'+ setSel + '/block*/*.sdf*'
    else:
    #check all sdf files in every block subdirectory in the set selected
    #in this way, I will avoid to generate further data in a different docking directory
        path = iter_path + '/docking/site_'+ siteSel + '/'+ setSel + '/block*/*.sdf*'
        path_labels = iter_path + '/*labels*'

    for f in glob.glob(path):
        files.append(f)

    print("num files in", path, ":", len(files))
    print("Files:", [os.path.basename(f) for f in files])
    if len(files) == 0:
        print('NO FILES IN: ', path)
        print('CANCEL JOB...')
        exit(1)
    if 'test' in setSel:
        new_name = 'testing'
    elif 'valid' in setSel:
        new_name = 'validation'
    elif 'train' in setSel:
        new_name = 'training'
    else:
        print("Could not generate new training set")
        exit()

    #write new file
    with open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/' + new_name + '_' + 'labels.txt', 'w') as ref:
        ref.write('r_i_docking_score' + ',' + 'ZINC_ID' + '\n')
    # Parallel running of the extract_glide_score() with each file path of the files array
    with closing(Pool(len(files))) as pool:
        pool.map(extract_glide_score, files)