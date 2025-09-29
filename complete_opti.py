
#Script by Looplye

#################################################################################################################
#                                                                                                               #
#       Needed: .pileup / .mpileup                                                                              #
#                                                                                                               #
#       Remember: put all the pileup files in one folder,                                                       #
#           the results folder will be created in this specific folder                                          #
#                                                                                                               #
#################################################################################################################

# ========== Parameters ==========

File = ""
Linux_file = ""
mac_file = ""
#Pathway of the file where all the .pileup/.mpileup to be analysed must be + where the file 'results' will be created


#researche_change
rate1 = 65 #minimum percentage of presence in all reads for it to be detected as a mutation
nbr_reads_min = 10 #minimum depth for it to analyse the localisation 

#regroupe_mut
sliding_window = 400 #size before analysis (will regroup all detected mutations in a nucelotides window of this number)
interval = 10 #will move of this number of nucleotides forward before doing the next analysis
margin = 800 #number of nucleotides not analysed in the beginning and at the end

#calcul_score
rate2 = 0 #minimum percentage of presence in all samples for the mutation to be compared to others

#main
rate3 = 0.05 #minimum score of differenciation for a window to be selected
rate4 = 0.8 #minimum score of differenciation for a window to be selected in the 2nd randomisation
SAVE = 1000 #will save a temporary file every SAVE (if the program stop abrutely, you'll have a backup and won't have to redo everything)
#random part (main & random)
nbr_tests = 1000000 #number of random selection of windows the program will test
nbr_windows = 5 #numer of windows selected each time



import platform, re, csv, tempfile, os, glob
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from tqdm import tqdm
from datetime import datetime
import pandas as pd
import numpy as np


# ========== Fonctions ==========

#Finds and counts the indels, finds the most recurring one, cleans the raw pileup file
def researche_change(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as f_in, \
         open(output_file, 'w', encoding='utf-8') as f_out:

        #defines and cleans the elements in the file
        for line in f_in:
            elements = line.strip().split('\t')[:-1]
            mutation = elements[4]
            nbr_reads = int(elements[3])
            
            #detects then counts the indels
            if ('+' in mutation or '-' in mutation) and nbr_reads >= nbr_reads_min:
                nbr_mut = mutation.count('+') + mutation.count('-')  
                ratio1 = (nbr_mut / nbr_reads) *100
                    
                #selection of the indels only if enough
                if ratio1 >= rate1:
                    elements.insert(5, str(ratio1))

                    #Finds the most recurring indel and replace all the indels by it
                    compteur = Counter(re.split(r'[.,]', mutation.upper()))
                    if compteur:
                        maj, count = compteur.most_common(1)[0]
                        elements[4] = maj

                        #Calculates the percentage of presence of the most recurring indel
                        ratio2 = count / nbr_reads * 100
                        elements.insert(6, str(ratio2))
                        
                        #Writes a file with the new info (2 ratios and most recurring indel)
                        f_out.write('\t'.join(elements) + '\n')

#Regroups all the indels found for one window together, for every windows
def regroup_mut(input_file, output_file):
    data_chr = defaultdict(dict)
    #set up the dictionnary

    #defines the elements in the file
    with open(input_file, 'r', encoding='utf-8') as f_in:
        for line in f_in:
            elements = line.strip().split('\t')
            chr_name = elements[0]
            position = int(elements[1])
            mutation = elements[4]

            data_chr[chr_name][position] = mutation
            
    #defines the beginning and the end of the chromosome that will be analysed
    with open(output_file, 'w', encoding='utf-8') as f_out:   
        for chr_name in sorted(data_chr):
            pos_dict = data_chr[chr_name]
            beginning = margin
            max_pos = max(pos_dict) - margin        
        
            #defines the choosen window
            while beginning <= max_pos:
                end = beginning + sliding_window -1
                total_mut = 0
                mut_detect = False

                #detects the indels & their position                
                for i in range(beginning, end + 1):
                    if i in pos_dict:
                        mut = pos_dict[i]
                        indel = re.findall(r'([\+\-])(\d+)', mut)   

                        if indel:
                            mut_detect = True

                        #regroups all the indels in one number
                        for symbol, number_str in indel :
                            number = int(number_str)
                            if symbol == '+':
                                total_mut += number
                            else:
                                total_mut -= number

                #writes the new file
                if mut_detect:
                    final_symbol = '+' if total_mut > 0 else ''
                    f_out.write(f"{chr_name}\t{beginning}-{end}\t{final_symbol}{total_mut}\n")

                #starts again for a new window
                beginning += interval

#Calls researche_change & regroup_mut for every file
def process_sample(chemin_fichier, tmpdir):
    #regroups all files
    nom = os.path.splitext(os.path.basename(chemin_fichier))[0]
    f_out_1 = os.path.join(tmpdir, f"{nom}.pileup")
    f_out_2 = os.path.join(tmpdir, f"{nom}_pre-score.pileup")

    #calls both functions
    researche_change(chemin_fichier, f_out_1)
    regroup_mut(f_out_1, f_out_2)

    return f_out_2, nom

#Creates a dictionnary for every windows in every sample & writes all the windows that contain indels
def find_window(input_file, output_file):
    #creates the dictionary
    data = {}
    for file, name in input_file:
        data[name] = {}
        
        #defines the elements in the file
        with open(file, 'r', encoding='utf-8') as f_in:
            for line in f_in:
                elements = line.strip().split('\t')
                chr_name = elements[0]
                window = elements[1]
                mut_str = elements[2]
                try:
                    mut_val = int(mut_str)
                except ValueError:
                    mut_val = mut_str
                
                #updates the dico
                data[name][(chr_name, window)] = mut_val

    #initials with a loop
    resultats = {}
    all_windows = set()
    for data_sample in data.values():
        all_windows.update(data_sample.keys())
    for window in all_windows:
        diff_count = 0
        total = 0

        #counts the number of variation window by window for each sample
        for name in data:
            val = data[name].get(window, 0)
            if val is None:
                continue
            total += 1
            if val != 0:
                diff_count += 1

        ratio = (diff_count / total)*100 if total > 0 else 0

        #stocks only windows with indels
        if diff_count > 0:  
            resultats[window] = (diff_count, total, ratio)

    #writes them down
    with open(output_file, 'w', encoding='utf-8') as f_out:
        for (chr_name, window), (diff_count, total, ratio) in sorted(resultats.items()):
            f_out.write(f"{chr_name}\t{window}\t{ratio:.2f}\t{diff_count} / {total}\n")

    return data

#Create score file
def calc_score(input_file, output_file, dic_mut):
    score_window = {}
    
    with open(input_file, 'r', encoding='utf-8') as f_in:
        lines = f_in.readlines()

        #defines elements
        for line in tqdm(lines, unit="window"):
            elements = line.strip().split('\t')
            chr_name = elements[0]
            window = elements[1]
            try:
                ratio = float(elements[2])
            except ValueError:
                continue
            
            #selection
            if ratio > rate2:
                mutations = {}

                #retrieve mutation from dictionnary for each window
                for name, wind_dic in dic_mut.items():
                    for (chrom, wind), mut in wind_dic.items():
                        if chrom == chr_name and wind == window:
                            mutations[name] = mut

                #format them
                val_valides = []
                for val in mutations.values():
                    if isinstance(val, int):
                        val_valides.append(val)
                    elif isinstance(val, str) and re.search(r'[\+\-]\d+', val):
                        val_valides.append(0)

                #score calculation
                nbr_val = len(set(val_valides))
                nbr_total = len(dic_mut)

                if nbr_total > 0:
                    score = nbr_val / nbr_total 
                else:
                    score = 0
                score_window[(chr_name, window)] = score

        #writes results
        with open(output_file, 'w', encoding='utf-8') as f_out:
            name_samp = list(dic_mut.keys())

            f_out.write("chromosome\tWindow\tScore")
            for name in name_samp:
                f_out.write(f"\t{name}")
            f_out.write("\n")

            #process potential errors / when no mutation was detected for some samples
            for (chr_name, window), score in sorted(score_window.items()):
                line_mut = []
                for name in name_samp:
                    if name not in dic_mut:
                        l_mut = f"ERROR : no sample named {name}"
                        print(f"ERROR : no sample named {name}")
                    elif (chr_name, window) not in dic_mut[name]:
                        l_mut = 0
                    else:
                        l_mut = dic_mut[name][(chr_name, window)]
                    line_mut.append(l_mut)

                line = f"{chr_name}\t{window}\t{score:.2f}\t" + "\t".join(str(m) for m in line_mut) + "\n"
                f_out.write(line)

#Creates first csv file "analysis"
def visual(fichier_entree, fichier_csv):
    #open both pileup and csv files
    with open(fichier_entree, 'r', encoding='utf-8') as f_in, \
         open(fichier_csv, 'w', newline='', encoding='utf-8') as f_out:

        #copy lines by lines
        csv_writer = csv.writer(f_out, delimiter=';')
        for ligne in f_in:
            elements = ligne.strip().split('\t')
            csv_writer.writerow(elements)

#Classify the samples with indels found
def group_numpy(call_matrix):

    #transcribes every column in tuple 
    keys = [tuple(call_matrix[:, i]) for i in range(call_matrix.shape[1])]

    #creates a dictionary with column number and their tuple 
    d = defaultdict(list)
    for col_id, indel in enumerate(keys):
        d[indel].append(col_id)
    #each d value => sample group with the same profil

    return list(d.values())

#Select random windows and send them to group_numpy to classify them
def random(windows, nbr_windows, calls, sample_names, seed):
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(windows), size=nbr_windows, replace=False)
    #random & size determination

    #sub group of windows chosen & their info (calls)
    sub_calls   = calls[idx]
    sub_windows = windows[idx]

    #classify samples & score calculation
    groups = group_numpy(sub_calls)     
    score  = len(groups) / sub_calls.shape[1]
    
    groups_named = [[sample_names[i] for i in g] for g in groups]

    return {
        'score': score,
        'windows': sub_windows.tolist(),
        'groups': groups_named
    }

#Display time when change step
def log(message):
    heure = datetime.now().strftime("%H:%M")
    print(f"[{heure}] {message}")

# ========== Main program ==========

#Main function
def main():
    #Uses the right pathway folder for the platform
    if platform.system() == "Linux":
        pileup_file = Linux_file
    elif platform.system() == "Darwin":
        pileup_file = mac_file
    else:
        pileup_file = File

    #Regroups all the pileup file and create the results folder
    all_pileup = glob.glob(os.path.join(pileup_file, "*.pileup")) + \
                 glob.glob(os.path.join(pileup_file, "*.mpileup"))

    results_folder = os.path.join(pileup_file, "results")
    os.makedirs(results_folder, exist_ok=True)

    #Creates a temporary folder where temporary files will be stored until used then deleted
    with tempfile.TemporaryDirectory() as tempdir:
        temp_file = []
        log("Analysis of the (m)pileup files...")

        #uses multi processor for a faster analysis (function process_sample)
        with ProcessPoolExecutor() as ex:
            #tqdm for the loading bar
            for res in tqdm(ex.map(partial(process_sample, tmpdir=tempdir), all_pileup),
                    total=len(all_pileup)):
                temp_file.append(res)

        #Centralising the files (functions find_window & calcul_score)
        log("Calculating the score...")
        comparison_file = os.path.join(tempdir, "comparison_temp.pileup")
    
        mut = find_window(temp_file, comparison_file)
    
        analysis_file = os.path.join(results_folder, "windows_analysis.pileup")
        calc_score(comparison_file, analysis_file, mut)

        #Creation of the first csv file "windows analysis"
        log("Creation of the analysis file...")
        analysis_file_csv = os.path.join(results_folder, "windows analysis.csv")
        visual(analysis_file, analysis_file_csv)

    #reads analysis file 
    log("Reading analysis file...")
    df = pd.read_csv(
        analysis_file,
        sep='\t',
        encoding='utf-8-sig',
        low_memory=False           
    )
    #defines the columns
    df.columns = df.columns.str.upper().str.strip()

    #keeps only the windows > rate3
    df_ok = df[df['SCORE'] > rate3].reset_index(drop=True)
    
    #defines the name
    sample_names = df_ok.columns[3:].tolist()

    #creates a samples data table in string
    calls = df_ok.iloc[:, 3:].to_numpy(dtype='U32') 
    #creates a chrom & name data table also in string   
    windows = (
            df_ok['CHROMOSOME'].astype(str) + ':' + df_ok['WINDOW'].astype(str)
        ).to_numpy(dtype='U32')

    log("Randomisation...")

    #generates random numbers
    RNG  = np.random.default_rng()
    #max taille
    seed = RNG.integers(0, 2**32, size=nbr_tests)
    
    #initialises
    results = []
    
    #multiprocessing
    with ProcessPoolExecutor() as ex:
        #partial prefill run_one
        func = partial(random, windows, nbr_windows, calls, sample_names)
        temp_file = os.path.join(results_folder, "score_temp.csv")

        #recover data if it exists
        if os.path.exists(temp_file):
            log("data recovery...")
            results = pd.read_csv(temp_file, sep=';').to_dict(orient="records")
            start_index = len(results)
            log(f"starting again from {start_index}\n   {nbr_tests - start_index} remaining")
        else:
            start_index = 0

        seeds_restants = seed[start_index:]

        #loading bar + for each r, a differend seed 
        for i, r in enumerate(tqdm(ex.map(func, seeds_restants), total=len(seeds_restants))):
            results.append(r)

            index_global = start_index + i + 1
            if index_global % SAVE == 0:
                pd.DataFrame(results).to_csv(temp_file, sep=';', index=False)
                log(f"temporary save : {i+1}")

    #writes final file "score.csv"
    df_res = pd.DataFrame(results)
    exit_file = os.path.join(results_folder, "score.csv")

    df_res.to_csv(exit_file, sep=';', index=False)           
    
    #deletes temp_file if/when error
    if os.path.exists(temp_file):
        os.remove(temp_file)

    log(f"Starting second randomisation with only windows from selections with score > {rate4}...")

    # Lecture du fichier score.csv
    df_score = pd.read_csv(exit_file, sep=';')

    # Extraction des fenêtres présentes dans les sélections où le score > rate4
    fenetres_filtrees = set()
    for _, row in df_score.iterrows():
        try:
            if float(row["score"]) > rate4:
                # row["windows"] est une liste sous forme de chaîne → on la convertit
                win_list = eval(row["windows"])  # attention : on part du principe que c'est sûr
                fenetres_filtrees.update(win_list)
        except Exception:
            continue

    log(f"{len(fenetres_filtrees)} unique windows kept after filtering.")

    # Conversion en tableau numpy
    mask_windows = np.isin(windows, list(fenetres_filtrees))
    windows_filtered = windows[mask_windows]
    calls_filtered = calls[mask_windows]

    if len(windows_filtered) < nbr_windows:
        log("Warning: not enough windows after filtering to run the second randomisation.")
    else:
        # Génération des seeds pour la deuxième randomisation
        seed2 = RNG.integers(0, 2**32, size=nbr_tests)
        results2 = []

        with ProcessPoolExecutor() as ex:
            func2 = partial(random, windows_filtered, nbr_windows, calls_filtered, sample_names)
            temp_file_2 = os.path.join(results_folder, "score_temp_filtered.csv")

            for i, r in enumerate(tqdm(ex.map(func2, seed2), total=len(seed2))):
                results2.append(r)

                index_global = start_index + i + 1
                if index_global % SAVE == 0:
                    pd.DataFrame(results).to_csv(temp_file_2, sep=';', index=False)
                    log(f"temporary save : {i+1}")

        # Sauvegarde du second fichier score
        exit_file2 = os.path.join(results_folder, f"score_filtered_{rate4}.csv")
        pd.DataFrame(results2).to_csv(exit_file2, sep=';', index=False)
        log(f"Second score file saved: {exit_file2}")

    log("Finish !")


#to not crash if using windows to run
if __name__ == "__main__":
    log("Starting the process...")
    from multiprocessing import freeze_support
    freeze_support()  
    main()


