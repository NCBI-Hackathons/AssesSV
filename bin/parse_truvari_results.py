from glob import glob
import os
import sys

path_to_files = sys.argv[1]
files = glob(path_to_files + "/giab*.txt")
print("FILES: ", files)
with open ("truvari_summary.txt", "w") as outfile:
    outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("sample", "FN", "FP", "TP", "f1", "precision", "recall"))
    for f in files:
        with open (f, "r") as report:
            lines=report.readlines()
            FN = lines[166].split()[1]
            FP = lines[167].split()[1]
            TP = lines[168].split()[1]
            f1 = lines[176].split()[1]
            precision = lines[177].split()[1]
            recall = lines[178].split()[1]
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(os.path.splitext(os.path.basename(f))[0], FN, FP, TP, f1, precision, recall))
