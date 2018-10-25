#!/usr/bin/env python3

"""
look at the argparse section
"""

import argparse
import os
import re
import glob
import sys
import json
import csv
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from lxml.html import parse
from urllib.request import urlopen
import matplotlib
import pandas as pd
import subprocess as sbp
matplotlib.use("Agg")  # Force matplotlib not to use any X-windows backend.
from matplotlib import pyplot as plt

# Settings part: make CDD and configure

def which(program):
    """
    :param program: any program you wish to find
    :return: a path to executable that are in the PATH variable
    """
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def make_config():
    # make a dict to write to a json file (config.json)
    file_units = ["domain"]
    soft_units = ["RepeatMasker", "RepeatModeler", "ORFfinder", "rpstblastn", "BuildDatabase", "bedtools", "makeprofiledb"]

    config_dict = dict.fromkeys(file_units + soft_units)

    # configuring the software paths
    for unit in soft_units:
        found_path = which(unit)
        if found_path is None:
            print("%s not found on your computer. Enter a valid path to an executable file" % unit)
            given_path = Path(input("enter the path > "))
            while not given_path.is_file():
                print("this path is invalid: file not found!")
                given_path = Path(input("enter the path > "))
            config_dict[unit] = str(given_path)
            print("Ok...")
        else:
            print("%s found here - %s" % (unit, found_path))
            config_dict[unit] = found_path
            print("OK...")

    # configure the paths to files which will be used in the pipeline
    for unit in file_units:

        # separate action for repeat type
        # if unit == "RepeatType":
        #     config_dict[unit] = input("What repeats you're interested in? [LINE, Gypsy, Penelope]... > ")
        #     print("OK...")

        # separate action for Conserved Domain Database
        # NB: this settings are valid only for RT domains
        if unit == "domain":
            domain_name = Path(input("Enter a name for Conserved Domain Database %s... > " % unit))
            config_dict[unit] = str(domain_name)
            print("OK...")

        # separate action for prefix (species name) to added to all other files
        elif unit == "prefix":
            prefix = input("Enter a prefix for your files (species name) > ")
            config_dict[unit] = prefix
            print("OK...")

        # and anything else
        else:
            unit_path = Path(input("Enter a path to %s... > " % unit))
            if unit_path.is_file() or str(unit_path) == "skip":
                config_dict[unit] = str(unit_path)
                print("OK...")
            else:
                print("ERROR. No such file or your path is wrong!\nQuit")
                sys.exit()

    # add fixed file names and paths to the config
    # config_dict["RepeatMasker Output"] = config_dict["genome"][:-3] + ".out"
    config_dict["RepeatModeler Output"] = "consensi.fa.classified"
    config_dict["GetORF Output"] = "hitdata.xml"
    config_dict["ORFinder Input"] = "matches_w_hits.fa"
    config_dict["CENSOR Output"] = "Unknown_classified.fa"
    # if you do not run makeCDD script, you have to add domain name
    config_dict["CDD"] = os.getcwd()+"/LocalCDD/" + config_dict["domain"]

    config_str = json.dumps(config_dict)
    with open("config.json", 'w') as outfile:
        outfile.write(config_str)
    print("OK. All files are configured.")


def make_local_cdd(dir_for_cd="LocalCDD"):
    """
    :param dir_for_cd: a name for a directory to store local CDD
    :return: no
    """
    domain_name = read_config("domain")
    # find smp files and make list of them, call it "smp_files_list.pn"
    workdir = os.getcwd()
    files = glob.glob("*.smp")

    # check if there is no such files
    if len(files) == 0:
        print("No smp files found in %s" % workdir)
        if os.path.isdir(dir_for_cd):
            ans = input("LocalCDD directory found. Maybe a database already exists? [y/n] > ")
            if ans == "y":
                print("Ok...")
                # add to config
                cdd_path_full = os.getcwd() + "/" + dir_for_cd + "/" + domain_name
                # os.chdir("../")
                add_config("CDD", cdd_path_full)
                print("Local Conserved Domain Database is made and full path to it is added to the config.")
                sys.exit(0)
            else:
                print("ERROR! Put *smp files in current working directory.")
                sys.exit(1)
        if not os.path.isdir(dir_for_cd):
            print("ERROR! Put *smp files in current working directory.")
            sys.exit(1)

    # check if dir for CDD exists and make it, if not exists, and move into
    if not os.path.isdir(dir_for_cd):
        os.mkdir(dir_for_cd)

    # check if smp table exists and move it
    cd_table = glob.glob("*.csv")
    if len(cd_table) == 0:
        print("no smp table found...")
        #print("assume there is only one CD type...")
    else:
        print("smp table found and added to the config...")
        new_path = workdir + "/" + dir_for_cd +"/" + cd_table[0]
        shutil.move(cd_table[0], new_path)
        add_config("cd_table", new_path)

    # move smp files to new dir
    os.system("mv *.smp %s" % dir_for_cd)
    os.chdir(dir_for_cd)

    with open("smp_files_list.pn", "w") as out:
        for f in files:
            out.write("%s\n" % f)
    print("A list of smp files has been compiled.")

    print("Database name - %s" % domain_name)

    # make local database
    cmd = "makeprofiledb -in smp_files_list.pn -out %s " % domain_name
    os.system(cmd)

    cdd_path_full = os.getcwd() + "/" + domain_name

    # go to level up to write config in right directory
    os.chdir("../")
    add_config("CDD", cdd_path_full)
    print("Local Conserved Domain Database is made and full path to it is added to the config.")

# GetCons part

def check_types(seq_file=""):
    """
    collects types of sequences in a file (if the file exists)
    :param seq_file: fasta file with sequences
    :return: a sorted list of repeat types in the file
    """

    if seq_file == "":
        seq_file = read_config("RepeatModeler Output")
    try:
        with open(seq_file) as f:
            repeats = [re.findall(r'#(.*)', rec.id)[0] for rec in SeqIO.parse(f, 'fasta')]
            repeats.sort()
            for rep in set(repeats):
                print("%s" % rep)
    except FileNotFoundError:
        print("File %s not found!" % seq_file)
        sys.exit()


def read_config(unit, home=False):
    """
    reads config file and returns a value of a requested key
    :param unit: key from config.json
    :param home: where is config.json: in 'home' dir or in CWD?
    :return: value of this key
    """
    # if home:
    #     config_location = os.getcwd() + "/config.json"
    # else:
    #     config_location = "config.json"
    config_location = "config.json"
    try:
        conf = open(config_location).readline()
        jconf = json.loads(conf)
    except FileNotFoundError:
        try:
            # check if config exists one level above (RepeatModeler and MakeCDD case)
            conf = open("../"+config_location).readline()
            jconf = json.loads(conf)
        except FileNotFoundError:
            try:
                conf = open("../../"+config_location).readline()
                jconf = json.loads(conf)
            except FileNotFoundError:
                print("ERROR! config.json not found!")
                sys.exit()
    return jconf[unit]


def add_config(unit, val, home=False):
    """
    add new item to the config file
    :param unit: new key to add
    :param val: a value for the key
    :param home: where is config.json: in 'home' dir or in CWD?
    :return: there is no return
    """

    # if home:
    #     config_location = os.getcwd() + "/config.json"
    # else:
    #     config_location = "config.json"
    config_location = "config.json"
    try:
        conf = open(config_location).readline()
        jconf = json.loads(conf)
        config_path = config_location
    except FileNotFoundError:
        try:
            # check if config exists one level above (RepeatModeler and MakeCDD case)
            conf = open("../"+config_location).readline()
            jconf = json.loads(conf)
            config_path = "../"+config_location
        except FileNotFoundError:
            try:
                conf = open("../../"+config_location).readline()
                jconf = json.loads(conf)
                config_path = "../../" + config_location
            except FileNotFoundError:
                print("ERROR! config.json not found!")
                print("config location - %s" % config_location)
                sys.exit()
    jconf[unit] = val
    conf = json.dumps(jconf)
    with open(config_path, 'w') as outfile:
        outfile.write(conf)
    # print("New item was added to the config - %s" % unit)


def repeat_collector(filename, word, merge=False):
    """
    repeat_collector extracts entries containing particular word in their headers (e.g. repeat type)
    from any fasta file (e.g. RepeatModeler output file - consensi.fa.classified)
    :param filename: name of the fasta file with sequences
    :param word: a word in sequence header you are interested in
    :param merge: if this run is for subsequent merge, default=False
    """
    try:
        f = open(filename)
    except FileNotFoundError:
        print("Error! File %s not found!" % filename)
        sys.exit()

    extracted_sequences = [record for record in SeqIO.parse(f, 'fasta') if word in record.id]

    if len(extracted_sequences) == 0:
        print("No %s sequences found in %s" % (word, filename))
        # rename initial file into final one
        if merge:
            cmd = "mv " + word + "_consensi.fa.classified" + " All_" + word + "_consensi.fa"
            os.system(cmd)
            print("Your repeat consensi file was renamed")
            # sys.exit()
        else:
            sys.exit()

    fileout = word + "_" + filename.split("/")[-1]
    with open(fileout, 'w') as output:
        SeqIO.write(extracted_sequences, output, 'fasta')
    if word == "Unknown":
        add_config("Unknown repeats file", fileout)
    print('%s repeats from RepeatModeler output have been collected' % word)


def censor_parser(url):
    """
    censor_parser reads and parses CENSOR html output
    and renames unclassified consensuses from RepeatModeler
    arg: unknown sequences file
    """

    # file = input("unknown sequences file? > ")
    file = read_config("Unknown repeats file")
    # insert URL from GIRI Censor output
    # 'http://www.girinst.org/cgi-bin/censor/show_results.cgi?id=79886&lib=root'
    # url = input("CENSOR output URL > ")

    parsed = parse(urlopen(url))
    doc = parsed.getroot()
    tables = doc.findall('.//table')
    text_tables = [tables[i].text_content().split(' ') for i in range(len(tables))]
    # collect consensus name and rep class/family it belongs to in tuples
    cons_type = [(item[12], item[-5]) for item in text_tables if len(item) >= 15]  # >= 3 is minimal

    # make new sequence file
    classified = []
    for rec in SeqIO.parse(file, 'fasta'):
        for item in cons_type:
            if rec.id == item[0]:
                newid = rec.id[:-7] + item[1]  # cut off the word 'Unknown'
                newrec = SeqRecord(rec.seq, id=newid, description='')
                classified.append(newrec)
    with open("Unknown_classified.fa", 'w') as f:
        SeqIO.write(classified, f, 'fasta')

    print('Repeats from web CENSOR output have been collected.')


def merger(file_a, file_b, fill):
    """
    merge file_a (repeats from repmod output) and file_b (repeats from censor output)
    :param file_a: repeats from collector_repmod
    :param file_b: repeats from collector_censor
    :param fill: a word to insert in out file name (typically - repeat name)
    :return: the name of new concatenated file
    """
    out = "All_" + fill + "_consensi.fa"
    command = "cat " + file_a + " " + file_b + " > " + out
    os.system(command)
    return out


def get_cons(mge_type='', censor=False, url='', unknown=False, recollect=False, merge=False, b_types=False, n_types=False, a_types=False, standard=False):
    """

    :param censor: use CENSOR or not
    :param url: URL address of CENSOR output
    :param unknown: collect unknown sequences
    :param recollect: collect your repeats from newly classified sequences
    :param merge:
    :param b_types:
    :param n_types:
    :param a_types:
    :param standard:
    :param mge_type: type of MGE to search. For some actions is unnecessary
    :return:
    """
    if standard:
        filenm = read_config("RepeatModeler Output")
        # collect repeats of interest
        repeat_collector(filenm, mge_type)
        # collect Unknown repeats
        repeat_collector(filenm, "Unknown")
        # print("Now you may run CENSOR -> http://www.girinst.org/censor/")
    elif censor:
        censor_parser(url)
        print("Censor output has been parsed and your repeats have been collected")
    elif unknown:
        filenm = read_config("RepeatModeler Output")
        reptype = "Unknown"
        repeat_collector(filenm, reptype)
        print("%s repeats have been collected" % reptype)
    elif recollect:
        # reptype = read_config("RepeatType")
        fname = read_config("CENSOR Output")
        repeat_collector(fname, mge_type, merge=True)
        print("collecting %s repeats from %s..." % (mge_type, fname))
        file1 = mge_type + "_" + read_config("RepeatModeler Output")
        file2 = mge_type + "_" + read_config("CENSOR Output")
        merger(file1, file2, fill=mge_type)
        print("%s repeats from newly classified sequences have been collected and merged with previous data set" % mge_type)
    elif merge:
        # file1 = input("file #1 to merge > ")
        # file2 = input("file #2 to merge > ")
        file1 = mge_type + "_" + read_config("RepeatModeler Output")
        file2 = mge_type + "_" + read_config("CENSOR Output")
        out = merger(file1, file2, fill=read_config("RepeatType"))
        print("Two data sets have been merged - see %s" % out)
    elif b_types:
        file = read_config("RepeatModeler Output")
        types = check_types(file)
        print("The following types of repeats found in %s:" % file)
        for type in sorted(set(types)):
            print(type, types.count(type), sep=" - ")
    elif n_types:
        file = read_config("CENSOR Output")
        types = check_types(file)
        print("The following types of repeats found in %s:" % file)
        for type in sorted(set(types)):
            print(type, types.count(type), sep=" - ")
    elif a_types:
        file1 = read_config("RepeatModeler Output")
        file2 = read_config("CENSOR Output")
        concat = merger(file1, file2, fill="Repeats")
        types = check_types(concat)
        print("The following types of repeats found in %s:" % concat)
        for type in sorted(set(types)):
            print(type, types.count(type), sep=" - ")
    else:
        print("Nothing to do :(")

# GetSeq part


def merge_overlap(df, ov=0):
    """
    to merge overlapping matches in one data frame record,
    runs inside ori_to_csv() function
    :param df: a DataFrame object with all the same values in 'strand' column
    :param ov: the overlap value; matches closer to each other than this 'ov'
    will be considered as overlapping
    :return: a pandas data frame object with merged matches
    """
    # see c_strand.csv as an example
    # df = pd.read_csv("c_strand.csv")
    start = df["start"]
    skewed_start = start[1:].append(pd.Series(0))
    diff = pd.Series(skewed_start.values - start.values, index=start.index)
    # df["diff"] = diff
    df["overlap"] = df["length"] > diff - ov
    # what blocks of overlapping matches we have?
    # '1' - means it's the 1st block; '2' - 2nd and so on;
    # 'False' in 'overlap' column means end of the block
    blocks = []
    n = 1
    for item in df['overlap']:
        if item:
            blocks.append(n)
        else:
            blocks.append(n)
            n += 1

    df['blocks'] = blocks
    # and now merge rows with overlapping matches:
    clean_groups = []
    for name, group in df.groupby('blocks'):
        # find the 'stop' of the last row in this table:
        z = int(group.tail(1)['stop'])
        # assign new value to the 'stop' of the 1st row of this table:
        group.set_value(min(group.index), 'stop', z)
        # remove all other rows in this table:
        # but first select these rows:
        these_rows = list(group.index[1:len(group) + 1])
        group = group.drop(these_rows)
        clean_groups.append(group)

    clean_df = pd.concat(clean_groups)
    # we don't need 'blocks' and 'overlaps' anymore:
    del clean_df["blocks"]
    del clean_df["overlap"]

    # update 'length' column:
    clean_df["length"] = clean_df["stop"] - clean_df["start"] + 1

    return clean_df


def ori_to_csv(rm_file, len_cutoff=0, expansion=0, make_plot=True, bins=200):
    """
    to read ori file, compute lengths of matches, drop off matches with
    lengths less than is set by 'len_value' parameter, expand coordinates, merge overlapping matches.
    Also prompts to print a histogram of lengths
    :param rm_file: RepeatMasker *ori.out file (text)
    :param len_cutoff: all matches with shorter length will be removed (integer)
    :param expansion: bps to expand matches (integer)
    :param make_plot: whether to make a density plot (logical)
    :param bins: number of bins in histogram
    :return: cleaner and shorter table of csv format (writes it to disk)
    """

    print("RepeatMasker table cleaning and transformation...")
    rm_tab = pd.read_csv(rm_file, delim_whitespace=True, header=None)

    # remove THESE columns
    rm_tab = rm_tab.drop([0, 1, 2, 3, 7, 9, 10, 11, 12, 13], axis=1)

    # headers = ["contig", "start", "stop", "strand"]
    length = rm_tab[6] - rm_tab[5] + 1
    rm_tab[9] = length

    # assign column names
    rm_tab.columns = ['contig', 'start', 'stop', 'strand', 'length']

    # drop off too short matches
    tab_sort = rm_tab[rm_tab['length'] >= len_cutoff]
    print("matches shorter than %i bp have been dropped" % len_cutoff)

    # expanding coordinates
    tab_sort["start"] = tab_sort["start"] - expansion
    tab_sort["stop"] = tab_sort["stop"] + expansion

    # split DF by contigs
    cleaned_df_lst = []

    # number of matches being processed,
    # frequency of message printing,
    # number of contigs
    x = 0
    num_con = len(tab_sort['contig'].unique())
    step = 500

    for name1, contig_group in tab_sort.groupby('contig'):
        x += 1
        # split unique contig DF by strand
        for name2, strand_group in contig_group.groupby('strand'):
            # get cleaned version of strand-specific DF
            cleaned_result = merge_overlap(strand_group)
            cleaned_df_lst.append(cleaned_result)
        if x % step == 0:
            print("Number of contigs processed: %i/%i" % (x, num_con))

    tab_clean = pd.concat(cleaned_df_lst)
    tab_clean.to_csv('file.fa.ori.out.cleaned.csv', index=False)

    # plot/histogram printing
    # ans = input("print a histogram? [y/n] >  ")
    if make_plot:
        plot3 = tab_clean['length'].hist(bins=bins)
        # plot = tab_clean['length'].plot(kind='kde')
        fig3 = plot3.get_figure()
        fig3.savefig(rm_file + "_hist.png")
        plt.close(fig3)


def cut_match_dict(genomic_seq, rm_ori_csv, begin=1, end=2, output="excised_matches.fa"):
    """
    to cut the matches from genomic sequences according to coordinates in RepMask output file
    counts matches per scaffold and rev-complements them if it is necessary
    :param genomic_seq: *fa/*fna file with whole genome assembly or with scaffolds where matches have been found
    :param rm_ori_csv: a *csv table from ori_to_csv()
    :param begin: start of a match
    :param end: end of a match
    :param output: output file name (default "excised_matches.fa")
    :return: writes fasta file with sequences of matches [extended & counted]
    """
    print("excision...")
    genome_dict = dict()

    # make genomic dictionary scaffold_name:SeqRecord
    for rec in SeqIO.parse(genomic_seq, "fasta"):
        genome_dict[rec.id] = rec

    # open cleaned rm table
    match_table = csv.reader(open(rm_ori_csv))
    s = int(begin)
    e = int(end)
    match_list = []

    # this is match counter
    count = 1

    # and name of 'previous' scaffold
    previous = ''

    # this skips the first row (header) of the CSV file
    next(match_table)

    #####################################################
    # main block of code:
    # here it counts matches in a scaffold, excises them
    # and rev-complements if it is necessary
    #####################################################
    for row in match_table:
        try:
            scaffold = genome_dict[row[0]]
        except KeyError:
            print("Error: scaffold not found: %s" % str(row[0]))
        else:

            # here it counts matches per one scaffold
            if row[0] == previous:
                count += 1
            else:
                count = 1
            previous = row[0]

            # and here it makes match's SeqRecord
            # define [and adjust] start and end coordinates
            start = int(row[s])
            if start <= 0:
                start = 1
            end = int(row[e])
            if end > len(scaffold.seq):
                end = len(scaffold.seq)

            # get sequence of a match
            match_seq = scaffold.seq[start - 1:end]

            # rev comp it if "C"
            if row[3] == "C":
                match_seq == match_seq.reverse_complement()
                descr = "-_%i_%i-%i" % (count, start, end)
            else:
                descr = "+_%i_%i-%i" % (count, start, end)

            accession = scaffold.id
            record = SeqRecord(match_seq, id=accession, description=descr)
            match_list.append(record)

    # write to a file
    with open(output, 'w') as handle:
        SeqIO.write(match_list, handle, 'fasta')
    genome_dict.clear()
    print("RM-matches excised...")
    add_config("GetSeq Output", output)


def translate_match(matches_file):
    """
    This script makes six translations and writes them into a file
    Default translation table = 1 (standard)
    :param matches_file: file with counted matches from previous function
    :return: file with translations for RPS-Blast
    """

    print("translating...")
    trans_list = list()
    for seq_record in SeqIO.parse(matches_file, 'fasta'):
        accession = seq_record.description
        rc_sequence = seq_record.seq.reverse_complement()

        # translating by all frames
        for frame in range(3):
            # for leading strand
            orf = seq_record.seq[frame:]
            protein = orf.translate()
            record = SeqRecord(protein, id=accession, description='+%i' % frame)
            trans_list.append(record)

            # for reverse complement strand
            rc_orf = rc_sequence[frame:]
            rc_protein = rc_orf.translate()
            rc_record = SeqRecord(rc_protein, id=accession, description='-%i' % frame)
            trans_list.append(rc_record)

    # write to file
    handle = open('translations.faa', 'w')
    SeqIO.write(trans_list, handle, 'fasta')


def bed_tools(rm_file, ref_genome, ori=False, expansion=0, output_prefix="excised_matches", bins=200):
    """
    invokes bedtools suite to excise repeat matches
    :param rm_file: RepeatMasker output file (*ori.out or *out)
    :param ref_genome: genome assembly file
    :param ori: if the RepeatMasker output is *ori.out or not (default: True)
    :param expansion: merge nearby (within x bp) repetitive elements into a single entry.
    :param output_prefix: prefix for the output file name (default "excised_matches")
    :param bins: number of bins in a histogram
    :return: there is no return :)
    """

    # every time add MGE's name to any output
    mge_name = read_config("RepeatType") + "_"
    output_prefix = mge_name + output_prefix
    if ori:
        rm_file = rm_file[:-4] + ".ori.out"
        print("Using *ori.out file with bedtools")
        ori2bed = "sed 's/^\s*//' %s | sed -E 's/\s+/\t/g' | cut -f5-9 | " \
                  "awk 'BEGIN { OFS=\"\t\" } { print $1, $2-1, $3, $4, $5}' | " \
                  "sed -E 's/\tC/\t-/g' | sed -E 's/\t\(/\tElement\t\(/g' | " \
                  "sort -k1,1 -k2,2n > temp.bed; mergeBed -s -d %s -i temp.bed > temp2.bed; mergeBed -i temp2.bed > " % (rm_file, str(expansion)) + "%s; rm temp.bed; rm temp2.bed" % (rm_file + ".bed")
        os.system(ori2bed)
        # sbp.call(["sh", "ori_to_bed.sh", rm_file, str(expansion)])
        sbp.call(["bedtools", "getfasta", "-s", "-fi", ref_genome, "-bed", rm_file + ".bed", "-fo", output_prefix + str(expansion) + ".fa"])
    else:
        if rm_file[-3:] == "bed":
            print("bed-file detected")
            sbp.call(["bedtools", "getfasta", "-s", "-fi", ref_genome, "-bed", rm_file + ".bed", "-fo", output_prefix + str(expansion) + ".fa"])
        else:
            print("Using RepeatMasker output with headers")
            out2bed = "tail -n +4 %s | sed 's/^\s*//' | sed -E 's/\s+/\t/g' | cut -f5-9 | " \
                      "awk 'BEGIN { OFS=\"\t\" } { print $1, $2-1, $3, $4, $5}' | sed -E 's/\tC/\t-/g' | " \
                      "sed -E 's/\t\(/\tElement\t\(/g' | sort -k1,1 -k2,2n > " \
                      "temp.bed; mergeBed -s -d %s -i temp.bed > temp2.bed; mergeBed -i temp2.bed > " % (rm_file, str(expansion)) + "%s; rm temp.bed; rm temp2.bed" % (rm_file + ".bed")

            os.system(out2bed)

            # sbp.call(["sh", "out_to_bed.sh", rm_file, str(expansion)])
            sbp.call(["bedtools", "getfasta", "-s", "-fi", ref_genome, "-bed", rm_file + ".bed", "-fo", output_prefix + str(expansion) + ".fa"])
    # add output to the config
    add_config("GetSeq Output", output_prefix + str(expansion) + ".fa")
    # make plot
    df = pd.read_table(rm_file + ".bed", header=None)
    length = df[2] - df[1] + 1
    plot2 = length.hist(bins=bins)
    fig2 = plot2.get_figure()
    fig2.savefig(output_prefix + str(expansion) + ".png")
    plt.close(fig2)
    # make some descriptive stats on lengths
    des_stat = length.describe()
    # write statistics to a file
    des_stat.to_csv(output_prefix + str(expansion) + ".stats", header=False, index=True, sep="\t")
    # print on the screen
    print("descriptive statistics about length of excised matches:")
    print(pd.DataFrame(des_stat))


def get_seq(pandas=False, merge=0, ori=False, rm_tab=""):
    """
    this function cuts repeat sequences from a genome according to coords from repeat masker table.
    it uses bedtools or pandas (the latter could be very slow).
    :param pandas: logical, to use pandas or not. Default False.
    :param merge: "merge all repeats within M bp into a single entry"
    :param ori: "explicitly set that RepeatMasker output is *ori.out"
    :param rm_tab: repeat masker table to use
    :return: no
    """
    if rm_tab == "":
        repmask_out = read_config("RepeatMasker Output")
    else:
        repmask_out = rm_tab
    genome = read_config("genome")
    # remove '.gz' suffix
    genome = genome[:-3]
    # d = os.getcwd()
    if not os.path.isfile(repmask_out):
        print("ERROR! File %s not found!" % repmask_out)
        sys.exit()
    elif not os.path.isfile(genome):
        print("ERROR! File %s not found!" % genome)
        sys.exit()

    if pandas:
        print("using pandas library")
        prfx = read_config("prefix")
        ori_to_csv(repmask_out, expansion=merge)
        cut_match_dict(genome, prfx + ".fa.cleaned.csv")
    elif not ori:
        bed_tools(repmask_out, genome, expansion=merge, ori=False)
    elif ori:
        bed_tools(repmask_out, genome, expansion=merge, ori=True)
    print("Getting sequences is done")


# GetORF part


# def rps_blast(in_file, cdd, e_value=0.01, threads=1, outprefix="matches_w_hits", xml_prefix="hitdata", stat=True, check=False):
#     """
#     runs rpsblast on sequences from GetSeq; the output set to xml.
#     then the script gets all CD hits and writes them to a file. Also plots a histogram of hit's lengths distribution
#     :param in_file: input file for RPS-BLAST - fasta file from GetSeq
#     :param cdd: path to local Conserved Domain Database (CDD), including index name
#     :param e_value: expectation value (E), default 0.01
#     :param threads: number of threads to use, default 1
#     :param outprefix: prefix of the output file name
#     :param xml_prefix: prefix of xml file (rpsblastn output)
#     :param stat: make descriptive statistics (only when check=False)
#     :param check: make just ORFs checking
#     :return: no return
#     """
#     # add MGE's name to any output
#     mge_name = read_config("RepeatType") + "_"
#     if not check:
#         outprefix = mge_name + outprefix
#
#     xml_prefix = mge_name + xml_prefix
#
#     print("run RPS-BLAST...")
#     # print("CDD location set as %s" % cdd)
#     xml_file = xml_prefix + str(e_value)[2:] + ".xml"
#
#     # for rpsblast 2.6.0+
#     # rpsblast -db RTCDD/RT -query excised_matches.fa -out test.fmt5 -evalue 0.01 -outfmt 5 -num_threads 12
#     # sbp.call(["rpsblast", "-i", in_file, "-p", "F", "-d", cdd,  "-e", str(e_value), "-m 7", "-l", "log", "-a", str(threads), "-o", xml_file])
#
#     sbp.call(["rpstblastn", "-db", cdd, "-query", in_file, "-out", xml_file, "-evalue", str(e_value), "-outfmt", "5", "-num_threads", str(threads)])
#     if not check:
#         # add to config only if it is not ORF checking
#         add_config(unit="GetORF Output", val=xml_file)
#
#     try:
#         hits_only = [item.query for item in NCBIXML.parse(open(xml_file)) if item.alignments]
#         if len(hits_only) == 0:
#             print("No sequences with hits found! Exit.")
#             sys.exit()
#         else:
#             print("Found %i sequences with hits." % len(hits_only))
#     except FileNotFoundError:
#         print("ERROR! File %s not found!" % xml_file)
#         sys.exit()
#
#     # make dictionary excised_match:SeqRecord
#     dna_dict = dict()
#     for record in SeqIO.parse(in_file, "fasta"):
#         dna_dict[record.description] = record
#
#     # collect dna sequences with hits to a new list
#     # note that last 3 symbols are cropped
#     hits_only_dna = [dna_dict[item] for item in hits_only]
#
#     if stat:
#         # make some stat and plots (not for check only)
#         hits_len = pd.Series([len(item) for item in hits_only_dna])
#         plot = hits_len.hist(color='g', alpha=0.5, bins=200)
#         fig = plot.get_figure()
#         fig.savefig(outprefix + "_e" + str(e_value)[2:] + ".png")
#         plt.close(fig)
#         des_stat = hits_len.describe()
#         # write statistics to a file
#         des_stat.to_csv(outprefix + "_e" + str(e_value)[2:] + ".stats", header=False, index=True, sep="\t")
#         # print on the screen
#         print("Statistics of excised matches' length:")
#         print(pd.DataFrame(des_stat))
#
#     # make output name
#     # if runs as part of check_orfs()
#     if not check:
#         matches_fname = outprefix + "_e" + str(e_value)[2:] + ".fa"
#         # add to config matches with hits, but not checked ORFs!
#         add_config("ORFinder Input", matches_fname)
#     else:
#         matches_fname = outprefix + ".fa"
#
#     with open(matches_fname, 'w') as handle:
#         SeqIO.write(hits_only_dna, handle, 'fasta')
#     # print("Done!")

def rps_blast(in_file, cdd, smp_table='', e_value=0.01, threads=1, outprefix="matches_w_hits", xml_prefix="hitdata", stat=True, check=False):
    """
    runs rpsblast on sequences from GetSeq; the output set to xml.
    then the script gets all CD hits and writes them to a file. Also plots a histogram of hit's lengths distribution
    :param in_file: input file for RPS-BLAST - fasta file from GetSeq
    :param cdd: path to local Conserved Domain Database (CDD), including index name
    :param smp_table: csv table of SMP files groupings
    :param e_value: expectation value (E), default 0.01
    :param threads: number of threads to use, default 1
    :param outprefix: prefix of the output file name
    :param xml_prefix: prefix of xml file (rpsblastn output)
    :param stat: make descriptive statistics (only when check=False)
    :param check: make just ORFs checking
    :return: no return
    """
    # add MGE's name to any output
    mge_name = read_config("RepeatType") + "_"
    if not check:
        outprefix = mge_name + outprefix

    xml_prefix = mge_name + xml_prefix

    print("run RPS-BLAST...")
    # print("CDD location set as %s" % cdd)
    xml_file = xml_prefix + str(e_value)[2:] + ".xml"

    # for rpsblast 2.6.0+
    # rpsblast -db RTCDD/RT -query excised_matches.fa -out test.fmt5 -evalue 0.01 -outfmt 5 -num_threads 12
    # sbp.call(["rpsblast", "-i", in_file, "-p", "F", "-d", cdd,  "-e", str(e_value), "-m 7", "-l", "log", "-a", str(threads), "-o", xml_file])

    # read a table of CD groupings (assume it is tab delimited)
    # cdgroups = "doms.txt"
    if smp_table == '':
        try:
            reader = csv.reader(open(read_config('cd_table')))
        except KeyError:
            print("No smp table found in the config file!")
            smptable = input("Type a path to smp table > ")
            reader = csv.reader(open(smptable))
    else:
        reader = csv.reader(open(smp_table))
    group_tab = [row for row in reader]
    group_label = set([row[1] for row in group_tab])
    group_dict = {}
    for l in group_label:
        group_dict[l]=[r[0] for r in group_tab if l in r[1]]


    sbp.call(["rpstblastn", "-db", cdd, "-query", in_file, "-out", xml_file, "-evalue", str(e_value), "-outfmt", "5", "-num_threads", str(threads)])
    if not check:
        # add to config only if it is not ORF checking
        add_config(unit="GetORF Output", val=xml_file)

    try:
        #hits_only = [item.query for item in NCBIXML.parse(open(xml_file)) if item.alignments]
        hits_only = []
        for record in NCBIXML.parse(open(xml_file)):
            if record.alignments:
                hits = []
                for align in record.alignments:
                    hits.append(align.title)
                # leave only the name of a CD file that produced a hit
                # based on assumption that all 'title' start from 'gnl|CDD|123456 '
                hits = [re.findall('^gnl\|CDD\|[0-9]* ([a-z,A-Z]*[0-9]*), ', hit)[0] for hit in hits]
                # now I need to check for each CD label if any of its SMP present in 'hits' list
                hits_counter = 0
                for key in group_dict.keys():
                    for smp in group_dict[key]:
                        if smp[:-4] in hits:
                            hits_counter += 1
                            break
                    continue

                if hits_counter == len(group_dict.keys()):
                    # all CDs were found in this record!
                    hits_only.append(record.query)

        if len(hits_only) == 0:
            print("No sequences with specified number of hits found!\nTry to search for one CD.\nExit.")
            sys.exit()
        else:
            print("Found %i sequences with %i hits." %(len(hits_only), len(group_label)))
    except FileNotFoundError:
        print("ERROR! File %s not found!" % xml_file)
        sys.exit()

    # make dictionary excised_match:SeqRecord
    dna_dict = dict()
    for record in SeqIO.parse(in_file, "fasta"):
        dna_dict[record.description] = record

    # collect dna sequences with hits to a new list
    # note that last 3 symbols are cropped
    hits_only_dna = [dna_dict[item] for item in hits_only]

    if stat:
        # make some stat and plots (not for check only)
        hits_len = pd.Series([len(item) for item in hits_only_dna])
        plot = hits_len.hist(color='g', alpha=0.5, bins=200)
        fig = plot.get_figure()
        fig.savefig(outprefix + "_e" + str(e_value)[2:] + ".png")
        plt.close(fig)
        des_stat = hits_len.describe()
        # write statistics to a file
        des_stat.to_csv(outprefix + "_e" + str(e_value)[2:] + ".stats", header=False, index=True, sep="\t")
        # print on the screen
        print("Statistics of excised matches' length:")
        print(pd.DataFrame(des_stat))

    # make output name
    # if runs as part of check_orfs()
    if not check:
        matches_fname = outprefix + "_e" + str(e_value)[2:] + ".fa"
        # add to config matches with hits, but not checked ORFs!
        add_config("ORFinder Input", matches_fname)
    else:
        matches_fname = outprefix + ".fa"

    with open(matches_fname, 'w') as handle:
        SeqIO.write(hits_only_dna, handle, 'fasta')


def orf_finder(input_file, s=0, g=1, l=1000, strand="plus", nested=True):
    """
    runs ORFfinder program with specified parameters.
    defaults are minimum length 1000 nt, plus strand and ignore nested ORFs
    :param input_file: fasta file with sequences that have CDD hits
    :param l: minimum length of ORF
    :param strand: output ORFs on specified strand only (both|plus|minus)
    :param nested: ignore nested ORFs (completely placed within another)
    :param s: ORF start codon to use. 0 = "ATG" only; 1 = "ATG" and alternative initiation codons; 2 = any sense codon; Default = 0
    :param g: genetic code to use (1-31, Default 1). See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.
    :return: no return, function writes to a file instead
    """
    print("run ORFinder...")

    # get MGE name
    mge_name = read_config("RepeatType") + "_"
    output_prefix = mge_name + "cds" + str(l)
    add_config("ORFinder Output", output_prefix+".fa")

    if input_file == '':
        input_file = read_config("ORFinder Input")

    if nested:
        n = "true"
    else:
        n = "false"
    stderr = open("ORFfinder.stderr", "w")
    sbp.call(["ORFfinder", "-in", input_file, "-s", str(s), "-g", str(g), "-ml", str(l), "-n", n, "-strand", strand, "-out", output_prefix + ".fa", "-outfmt", "1"], stderr=stderr)
    stderr.close()


def check_orfs(l, cdd, cd_table='', threads=1, eval=0.01):
    """
    a function to check ORFs; it runs rps_blast() function with some different parameters to check already found ORFs
    :param l: minimum length of ORF (should be the same as in orf_finder() function)
    :param cdd: a path lo local conserved domain database (should be the same as in orf_finder() function)
    :param threads: number of threads to use
    :return: nothing
    """
    # infilename = "cds" + str(l) + ".fa"
    infilename = read_config("ORFinder Output")  # RepType_cds1000.fa
    domain = read_config("domain")
    output_prefix = infilename[:-3] + "_with_" + domain + "_e" + str(eval)[2:]
    print("Checking found ORFs with RPS-Blast...")
    rps_blast(infilename, cdd, smp_table=cd_table, e_value=eval, threads=threads, outprefix=output_prefix, xml_prefix="orf_hitdata", stat=True, check=True)
    add_config("ORFinder Output Checked", output_prefix + ".fa")


def get_orf(sequence, cd_data, min_len, orf_in, strnd, start_codon, cd_tab='', eval=0.01, threads=1, gencode=1):
    """
    main function
    :param sequence: input file for RPS-BLAST - fasta file from GetSeq
    :param cd_data: path to local Conserved Domain Database (CDD), including index name
    :param min_len: minimum length of ORF to search
    :param orf_in: fasta file with sequences that have CDD hits
    :param strnd: output ORFs on specified strand only (both|plus|minus)
    :param start_codon: ORF start codon to use. 0 = "ATG" only; 1 = "ATG" and alternative initiation codons; 2 = any sense codon; Default = 0
    :param eval: expectation value (E), default 0.01
    :param threads: number of threads to use
    :param gencode: genetic code to use (1-31, Default 0). See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.
    :return:
    """
    rps_blast(in_file=sequence, e_value=eval, cdd=cd_data, smp_table=cd_tab, threads=threads)
    orf_finder(orf_in, s=start_codon, l=min_len, strand=strnd, g=gencode)
    check_orfs(l=min_len, cdd=cd_data, cd_table=cd_tab, threads=threads, eval=eval)

# Extender part


def make_bed(cds_records, right_expand=500, left_expand=500):
    """
    a function to make a list of lines formatted for bed file
    :param cds_records: a list of coding sequences
    :param right_expand: elongation value for right (3') flanking region
    :param left_expand: elongation value for left (5') flanking region
    :return: a list of lines formatted for bed file
    """
    # collect length of scaffolds
    genome_file = read_config("genome")
    if ".gz" in genome_file:
        genome_file = genome_file[:-3]
    scaff_lengths = [(rec.id, len(rec.seq)) for rec in SeqIO.parse(open(genome_file), "fasta")]
    scaff_dict = dict()
    for item in scaff_lengths:
        scaff_dict[item[0]] = item[1]

    bed_lines_list = []
    for record in cds_records:
        header = record.id  # 'lcl|NC_018898.1:9535131-9536790():62-1489'
        scaff_name = re.findall(r'lcl\|(.*?):', header)[0]
        positions = re.findall(r':(.*):', header)  # ['9535131-9536790()']
        positions = positions[0][:-2].split("-")
        match_left, match_right = positions[0], positions[1]  # '9535131' & '9536790'

        # find start and end of an ORF
        orf_position = header.split(":")[2]
        start, end = int(orf_position.split("-")[0]) + int(match_left), int(orf_position.split("-")[1]) + int(match_left)

        left_flank = start - left_expand
        right_flank = end + right_expand

        if left_flank < 0:
            left_flank = 0
        if right_flank > scaff_dict[scaff_name]:
            right_flank = scaff_dict[scaff_name]

        # scaffold name
        scaffold_name = re.findall(r'lcl\|(.*?):', header)[0]  # search in non-greedy manner

        # make a list with each item representing a line in a bed file
        # the question is what strand to choose :(
        line = scaffold_name + "\t" + str(left_flank) + "\t" + str(right_flank) + "\t" + "+"
        bed_lines_list.append(line)
    return bed_lines_list


def extender(genome_file, cds_file, left_elongation_value, right_elongation_value):
    """

    :param genome_file:
    :param cds_file:
    :param left_elongation_value:
    :param right_elongation_value:
    :return:
    """
    if left_elongation_value == 0 and right_elongation_value == 0:
        # then do nothing
        pass
    else:
        cds_w_rt = [rec for rec in SeqIO.parse(open(cds_file), "fasta")]
        bed_lines = make_bed(cds_w_rt, left_expand=left_elongation_value, right_expand=right_elongation_value)

        # write bed formatted file
        output_prefix = "extended_L" + str(left_elongation_value) + "R" + str(right_elongation_value)
        output_bed = output_prefix + ".bed"
        with open(output_bed, "w") as bed:
            for item in bed_lines:
                bed.write("%s\n" % item)
        if ".gz" in genome_file:
            genome_file = genome_file[:-3]
        # excision process
        sbp.call(["bedtools", "getfasta", "-s", "-fi", genome_file, "-bed",  output_bed, "-fo", cds_file[:-3] + "_" + output_prefix + ".fa"])


# RepeatModeler & RepeatMasker part


def rmodeler(filename, cpus):
    # get its source, e.g. absolute path to it
    filename_src = os.getcwd() + "/" + filename

    # create a directory for RepeatModeler runs
    if not os.path.isdir("RepModOut"):
        os.mkdir("RepModOut")
        # and go into it
        os.chdir("RepModOut")
    else:
        # and go into it
        os.chdir("RepModOut")

    # get destination for filename to move it
    filename_dest = os.getcwd() + "/" + filename

    # move it to CWD
    try:
        os.rename(filename_src, filename_dest)
    except FileNotFoundError:
        print("ERROR! File %s not found!" % filename_src)
        sys.exit()

    # check if the file is gzipped
    if ".gz" in filename:
        # print("File %s seems to be gzipped..." % filename)
        # gunzip the filename
        os.system("gunzip %s" % filename)
        ungzipped_filename = filename[:-3]
    else:
        ungzipped_filename = filename

    # rename genome's scaffolds properly (otherwise RepeatModeler will crash)
    cmd = "awk '/^>/{print \">sc\" ++i; next}{print}' < %s > ref.fa" % ungzipped_filename
    os.system(cmd)

    # build database
    cmd = "%s -name %s.db -engine ncbi ref.fa" % (read_config("BuildDatabase"), ungzipped_filename)
    os.system(cmd)

    # run RepeatModeler
    print("run RepeatModeler on %s CPUs..." % cpus)
    cmd = "%s -engine ncbi -pa %s -database %s.db > RepMod.out" % (read_config("RepeatModeler"), cpus, ungzipped_filename)
    os.system(cmd)

    # move file upwards (where it was initially)
    os.system("mv %s ../" % ungzipped_filename)
    # final message
    print("RepeatModeler finished")
    # find output file path, move it to main working directory and add it to the config
    repmod_outfile = glob.glob("RM*/consensi.fa.classified")[0]
    # copy consensi file to CWD
    os.system("cp %s ../" % repmod_outfile)
    print("consensi.fa.classified copied to the current working directory")
    # go up! to initial working dir, where the whole pipeline was started
    os.chdir("../../")


def rmasker(filename, cpus, lib):
    """
    function to run RepeatMasker
    :param filename: name of a genome assembly file
    :param cpus: number of CPUs to use
    :return: no return
    """

    # genome file is already ungzipped
    # but it's name still contains 'gz' suffix
    filename = filename[:-3]
    mge_name = read_config("RepeatType")
    exe = read_config("RepeatMasker")
    # lib = "All_" + read_config("RepeatType") + "_consensi.fa"
    out = mge_name + "_" + filename + "_RMout.txt"

    # check if the genome exists
    if not os.path.isfile(filename):
        print("ERROR! File %s not found in the working directory %s" % (filename, os.getcwd()))
        sys.exit()

    # check if custom library file exists
    if not os.path.isfile(lib):
        print("ERROR! File %s not found in the working directory %s" % (lib, os.getcwd()))
        sys.exit()

    # the actual bash command is as follows:
    # ~/RepeatMasker/RepeatMasker -lib $species"_AllPLEconsensi.fa" -pa 24 -s -no_is -nolow $species".fa" -u -gff > $species"RepMask.stdout" &&
    cmd = "%s -lib %s -pa %s -s -no_is -nolow %s -u -gff > %s" % (exe, lib, cpus, filename, out)
    print("Full command:\n%s" % cmd)
    os.system(cmd)
    # Sadly RepeatMasker itself names the output files according to the scheme: 'genome.fna + out*' or 'genome.fna + ori.out*'
    # these several files are (plus to the genome basename): *.fai, *.masked, *.out, *.out.gff, *.tbl, *ori.out
    # find these files and rename them
    rm_files_list = ["%s.masked" % filename, "%s.out" % filename, "%s.out.gff" % filename, "%s.tbl" % filename, "%s.ori.out" % filename]
    rm_files = list(map(lambda x: glob.glob(x), rm_files_list))
    n = 0
    for item in rm_files:
        try:
            new_name = mge_name + "_" + item[0]
            cmd = "mv %s %s" % (item[0], new_name)
            os.system(cmd)
            n += 1
        except IndexError:
            print("Warning %s not found! Cannot rename." % rm_files_list[n])
            n += 1

    add_config("RepeatMasker Output",  "%s.out" % (mge_name + "_" + filename))
    # final message
    print("RepeatMasker finished")


def translator(seq_obj, gcode=1):
    """
    function to translate the CDS
    :param seq_obj: a list of Seq objects from SeqIO.parse
    :param gcode: genetic code to use (have to be the same as in ORFfinder; default=1
    :return: a list of translated seqs
    """
    # recode NCBI's names of genetic code tables to Biopython's
    gencode_table ={'Alternative Flatworm Mitochondrial': 14, 'Alternative Yeast Nuclear': 12, 'Archaeal': 11, 'Ascidian Mitochondrial': 13, 'Bacterial': 11, 'Blepharisma Macronuclear': "",
     'Candidate Division SR1': 25, 'Chlorophycean Mitochondrial': 16, 'Ciliate Nuclear': 6, 'Coelenterate Mitochondrial': 4, 'Dasycladacean Nuclear': 6,
     'Echinoderm Mitochondrial': 9, 'Euplotid Nuclear': 10, 'Flatworm Mitochondrial': 9, 'Gracilibacteria': 25, 'Hexamita Nuclear': "", 'Invertebrate Mitochondrial': 5,
     'Mold Mitochondrial': 4, 'Mycoplasma': 4, 'Pachysolen tannophilus Nuclear Code': 26, 'Plant Plastid': 11, 'Protozoan Mitochondrial': 4, 'Pterobranchia Mitochondrial': 24,
     'Scenedesmus obliquus Mitochondrial': 22, 'Spiroplasma': 4, 'Standard': 1, 'Thraustochytrium Mitochondrial': 23, 'Trematode Mitochondrial': 21, 'Vertebrate Mitochondrial': 2,
     'Yeast Mitochondrial': 3}
    # reverse the dict
    gencode_table_inv = {v: k for k, v in gencode_table.items()}
    # recode
    biopython_code = gencode_table_inv[gcode]
    # translate CDS
    translated_seq = [SeqRecord(rec.seq.translate(table=biopython_code, to_stop=True), id=rec.id) for rec in seq_obj]
    return translated_seq


# The MGERT pipeline itself


def pipe(genome_file, mge_type, dom_table="", lib="", rm_tab="", seq_for_dom='', stage=1, threads=1, censor=False, pandas=False, ori=False, merge=0, l=1000, e=0.01, c=0, strnd="plus", g=1, le=500, re=500):
    """
    a function to run the whole pipeline
    :param genome_file: a genome assembly file
    :param mge_type: type of mobile element to search
    :param dom_table: csv table with smp files and groupings
    :param lib: fasta file of consenus sequences of MGE
    :param rm_tab: a repeat masker table
    :param seq_for_dom: file with sequences too search for domains
    :param stage: the step of the analysis to start with. 1 - rmodeler; 2 - get_cons; 3 - get_seq; 4 - get_orf
    :param threads: number of threads to use (for RepeatModeler, RepeatMasker, ORFfinder, rpstblastn)
    :param censor: use censor output data or not; default False
    :param pandas: use pandas or not; default False
    :param ori: use ori output or not; default False
    :param merge: merge all repeats within M bp into a single entry
    :param l: minimum length of ORF to report
    :param e: minimum e-value of rpstblastn output
    :param c: ORF start codon to use. 0 = "ATG" only; 1 = "ATG" and alternative initiation codons; 2 = any sense codon; Default = 0
    :param strnd: output ORFs on specified strand only (both/plus/minus)
    :param g: genetic code to use (1-31, Default 0). See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.
    :param le:
    :param re:
    :return: no
    """
    # check if the genome file has right extensions
    if ".fna.gz" in genome_file:
        # genome assembly is gzipped
        # add genome name to the config to use it in the following steps
        add_config("genome", genome_file)
        # we assume, that the genome file is called `genome.fna.gz`, so we can derive the prefix by cropping last seven characters
        add_config("prefix", genome_file[:-7])
        dirname = read_config("prefix")
        # make name for RepeatMasker output
        add_config("RepeatMasker Output", genome_file[:-3] + ".out")
    elif ".fna" in genome_file and ".gz" not in genome_file:
        add_config("prefix", genome_file[:-4])
        add_config("genome", genome_file + ".gz")
        dirname = read_config("prefix")
        add_config("RepeatMasker Output", genome_file + ".out")
    else:
        print("Error! Seems the genome assembly has wrong name.\nCheck if it has 'fna' or 'fna.gz' extension")
        sys.exit()

    if stage == 1:
        if not os.path.isdir(dirname):
            # make dir with name like taken from 'prefix' in config
            os.mkdir(dirname)
        # move genome to this new directory
        os.system("mv %s %s" % (genome_file, dirname))
        os.chdir(dirname)

        # run RepeatModeler
        print("1/5. Start RepeatModeler pipeline...\nNumber of CPUs - %s" % str(threads))
        rmodeler(genome_file, threads)
        # when rmodeler finished, it chdir back
        stage += 1  # stage = 2
        # print("%s search started..." % mge_type)

    # add MGE type to the config
    add_config(unit="RepeatType", val=mge_type)

    if stage == 2:
        # change "RepeatModeler Output" in the config
        # add_config(unit="RepeatModeler Output", val=mge_type + "_consensi.fa.classified")
        # go to 'prefix' dirname again
        try:
            os.chdir(dirname)
        except FileNotFoundError:
            print("...")

        if lib == "":
            # no library specified - collect repeats consensi
            # next step -> GetCons
            print("2/5. Collecting user defined consensuses...")
            get_cons(mge_type=mge_type, standard=True)
            if censor:
                url = input("Censor output URL > ")
                get_cons(censor=True, url=url)
                get_cons(mge_type=mge_type, recollect=True)
            else:
                # if no censor, then no recollect
                # rename initial file into final one
                # reptype = read_config("RepeatType")
                cmd = "mv " + mge_type + "_consensi.fa.classified" + " All_" + mge_type + "_consensi.fa"
                os.system(cmd)
            print("All consensi found.\nStart RepeatMasker with the consensi as a library...\nNumber of CPUs - %s" % str(threads))
            # run RepeatMasker
            # if consensi library not specified, take it from RepeatModeler
            lib = "All_" + read_config("RepeatType") + "_consensi.fa"
            rmasker(genome_file, str(threads), lib=lib)
        else:
            # if user specified the library use it instead
            rmasker(genome_file, str(threads), lib=lib)
        stage += 1

    if stage == 3:

        try:
            os.chdir(dirname)
        except FileNotFoundError:
            print("...")

        # run GetSeq
        print("3/5. Defining coordinates in the assembly...")

        if pandas and not ori:
            get_seq(pandas=True, merge=merge, rm_tab=rm_tab)
            # os.system("GetSeq.py -p -m %s" % str(merge))
        elif pandas and ori:
            get_seq(pandas=True, ori=True, merge=merge, rm_tab=rm_tab)
            # os.system("GetSeq.py -p -o -m %s" % str(merge))
        elif not pandas and ori:
            get_seq(pandas=False, ori=True, merge=merge,rm_tab=rm_tab)
            # os.system("GetSeq.py -o -m %s" % str(merge))
        elif not pandas and not ori:
            get_seq(pandas=False, ori=False, merge=merge,rm_tab=rm_tab)
            # os.system("GetSeq.py -m %s" % str(merge))
        stage += 1

    if stage == 4:

        try:
            os.chdir(dirname)
        except FileNotFoundError:
            print("...")

        print("4/5. Collecting ORFs...")
        # os.system("GetORF.py -l %s -e %s -g %s" % (str(l), str(e), str(g)))
        if seq_for_dom == '':
            seq = read_config("GetSeq Output")
        else:
            seq = seq_for_dom

        my_cdd = read_config("CDD")
        # ORFinder Input may be absent
        try:
            orf_infile = read_config("ORFinder Input")
        except KeyError:
            orf_infile = ''

        if dom_table == '':
            # no smp table provided
            print("You didn't provide smp table.\nTrying to look in the config file...")
            try:
                dom_table = read_config("cd_table")
                print("OK...")
            except KeyError:
                print("smp table not found in the config file.")
                dom_table = input("Please provide valid path to the file > ")


        get_orf(sequence=seq, cd_data=my_cdd, cd_tab=dom_table, orf_in=orf_infile, min_len=l, eval=e, start_codon=c, strnd=strnd, threads=threads, gencode=g)
        # from here there should be an entry "ORFinder Output Checked" in config
        # so retrieve its value
        checked_orfs = [rec for rec in SeqIO.parse(read_config("ORFinder Output Checked"), "fasta")]
        print("...translating ORFs...")
        translated_orfs = translator(checked_orfs, gcode=g)
        # write them; filename is made by adding "a" to the suffix of checked ORFs
        f = read_config("ORFinder Output Checked")+"a"
        SeqIO.write(translated_orfs, f, "fasta")
        stage += 1

    if stage == 5:
        try:
            os.chdir(dirname)
        except FileNotFoundError:
            print("...")

        print("5/5. Adding flanking regions...")
        cds = read_config("ORFinder Output Checked")
        extender(genome_file=genome_file, cds_file=cds, left_elongation_value=le, right_elongation_value=re)


if __name__ == '__main__':
    my_message = "Mobile Genetic Elements Retrieving Tool (MGERT).\n" \
                 "MGERT is a pipeline for retrieving coding sequences of mobile genetic elements from genomic assemblies\n" \
                 "To make initial configuration run configuration script with: `MGERT.py --configure`\n"

    parser = argparse.ArgumentParser(description=my_message,  formatter_class=argparse.RawTextHelpFormatter, prog="MGERT", usage='%(prog)s -a genome.fna.gz -T Penelope [options]')
    parser.add_argument("--configure", action="store_true", help="run the configuration script")
    parser.add_argument("--make_cdd", action="store_true", help="make local CDD")
    parser.add_argument("-cd", "--cd_table", type=str, metavar="[domains.csv]", help="comma delimited file with smp files and their groupings.", default="")
    parser.add_argument("-a", "--assembly", type=str, metavar="[genome.fa.gz]", help="specify a genome assembly file", default="")
    parser.add_argument("-T", "--mge_type", type=str, metavar="[Penelope/BovB/RTE/CR1/L1/LINE etc]", help="specify the type of MGE to search", default="")
    parser.add_argument("-f", "--from_stage", type=str, metavar="[cons/coords/orfs/flanks]", help="specify the step from which the pipeline should start. 'consensus' - get consensus sequences; "
                                                                                 "'coords' - get sequences; "
                                                                                 "'orfs' - get ORFs; flanks - extend CDS.\nDefault 'rmod'", default="rmod")
    #parser.add_argument("-k", "--check_types", action="store_true", help="Print out all the types of MGE found in the RepeatModeler output")
    parser.add_argument("-k", "--check_types", type=str, metavar="[consensi file]", help="Print out all the types of MGE found in the RepeatModeler output")
    parser.add_argument("-t", "--threads", type=int, metavar="[integer]", help="set number of threads", default=6)
    parser.add_argument("-C", "--censor", type=str, metavar="[URL]", help="pass URL of Censor classification results to MGERT")
    parser.add_argument("-o", "--ori", action="store_false", help="if specified MGERT will use the *.ori file to fetch the coordinates instead of *_rm.out file", default=False)
    # parser.add_argument("-d", "--pandas", action="store_false", help="run GetSeq with pandas (can be very slow)", default=False)
    parser.add_argument("-m", "--merge", type=int, metavar="M", help="merge all hits within M bp into a single entry", default=500)
    parser.add_argument("-e", "--e_value", type=float, metavar="[real]", help="set expectation value (E), default 0.01", default=0.01)
    parser.add_argument("-c", "--start_codon", type=int, metavar="[integer]", help="ORF start codon to use. 0 = 'ATG' only; 1 = 'ATG' and alternative initiation "
                                                                                                  "codons; 2 = any sense codon; Default = 0", default=0)
    parser.add_argument("-l", "--min_length", type=int, metavar="[integer]", help="set minimum length of ORF, default 1000 nt", default=1000)
    parser.add_argument("-s", "--strand", type=str, metavar="[plus/minus/both]", help="output ORFs on specified strand only. Default 'plus'", default="plus")
    parser.add_argument("-g", "--genetic_code", type=int, metavar="[integer]", help="genetic code to use (1-31, Default 1). "
                                                                                    "See http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.", default=1)
    parser.add_argument("-le", "--left_end", type=int, metavar="[500]", help="length of ORFs' left flanking region. Default 500 bp", default=500)
    parser.add_argument("-re", "--right_end", type=int, metavar="[500]", help="length of ORFs' right flanking region. Default 500 bp", default=500)
    parser.add_argument("-L", "--lib", type=str, metavar="[fasta file]", help="library for RepeatMasker (fasta format). Use with `-f cons` only.\n"
                                                                              "When consensus library is not specified it will be automatically composed from RepeatModeler output", required=False, default="")
    parser.add_argument("-rm", "--rm_table", type=str, metavar="[RepeatMasker table]",
                        help="specify repeat masker table to use, default none. Use with `-f coords` option only", required=False, default="")
    parser.add_argument("-sq", "--sequence", type=str, metavar="[sequence.fasta]",
                        help="specify file name of sequences where to look for domains. Use with `-f orf` option only", required=False, default="")
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 0.3.9')

    args = parser.parse_args()

    # if no args have been passed, then there will be error message and short description of the available options
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
        # if one specified configuration

    if args.check_types:
        check_types(seq_file=args.check_types)
        sys.exit()

    if args.configure:
        make_config()
        sys.exit()
    elif args.make_cdd:
        make_local_cdd(dir_for_cd="LocalCDD")
        sys.exit()

    # check if user didn't specify genome assembly and MGE type
    if args.assembly == "":
        genome = input("Specify genome assembly file > ")
    else:
        genome = args.assembly

    if args.mge_type == "":
        mge = input("Specify MGE type to search [Penelope/BovB/RTE/CR1/L1/LINE etc] > ")
    else:
        mge = args.mge_type

    # recode the stage parameter
    if args.from_stage == "rmod":  # default
        stage = 1
    elif args.from_stage == "cons" or args.from_stage == "consensus":
        stage = 2
    elif args.from_stage == "coords":
        stage = 3
    elif args.from_stage == "orfs" or args.from_stage == "orf":
        stage = 4
    elif args.from_stage == "flanks":
        stage = 5
    else:
        print("Wrong stage name!")
        sys.exit()

    if args.censor:
        # print("Standard mode, with CENSOR")
        # get genome file
        # genome = read_config("genome", home=True)
        pipe(genome_file=genome, mge_type=mge, stage=stage, seq_for_dom=args.sequence, threads=args.threads, censor=True, ori=args.ori, merge=args.merge,
             l=args.min_length, e=args.e_value, c=args.start_codon, strnd=args.strand, g=args.genetic_code, le=args.left_end, re=args.right_end, rm_tab=args.rm_table)
    elif args.lib:
        # print("Standard mode, with specified library")
        # get genome file
        # genome = read_config("genome", home=True)
        pipe(genome_file=genome, mge_type=mge, stage=2, seq_for_dom=args.sequence, threads=args.threads, censor=False, ori=args.ori, merge=args.merge,
             l=args.min_length, e=args.e_value, c=args.start_codon, strnd=args.strand, g=args.genetic_code, le=args.left_end, re=args.right_end, lib=args.lib, rm_tab=args.rm_table)

    else:
        # print("Standard mode")
        # get genome file
        # genome = read_config("genome", home=True)
        pipe(genome_file=genome, mge_type=mge, stage=stage, seq_for_dom=args.sequence, threads=args.threads, censor=False, ori=args.ori, merge=args.merge,
             l=args.min_length, e=args.e_value, c=args.start_codon, strnd=args.strand, g=args.genetic_code, le=args.left_end, re=args.right_end, rm_tab=args.rm_table)
