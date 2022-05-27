#!/usr/bin/env python
# coding: utf-8




import requests
import os
import os.path
import io
import subprocess
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from numpy import median
from collections import Counter
import dataframe_image as dfi
import argparse
import sys
from Bio import Phylo
from Bio.Phylo.Consensus import *
import datetime
e = datetime.datetime.now()
#paht to scrpit location

here = os.getcwd()


# In[4]:


#create file if doesnt exist 
def create_file(path):
    exist = os.path.exists(path)
    if not exist:
        os.makedirs(path)


# In[5]:


def make_list_from_file(file):
    f = open(file)
    g = f.read().split()
    zw = []
    pom = []
    
    for j in g:
        if len(j.split()) > 0:
            pom.append(j)
    #print(pom)
    for i in range(0, len(pom)-1, 2):
        zw.append(pom[i] + ' ' + pom[i+1])
    #print(zw)
    
    return zw


# In[6]:


#create file with names of organism in row that are avaliable in uniport
#create file with names and them proteom id
def download_id_proteoms_by_taxa():
    url = requests.get('https://www.uniprot.org/proteomes/?query=taxonomy%3AAves&sort=score&format=tab')
    a = url.content.decode()
    data = io.StringIO(a)
    table = pd.read_csv(data, sep = '\t')

    s = {}
    
    for i in range(len(table)):
        #take name as 2 words
        name_org = table.loc[i, 'Organism']
        name_org = name_org.split()
        name_org = ' '.join(map(str, name_org[:2]))
        
        #check size of proteom
        size  = table.loc[i, 'Protein count']
        
        if size >= 5000 and name_org not in s.keys():
            s[name_org] = table.loc[i, 'Proteome ID']
            
    ava_name = list(s.keys())
    ava_name.sort()
    
    with open(here + "/avaliable_organism", 'w') as f, open(here + "/proteom_id", 'w') as f1:
        for i in ava_name:
            f.write(i + '\n')
            f1.write(i + '\t' + s[i] + '\n')


# In[4]:


def download_proteom(proteom_id, path_to_save):
    url = 'https://www.uniprot.org/uniprot/?query=proteome:' + proteom_id + '&format=fasta&compress=no'
    re = requests.get(url) 
    d = re.content.decode()
    with open(path_to_save, 'w') as f:
        f.write(d)


# In[8]:


#dowload proteoms for organism 
def dowload_proteoms_from_names_list(l = [], file_to_save = here + '/POB_sequences'):
    create_file(file_to_save)
    
    f = pd.read_csv('proteom_id', delimiter = '\t', header = None)
    #print(f)
    name_id = dict(zip(list(f.iloc[:,0]), list(f.iloc[:,1]) ))
    wrong = []
    good = []
    #print(name_id.keys())
    for i in l:
        if i not in name_id.keys():
            wrong.append(i)
            continue
            
        #name = i.replace(' ', '_')
        name = i
        print('Downloading proteom for ' +  i + '...' )
        download_proteom(name_id[i], file_to_save + '/' + name  )
        good.append(name)
    
    if len(wrong) > 0:
        print('Could not download proteom(s) for: ' + ', '.join(wrong) + '.' ) 
        print('Please check correctness and try again.')
        return [wrong, good]
    else:
        print('Done')
        return ['Done', good]


# In[9]:


def download_from_list(names_list):
    ready = []
    pom = dowload_proteoms_from_names_list(names_list)
    ready.extend(pom[1])
    while pom[0] != 'Done':
        to_do = []    
        for i in pom[0]:
            print('Write correct name instead ' + i)
            x = input()
            to_do.append(x)
        pom = dowload_proteoms_from_names_list(to_do)
        ready.extend(pom[1])
    #print(ready)

    with open(here + '/' + 'actual_analysis_organisms_names', 'w') as f:
        for i in ready:
            f.write(i + '\n')
    print('Created file with corrected, used names: ' + here + '/' + 'actual_analysis_organisms_names' )
    return ready


# In[10]:


#names_list=[
#    'Amazona guildingii',
#    'Aphelocoma coerulescens', 
#    'Bubo bubo']


# In[55]:

def auto_rep(list_of_names,path_to_save):
    merger(list_of_names,path_to_save)
    path_to_phyml = open("Path.txt","r")
    path_to=''
    for line in path_to_phyml:
        path_to = line
    os.system('mmseqs createdb DB.txt DB && mmseqs linclust DB DB_clu tmp && mmseqs createsubdb DB_clu DB DB_clu_rep && mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta && mmseqs createtsv DB DB DB_clu DB_clu.tsv && mmseqs createseqfiledb DB DB_clu DB_clu_seq && mmseqs result2flat DB DB DB_clu_seq DB_clu_seq.fasta')
    x=0
    while x==0:
        if os.path.exists("DB_clu_seq.fasta"):
            x=1
            with open("DB_clu_seq.fasta", 'r') as file:
                clusters_fasta = {}
                for line in file:
                    if '>' in line and '|' not in line:
                        clu = line[1:-1]
                        clusters_fasta[clu] = []
                    elif '>' in line and '|' in line:
                        clusters_fasta[clu].append({line.split('|')[1]: {
                            'org': line[line.index('OS=') + 3:line.index('OX=') - 1].replace(' ', '_'),
                            'seq': file.readline()[:-1]}})

            #print('liczba klastrow ', len(clusters_fasta))

            for c in clusters_fasta:
                if len(clusters_fasta[c]) > 2:
                    #print(f'w klastrze {c} jest {len(clusters_fasta[c])} sekwencji')
                    with open(f'{c}.cluster', 'w') as file:
                        for s in clusters_fasta[c]:
                            file.write('>' + s[list(s.keys())[0]]['org'] + '|' + list(s.keys())[0] + '\n' + s[list(s.keys())[0]][
                                'seq'] + '\n')
            os.system('for f in ./*.cluster; do muscle -in $f -phyiout ./$f.phypob; done && for f in ./*.phypob; do '+path_to+' -i $f -d aa; done && cat *.phy > Merge.txt')
            mc("Merge.txt")
            subprocess.run(['rm *.phy'])
            subprocess.run(['rm *.phypob'])
            subprocess.run(['rm *.cluster'])


def merger(list_of_names,path_to_save):
    if not os.path.exists(path_to_save):
        os.makedirs(path_to_save)
    file_name = "DB"
    completeName = os.path.join(path_to_save, file_name+".txt")
    seq_path = here + '/POB_sequences/'
    merge=''
    for i in list_of_names:
        with open(seq_path + i,"r") as file:
            for line in file:
                merge += line
    with open(completeName, "w") as f:
        f.write(merge)
        f.close()


# In[73]:


#klastry z powtorkami
def clusters_with_rep(file_path):
    with open(file_path, 'r') as file:
        clusters_fasta = {}
        for line in file:
            if '>' in line and '|' not in line:
                clu = line[1:-1]
                clusters_fasta[clu] = []
            elif '>' in line and '|' in line:
                clusters_fasta[clu].append({line.split('|')[1]: {
                    'org': line[line.index('OS=') + 3:line.index('OX=') - 1].replace(' ', '_'),
                    'seq': file.readline()[:-1]}})

    #print('liczba klastrow ', len(clusters_fasta))

    for c in clusters_fasta:
        if len(clusters_fasta[c]) > 2:
            #print(f'w klastrze {c} jest {len(clusters_fasta[c])} sekwencji')
            with open(f'{c}.cluster', 'w') as file:
                for s in clusters_fasta[c]:
                    file.write('>' + s[list(s.keys())[0]]['org'] + '|' + list(s.keys())[0] + '\n' + s[list(s.keys())[0]][
                        'seq'] + '\n')




# In[88]:


#klastry bez powtorek
def clusters_without_rep(file_path):
    print('Type your minimum species number per cluster')
    minimum_number = int(input())
    with open(file_path, 'r') as file:
        clusters_fasta = {}
        for line in file:
            if '>' in line and '|' not in line:
                clu = line[1:-1]
                clusters_fasta[clu] = []
            elif '>' in line and '|' in line:
                clusters_fasta[clu].append({
                    line.split('|')[1]: {'org': line[line.index('OS=') + 3:line.index('OX=') - 1].replace(' ', '_'),
                        'seq': file.readline()[:-1]}})

    #print('liczba klastrow ', len(clusters_fasta))

    for c in clusters_fasta:
        if len(clusters_fasta[c]) >= minimum_number:
            org_list = {}
            for s in clusters_fasta[c]:
                if s[list(s.keys())[0]]['org'] not in org_list:
                    # file.write(
                    #     '>' + s[list(s.keys())[0]]['org'] + '|' + list(s.keys())[0] + '\n' + s[list(s.keys())[0]][
                    #         'seq'] + '\n')
                    org_list[s[list(s.keys())[0]]['org']] = [s[list(s.keys())[0]]['org'], list(s.keys())[0],
                                                         s[list(s.keys())[0]]['seq']]
            if len(org_list) == minimum_number:
                with open(f'{c}.cluster', 'w') as file:
                    for k in org_list:
                        file.write('>' + org_list[k][0] + '|' + org_list[k][1] + '\n' + org_list[k][2] + '\n')
            #print(f'w klastrze {c} jest {len(org_list)} sekwencji')




# In[ ]:


##bootstrapping
def bootstrapping(input_file):
    print('Type your bootstrap support value')
    threshold=float(input())
    f1=open(input_file,'r')
    f2=open('bootstrap_output_file','w+')
    for i in f1:
        a=str(i)
        lista = a.split('/')
        wartosc = []
        for elem in lista:
            try:
                i = int(elem[0:2])
                wartosc.append(i)
            except ValueError:
                pass
        n=(sum(wartosc)/len(wartosc))
        if n >=threshold:
            f2.write(a+'\n')
    f1.close()
    f2.close()


# In[ ]:


#MC tree
def mc(file_path):
    print('Type your majority threshhold')
    majority_thresh = float(input())
    trees = list(Phylo.parse(file_path, "newick"))
    #strict_tree = strict_consensus(trees)
    #Phylo.write([strict_tree], "Strick_tree.phy", "newick")
    majority_tree = majority_consensus(trees, majority_thresh)
    Phylo.write([majority_tree], 'Majority.phy', "newick")

    mt1 = majority_tree

    fig = plt.figure(figsize=(20,30))
    axes = fig.add_subplot(1, 1, 1)
    plt.title('Drzewo MC', size = 20)
    Phylo.draw(mt1,  axes = axes)
    fig.savefig('MC_consensus_tree.png')


# In[15]:


def average_protein_length(file):
    poml = []
    record = SeqIO.parse(file, "fasta")
    record_list = list(record)

    for j in record_list:
        poml.append(len(j.seq))
    aver = sum(poml) / len(record_list)
    var = sum((x-aver)**2 for x in poml) / len(poml)
    
    return [aver, var]


# In[16]:


def make_average_table(list_of_names, path_to_save):
    seq_path = here + '/POB_sequences/'
    aver = []
    var = []
    for i in list_of_names:
        lists = average_protein_length(seq_path + i)
        aver.append(lists[0])
        var.append(lists[1])

    allposn = pd.DataFrame({'Name' : list_of_names, 'Average' : aver})
    dfi.export(allposn, path_to_save)


# In[17]:


def make_df(file, i):
    length = 0
    name = []
    le = []
    record = SeqIO.parse(file, "fasta")
    record_list = list(record)
    for j in record_list:
        le.append(len(j.seq))
        name.append(i)
    
    
    return [le, name]


# In[18]:


def make_average_plot(list_of_names, path_to_save):
    seq_path = here + '/POB_sequences/'
    name_list = []
    le = []
    for i in list_of_names:
        lists = make_df(seq_path + i, i)
        name_list+=lists[1]
        le+=lists[0]

    allposn = pd.DataFrame({'Name' : name_list, 'Length' : le})

    plt.figure(figsize=(10,6))
    sns.barplot(x = 'Name', y= 'Length', data = allposn, ci = 68, dodge = False, order =list(allposn.groupby(['Name']).mean().sort_values("Length").index))
    plt.xticks(rotation = 50, ha='right', size = 12, fontweight = 'bold')
    plt.title("Average protein length", size = 29)
    plt.xlabel('Organism name')
    plt.ylabel('Average length')
    plt.savefig(path_to_save, bbox_inches='tight', facecolor = 'white')


# In[19]:


def make_median_plot(list_of_names, path_to_save):
    seq_path = here + '/POB_sequences/'
    name_list = []
    le = []
    for i in list_of_names:
        lists = make_df(seq_path + i, i)
        name_list+=lists[1]
        le+=lists[0]

    allposn = pd.DataFrame({'Name' : name_list, 'Length' : le})
    
    plt.figure(figsize=(10,6))
    sns.barplot(x = 'Name', y= 'Length',estimator = median, ci = 68, data = allposn, dodge = False, order =list(allposn.groupby(['Name']).median().sort_values("Length").index))
    plt.xticks(rotation = 50, ha='right', size = 12, fontweight = 'bold')
    plt.title("Median protein length", size = 26)
    plt.xlabel('Organism name')
    plt.ylabel('Average length')
    plt.savefig(path_to_save, bbox_inches='tight', facecolor = 'white')


# In[20]:


def amino_acid_content(file):
    aa = ['C', 'D', 'S', 'Q', 'K', 'I', 'P', 'T', 'F', 'N', 'G', 'H', 'L', 'R', 'W', 'A', 'V', 'E', 'Y', 'M']
    pom = Counter()
    delete = []
    record = SeqIO.parse(file, "fasta")
    record_list = list(record)
    for j in record_list:
        pom += Counter(j.seq)
    pom = dict(pom)
    alll = sum(pom.values())
    for i in pom:
        pom[i] = (pom[i]/alll)*100
        if i not in aa:
            delete.append(i)
    for i in delete:
        pom.pop(i)
     
    return pom


# In[21]:


def amino_content_from_list(list_of_names, tf, path_out):
    seq_path = here + '/POB_sequences/'
    w = 0
    for i in list_of_names:
        w += 1
        k = amino_acid_content(seq_path + i)

        if w == 1:
            fi = pd.DataFrame({'Amino acid' : k.keys(), i : k.values() })
        else:
            fi2 = pd.DataFrame({'Amino acid' : k.keys(), i : k.values() })
            fi = pd.merge(fi, fi2, on = 'Amino acid')

    if tf == 'table':
        dfi.export(fi, path_out)
    if tf == 'figure':
        fi.plot(x = 'Amino acid', y = list_of_names, kind = 'bar', figsize = (15,8))
        plt.title("Percentage content of all amino acids ", size = 26)
        plt.xticks(rotation = 0, ha='right', size = 12, fontweight = 'bold')
        plt.xlabel('Amino acid')
        plt.ylabel('Percentage content')
        plt.savefig(path_out, bbox_inches='tight', facecolor = 'white')


# In[22]:


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog = 'POB', formatter_class=argparse.RawDescriptionHelpFormatter, description='''POB - Phylogeny of birds
----------------------------------------------------------------------------
Available organisms names you can check in file avaliable_organism.                         
If you do not have that file please run first only --update option.
After run you can find file: actual_analysis_organisms_names, that contains correct names of organism that you had used. Now you can use it now as input file with option 
--no_download to avoid re-downloading.''')

    parser.add_argument('-in_list', type=str,
                        help='Input as string e.g. "cat,dog"')

    parser.add_argument('-in_file', type=str,
                        help='Input as file - give path to file with one name in each row e.g. cat<new_line>dog')

    parser.add_argument('-median', type=str,
                        help='Path to output - median length barplot.')

    parser.add_argument('-aver_fig', type=str,
                        help='Path to output - average lenth barplot.')
    
    parser.add_argument('-aver_table', type=str,
                        help='Path to output - average length table')

    parser.add_argument('-aa_content_table',type=str,
                        help='Path to output - amino acid content table')

    parser.add_argument('-aa_content_figure',type=str,
                        help='Path to output - amino acid content plot')
    
    parser.add_argument('-merge',type=str,
                        help='Merge sequence files to clustering')
    
    parser.add_argument('-clusters_rep',type=str,
                        help='Divide clusters into separate files with paralogs')
    
    parser.add_argument('-clusters_no_rep',type=str,
                        help='Divide clusters into separate files without paralogs')
    
    parser.add_argument('-bootstrap',type=str,
                        help='Filter trees after bootstraping with given threshold')
    
    parser.add_argument('-mc',type=str,
                        help='Makes MC tree with given majority percent value')
    
    parser.add_argument('-auto',type=str,
                        help='Makes MC tree with repetitions automatically for inexpirienced users.')

    parser.add_argument('--update', action = 'store_true',
                        help='Check and update available proteoms file. / Create available_organims and id_protom_id files.')

    parser.add_argument('--no_download', action = 'store_true',
                        help='Do not download data (WARNING! in that case you have to be sure that all names are corrected. Use actual_analysis_organisms_names as input.')
    
    
    args = parser.parse_args()

    if args.update:
        download_id_proteoms_by_taxa()
        sys.exit()

    exist = os.path.exists(here + '/proteom_id' )
    if not exist:
        download_id_proteoms_by_taxa()
            
  #  if args.in_list:
  #      lists = args.in_list.split(',')
  #  elif args.in_file:
  #      lists = make_list_from_file(args.in_file)
  #  else:
  #      print('You have to give names of organisms as list or file.')
  #      print('Check help.')
  #      sys.exit()
  #  #print(lists)
  #  downloaded = download_from_list(lists)



    if args.no_download:
        if args.in_list:
            downloaded = args.in_list.split(',')
        elif args.in_file:
            downloaded = make_list_from_file(args.in_file)
        else:
            print('You have to give names of organisms as list or file.')
            print('Check help.')
            sys.exit()
                               
    else:
        if args.in_list:
            lists = args.in_list.split(',')
        elif args.in_file:
            lists = make_list_from_file(args.in_file)
        else:
            print('You have to give names of organisms as list or file.')
            print('Check help.')
            sys.exit()

        downloaded = download_from_list(lists)
    if args.median:
        make_median_plot(downloaded, args.median)
    if args.aver_fig:
         make_average_plot(downloaded, args.aver_fig)
    if args.aa_content_table:
        amino_content_from_list(downloaded, 'table', args.aa_content_table)
    if args.aa_content_figure:
        amino_content_from_list(downloaded, 'figure', args.aa_content_figure)
    if args.aver_table:
         make_average_table(downloaded, args.aver_table)
    if args.merge:
        merger(downloaded, args.merge)
    if args.clusters_rep:
        clusters_with_rep(args.clusters_rep)
    if args.clusters_no_rep:
        clusters_without_rep(args.clusters_no_rep)
    if args.bootstrap:
        bootstrapping(args.bootstrap)
    if args.mc:
        mc(args.mc)
    if args.auto:
        auto_rep(downloaded,args.auto)

