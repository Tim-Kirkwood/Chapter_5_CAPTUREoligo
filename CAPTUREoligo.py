# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 18:23:57 2021

@author: u03132tk
"""
from snapgene_reader import snapgene_file_to_seqrecord
import sys
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import primer3 as p3
import pandas as pd
import plotly.express as px
def write_histogram (x_list,name_list, colour_list, file):
        fig = go.Figure()
        #count = 1
        fig = make_subplots(rows=2, cols=1, row_heights=[0.9, 0.1], shared_xaxes=True)
        for x_array, x_name, x_color in zip(x_list, name_list, colour_list):
            fig.add_trace(go.Histogram(x=x_array, name = x_name, marker = {'color':x_color}), row = 1, col = 1)
            fig.add_trace(go.Box(x=x_array, 
                                        marker_symbol='line-ns-open', 
                                        marker_color=x_color,
                                        boxpoints='all',
                                        jitter=0,
                                        fillcolor='rgba(255,255,255,0)',
                                        line_color='rgba(255,255,255,0)',
                                        hoveron='points',
                                        name=f'{x_name}_rug'),
                          row = 2, 
                          col = 1)
            #count+=1
        
        # Overlay both histograms
        fig.update_layout(barmode='overlay')
        # Reduce opacity to see both histograms
        fig.update_traces(opacity=0.75)
        fig.write_html(file)

def make_dir_if_absent(dir_path):
    #if you want to extend to other directories add option to make original_cwd another directory rather than cwd with os.chdir, but make sure you set everythin back to the original cwd regardless
    original_cwd=os.getcwd()
    #TODO make this more flexible 
    if original_cwd not in dir_path:
        sys.exit(f'Current working directory is:\n{original_cwd}.\n\nThis is not part of the directory being built: \n{dir_path}\n\nThis will cause naming issues - use os.chdir to change cwd if you cannot use this function within the current working directory.')
    if dir_path.index(original_cwd) != 0:
        sys.exit('dir_path indexing error')
    
    original_dirs = [item for item in original_cwd.split('\\') if len(item)>0]
    #0 indexed but len is 1-indexed so >=.  from first tests we know original dirs matches start of the new dirs, but we dont want to include the duplicates.  just testing for presence\absence in old_dirs means you cant have a dir with same name subdir.
    new_dirs=[item for index, item in enumerate(dir_path.split('\\')) if len(item)>0 and index >= len(original_dirs)]
    #print (f'original_dirs: {original_dirs}\nnew_dirs: {new_dirs}\ndir_path: {dir_path}')
    for folder in new_dirs:
        path=fr'{os.getcwd()}\{folder}'
        #print (f'new_path:  {path}\n')
        try:
            os.mkdir(path)
        except FileExistsError:
            None#print (f'exists already:  {path}\n')
        os.chdir(path)
    os.chdir(original_cwd)
    return dir_path

def write_diagrams(dictionary_list, filepath):
    parent_dir = make_dir_if_absent(filepath)
    dictionary_mapping = dictionary_list
    

    
    write_histogram([[p['homology_g_stretch'] for p in dictionary_mapping]], 
                    ['homology_g_stretch'],
                    ['red'],
                    fr'{parent_dir}\homology_g_stretch.html')

    write_histogram([[p['spacer_gc'] for p in dictionary_mapping]], 
                    ['spacer_gc'],
                    ['red'],
                    fr'{parent_dir}\spacer_gc.html')
    write_histogram([[p['homology_gc'] for p in dictionary_mapping]], 
                    ['homology_gc'],
                    ['red'],
                    fr'{parent_dir}\homology_gc.html')
    write_histogram([[p['general_score'] for p in dictionary_mapping]], 
                    ['general_score'],
                    ['red'],
                    fr'{parent_dir}\general_score.html')
    write_histogram([[p['target_score'] for p in dictionary_mapping]], 
                    ['target_score'],
                    ['red'],
                    fr'{parent_dir}\target_score.html')
    write_histogram([[p['homology_hairpin'] for p in dictionary_mapping]], 
                    ['homology_hairpin'],
                    ['red'],
                    fr'{parent_dir}\homology_hairpin.html')
    
    
        
        
        
#TODO now write target protospacers to a fasta file, use in align twoi sequences vs dem, then export all to excel.  
#TODO ONLY WRITE SPACERS AND HOMO REGIONS WITH OK GC AND STRETCH PARAMS
#sort by alignet lenght, and you should see all of the length 8 one thes and then discount them hopefully,
#TODO you need to check for specifity by removing bases from one end, which is what you are doing now.  You then need to repeat from the other end as well, up to the score you get originally, to see if theres a longer one there.   maybe not necessery if you are swithcing 5\3 prime etc?
def check_protospacer(protospacer : Seq, spacer_from_top_strand : bool, spacer_from_bottom_strand : bool, 
                      top_dna_sequence_5to3 : Seq, bottom_dna_sequence_5to3 : Seq, allow_absent_protospacer : bool) -> bool:
    '''
    

    Parameters
    ----------
    protospacer : Seq
        DESCRIPTION.
    spacer_from_top_strand : bool
        DESCRIPTION.
    spacer_from_bottom_strand : bool
        DESCRIPTION.
    top_dna_sequence_5to3 : Seq
        DESCRIPTION.
    bottom_dna_sequence_5to3 : Seq
        DESCRIPTION.
    allow_absent_protospacer : bool
        DESCRIPTION.

    Returns
    -------
    bool
        Is the protospacer unique, found in the correct strand, and (optionally) found at all.

    '''
    reverse_protospacer_present =  protospacer[::-1] in top_dna_sequence_5to3 or protospacer[::-1] in bottom_dna_sequence_5to3
    if reverse_protospacer_present:
        return False
    if spacer_from_top_strand:
        unique_protospacer = top_dna_sequence_5to3.count(protospacer) == 1 and bottom_dna_sequence_5to3.count(protospacer) == 0
    elif spacer_from_bottom_strand:
        unique_protospacer = top_dna_sequence_5to3.count(protospacer) == 0 and bottom_dna_sequence_5to3.count(protospacer) == 1
    no_protospacer = top_dna_sequence_5to3.count(protospacer) == 0 and bottom_dna_sequence_5to3.count(protospacer) == 0
    if no_protospacer and not allow_absent_protospacer:
        sys.exit(f'{protospacer} is not present in any sequence orientations')
    return unique_protospacer or no_protospacer
 
def score_protospacer(protospacer : Seq, 
                      spacer_from_top_strand : bool, spacer_from_bottom_strand : bool, 
                      top_dna_sequence_5to3 : Seq, bottom_dna_sequence_5to3 : Seq, 
                      allow_absent_protospacer : bool) -> int:
    '''
    

    Parameters
    ----------
    protospacer : Seq
        DESCRIPTION.
    spacer_from_top_strand : bool
        DESCRIPTION.
    spacer_from_bottom_strand : bool
        DESCRIPTION.
    top_dna_sequence_5to3 : Seq
        DESCRIPTION.
    bottom_dna_sequence_5to3 : Seq
        DESCRIPTION.
    allow_absent_protospacer : bool
        DESCRIPTION.

    Returns
    -------
    max_score : int
        Maximum sized substring of protspacer that is non-unique (checks 5' -> 3' and vice-versa).

    '''
    max_score = 0 
    no_higher_scores = False
    for start_i in range(0,len(protospacer)):
        for end_i in range(len(protospacer),5, -1):
            if start_i < end_i and max_score<end_i - start_i:
                
                partial_protospacer = protospacer[start_i : end_i]
               
                ok_protospacer = check_protospacer(partial_protospacer, spacer_from_top_strand, 
                                                   spacer_from_bottom_strand, top_dna_sequence_5to3, 
                                                   bottom_dna_sequence_5to3, allow_absent_protospacer)
                if ok_protospacer:
                    continue
                #you have a hit for the sub-spacer from start_i to end_i
                score = end_i - start_i
                if score>max_score:
                    max_score = score
                #move to next start_i
                break
            else:
                no_higher_scores = True
                break
        if no_higher_scores:
            break

    return max_score
                
def define_source (spacer_from_top_strand : bool, spacer_from_bottom_strand : bool, 
                   top_5to3 : Seq, bottom_5to3 : Seq) -> str:
    '''
    

    Parameters
    ----------
    spacer_from_top_strand : bool
        DESCRIPTION.
    spacer_from_bottom_strand : bool
        DESCRIPTION.
    top_5to3 : Seq
        DESCRIPTION.
    bottom_5to3 : Seq
        DESCRIPTION.

    Returns
    -------
    str
        String from which to take the protospacer.

    '''        
    
    if spacer_from_top_strand:
        #print ('the spacer is from the top strand')
        var_out = top_5to3.lower()
        #print (f'top_5to3: {type(top_5to3)}\nbottom_5to3: {type(bottom_5to3)}\nvar_out: {type(var_out)}\n')
        return var_out
    elif spacer_from_bottom_strand:
        #print ('the spacer is from the bottom strand')
        var_out = bottom_5to3.lower()
        #print (f'top_5to3: {type(top_5to3)}\nbottom_5to3: {type(bottom_5to3)}\nvar_out: {type(var_out)}\n')
        return var_out
    #return source_seq

def process_locus (locus : int, 
                   #its a number referring to edge of capture region i think
                   #top 5to 3 is being referenced from outside of function
                   search_left : bool, search_right : bool, 
                   spacer_from_top_strand : bool, spacer_from_bottom_strand : bool, 
                   search_span : int) -> range:
    '''
    

    Parameters
    ----------
    locus : int
        DESCRIPTION.
    search_left : bool
        DESCRIPTION.
    search_right : bool
        DESCRIPTION.
    spacer_from_top_strand : bool
        DESCRIPTION.
    spacer_from_bottom_strand : bool
        DESCRIPTION.
    search_span : int
        DESCRIPTION.

    Returns
    -------
    range
        range obj containing indexes within arbritray string to search for protospacers 
        based on input strand, locus, BGC insert flank, size to search in.

    '''
    if search_left:
        if spacer_from_top_strand:
            
            #source_seq = top_5to3
            processed_locus = locus
            #print ('search left, spacer_from_top_strand, the locus is unchanged\n')
            search_range = range(processed_locus-search_span, processed_locus)
        elif spacer_from_bottom_strand:
            #source_seq = bottom_5to3#
            #print ('search left, spacer_from_bottom_strand, the locus is changed\n')
            processed_locus = len(top_5to3) - locus
            search_range = range(processed_locus, processed_locus+search_span)
    elif search_right:
        if spacer_from_top_strand:
            #source_seq = top_5to3
            #print ('search right, spacer_from_top_strand, the locus is unchanged\n')
            processed_locus = locus
            search_range = range(processed_locus, processed_locus+search_span)
        elif spacer_from_bottom_strand:
            #source_seq = bottom_5to3
            #print ('search right, spacer_from_bottom_strand, the locus is changed\n')
            processed_locus = len(top_5to3) - locus
            search_range = range(processed_locus - search_span, processed_locus)
    return search_range

def find_sites (spacer_from_top_strand:bool, spacer_from_bottom_strand:bool, pam_source_index:int, pam_sequence:str, top_5to3:str):
    if spacer_from_top_strand:
        
        #pam_start_index = len(top_5to3) - pam_source_index
        pam_end_index = pam_source_index + len(pam_sequence)
        pam_strand_cut_site = pam_end_index + 18.5
        non_pam_strand_cut_site = pam_end_index + 23.5
        
        return {'pam_strand_cut_site' : pam_strand_cut_site,
                'non_pam_strand_cut_site' : non_pam_strand_cut_site,
                'pam_start_index' : pam_source_index,
                'pam_end_index' : pam_end_index,
                'homologous_sequence_start' : pam_source_index - 12,
                'homologous_sequence_end' : pam_source_index + len(str(pam_sequence)) + 24}
    elif spacer_from_bottom_strand:
        pam_start_index = len(top_5to3) + 1 - (pam_source_index + len(pam_sequence))
        pam_end_index = len(top_5to3) - pam_source_index 
        pam_strand_cut_site = pam_start_index - 18.5
        non_pam_strand_cut_site = pam_start_index - 24.5 
        return {'pam_strand_cut_site' : pam_strand_cut_site,
                'non_pam_strand_cut_site' : non_pam_strand_cut_site,
                'pam_start_index' : pam_start_index,
                'pam_end_index' : pam_end_index,
                'homologous_sequence_start' : pam_start_index - 24 - 1,
                'homologous_sequence_end' : pam_end_index + 12}

def base_stretch(string, base):
    stretch=0
    max_stretch=0
    string_l = string.lower()
    base_l = base.lower()
    for bp in string_l:
        if bp == base_l:
            stretch+=1
        else:
            if stretch > max_stretch:
                max_stretch = stretch
            stretch=0
    #if last base is part of biggest stretch it wont be added to max_stretch
    if stretch > max_stretch:
        return stretch
    else:
        return max_stretch

def flexible_operator(item, operator, comparitor): 
    if operator == '>':
        return item > comparitor
    elif operator == '>=':
        return item >= comparitor
    if operator == '<':
        return item < comparitor
    elif operator == '<=':
        return item <= comparitor
    if operator == '==':
        return item == comparitor
    elif operator == '!=':
        return item != comparitor

def filter_dict_obj(dictionary, threshold_map, operator_list):#tests for threshold exceeding - >0 = failed
    #just define the key mapping dont be lazy                
    counts = []
    target_input_keys = {}
    for threshold_key in threshold_map.keys():
        for input_key in dictionary.keys():
            if input_key in threshold_key:
                target_input_keys[input_key] = threshold_key
                break
    for (input_key, map_key), operator in zip(target_input_keys.items(), operator_list):
        #find threshold_map_val corresponding to dictionary key
        
        #map_key = target_input_keys[input_key]
        #print (input_key, map_key)
        if flexible_operator(dictionary[input_key], operator, threshold_map[map_key]):
            counts += [[input_key, 1]]
    return counts

def compute_protospacer (locus, #locus is not in here also make it automaticaaly go to 0 if exceed contig
                         top_5to3, bottom_5to3, 
                         target_top_5to3, target_bottom_5to3, 
                         pam_constant, pam_variable, 
                         protospacer_len, 
                         spacer_from_top_strand, spacer_from_bottom_strand, 
                         search_left, search_right, search_span, 
                         constant_primer_region_5to3, 
                         terminus, 
                         max_general_score,
                         max_homology_g_stretch,
                         max_homology_gc,
                         min_homology_hairpin,
                         max_spacer_gc, 
                         max_target_score,
                         check_spacer_for_any_pam = False, 
                         check_spacer_for_defined_pam = False,
                         check_protospacer_for_defined_terminal_pam = False,
                         check_protospacer_for_any_terminal_pam = False,
                         pam_db = ['TTA', 'TTC', 'TTG', 'CTA', 'CTC', 'CTG']): #YTN PAM https:\\www.nature.com\articles\nature17945
    protospacers = []
    failed_protospacers = []
    #pick top or bottom strand to extract protospacers from 
    source_seq = define_source (spacer_from_top_strand, spacer_from_bottom_strand, top_5to3, bottom_5to3)
    #print (f'source_seq is {len(top_5to3)} bp\nsource seq is same as top strand:  {str(source_seq).lower() == str(top_5to3).lower()}\nsource seq is same as bottom strand:  {str(source_seq).lower() == str(bottom_5to3).lower()}\n')
    #get range of indexes within source_seq to extract protospacers from
    search_range = process_locus (locus, search_left, search_right, spacer_from_top_strand, spacer_from_bottom_strand, search_span)
    #print ('search range ok')
    #check each potent-ial PAM
    
    #TODO - have a box saying 'general PAM in spacer', 'specific pam in spacer' and make this a filterable option
    for pam_source_index in search_range:
        correct_constant_region = source_seq[pam_source_index:pam_source_index+len(pam_constant)] == pam_constant
        correct_variable_region = source_seq[pam_source_index + len(pam_constant)] in pam_variable
        correct_pam = correct_constant_region and correct_variable_region
        if correct_pam:
            #print ('correct pam')
            #check protospacer
            
            len_pam = len(pam_constant) + 1 
            protospacer = source_seq[pam_source_index + len_pam 
                                     : pam_source_index + len_pam + protospacer_len]
            protospacer_bases_are_canonincal = all([base.lower() in 'atcg' for base in protospacer])
            if not protospacer_bases_are_canonincal:#dont have Ns
                continue
            if check_spacer_for_any_pam:
                pam_in_fw_spacer = any([pam_sequence in protospacer for pam_sequence in pam_db])
                #pam_in_rev_spacer = any([pam_sequence in protospacer.reverse_complement() for pam_sequence in pam_db])
                if pam_in_fw_spacer:# or pam_in_rev_spacer:
                    continue
            
            pam_sequence = source_seq[pam_source_index : pam_source_index + 1 + len(pam_constant)]
            
            if check_spacer_for_defined_pam:
                if pam_sequence in protospacer:# or pam_sequence in protospacer.reverse_complement():
                    continue
            if check_protospacer_for_any_terminal_pam:
                any_terminal_pam_present = False
                for db_pam in pam_db:
                    protospacer_ends = [protospacer[0 : len(db_pam)], 
                                        protospacer[-len(db_pam) : ]]
                    if db_pam in protospacer_ends:
                        any_terminal_pam_present = True
                        break
                if any_terminal_pam_present:
                    continue
            if check_protospacer_for_defined_terminal_pam:
                protospacer_ends = [protospacer[0 : len(pam_sequence)], 
                                    protospacer[-len(pam_sequence) : ]]
                if pam_sequence in protospacer_ends:
                    continue
            #print (f'protospacer is {protospacer}')
            ok_protospacer = check_protospacer(protospacer, spacer_from_top_strand, spacer_from_bottom_strand, 
                                               top_5to3, bottom_5to3, allow_absent_protospacer = False)
            #print ('protospacer is ok')
            if ok_protospacer:
                
                #get scores for various quality metrics
                
                #uniqueness score for whole genome
                general_score = score_protospacer(protospacer, spacer_from_top_strand, spacer_from_bottom_strand, top_5to3, bottom_5to3, allow_absent_protospacer = False)
                
                #uniqueness score within CAPTURE fragment - should be as low as possible to prevent internal cleavage
                target_score = score_protospacer(protospacer, spacer_from_top_strand, spacer_from_bottom_strand, target_top_5to3, target_bottom_5to3, allow_absent_protospacer = True)
                
                #get important sites
                #source https:\\eu.idtdna.com\pages\technology\crispr\crispr-genome-editing\Alt-R-systems\crispr-cas12a
                
                sites = find_sites(spacer_from_top_strand, spacer_from_bottom_strand, pam_source_index, pam_sequence, top_5to3)
                
                # check overhang - all being done relative to top strand not source seq, to facilate manual checking by user in snapgene
                homologous_region = [sites['homologous_sequence_start'], sites['homologous_sequence_end']]
                homologous_seq_top_5to3 = top_5to3[homologous_region[0]:homologous_region[1]]
                homologous_seq_bottom_5to3 = homologous_seq_top_5to3.reverse_complement()
                
                #gc_homologous_region = GC(homologous_seq_top_5to3)
                
                #if homology_g_stretch<4 and 20<gc_homologous_region<72:
                if terminus == 'start':
                    primer = homologous_seq_bottom_5to3 + constant_primer_region_5to3
                    #homo_seq = homologous_seq_bottom_5to3
                elif terminus == 'end':
                    primer = homologous_seq_top_5to3 + constant_primer_region_5to3
                    #homo_seq = homologous_seq_top_5to3
                    
                #TODO this isnt working - it misses Gs that are on opposite strand in some cases.
                homology_g_stretch = max([base_stretch(homologous_seq_bottom_5to3, 'G'),
                                          base_stretch(homologous_seq_top_5to3, 'G')])
                #format:  forward gRNA should always be 'constant region' then 'protospacer on strand pointing out of the insert 5 -> 3'
                gRNA_f = 'AATTAATACGACTCACTATAGGGAATTTCTACTGTTGTAGAT' + str(protospacer)
                gRNA_r = str(Seq(gRNA_f).reverse_complement())
                protospacer_data = {'pam_start_index' : sites['pam_start_index'], 
                                    'pam_end_index' : sites['pam_end_index'],
                                    'pam_strand_cut_index (+- 0.5)' : sites['pam_strand_cut_site'],
                                    'non_pam_strand_cut_index': sites['non_pam_strand_cut_site'],
                                    'homologous_region_top_5to3' : homologous_region,
                                    'homologous_seq_top_5to3' : str(homologous_seq_top_5to3),
                                    'homologous_seq_bottom_5to3' : 	str(homologous_seq_bottom_5to3),
                                    'homology_g_stretch' : homology_g_stretch,
                                    'spacer_gc' : GC(protospacer),
                                    'homology_gc' : GC(homologous_seq_top_5to3),#same top and bottom
                                    'homology_hairpin' : max([p3.calcHairpin(str(homologous_seq_top_5to3), 
                                                                             mv_conc, 
                                                                             dv_conc, 
                                                                             dntp_conc).dg,
                                                              p3.calcHairpin(str(homologous_seq_bottom_5to3), 
                                                                             mv_conc, 
                                                                             dv_conc, 
                                                                             dntp_conc).dg]),
                                    'pam' : str(pam_sequence), 
                                    'protospacer' : str(protospacer), 
                                    'general_score' : general_score, 
                                    'target_score' : target_score,
                                    'primer':str(primer),
                                    'gRNA_f' : gRNA_f,
                                    'gRNA_r' : gRNA_r}
                
                threshold_map = {'max_general_score' : max_general_score,
                                 'max_homology_g_stretch' : max_homology_g_stretch,
                                 'max_homology_gc' : max_homology_gc,
                                 'min_homology_hairpin' : min_homology_hairpin,
                                 'max_spacer_gc' : max_spacer_gc, 
                                 'max_target_score' : max_target_score}
                filter_protospacer = filter_dict_obj(protospacer_data, threshold_map, ['>','>','>','<','>','>'])
                if filter_protospacer == []:   
                    protospacers.append(protospacer_data)
                else:
                    failed_protospacers += filter_protospacer
    return {'protospacers' : protospacers,
            'failed_protospacers' : failed_protospacers}



    
def concatenate_records (filepath_gbk, gbk_directory, out_directory = None):
    gbk_directory = gbk_directory.replace("/","\\")
    records = []
    for index, file in enumerate(os.listdir(gbk_directory)):
        with open(fr'{gbk_directory}\{file}', 'r') as file:
            seqrecords = list(SeqIO.parse(file, 'genbank'))
            if len(seqrecords) == 0:
                sys.exit('error reading file - no seqrecords')
            records += seqrecords
            #concatenate contigs
    concatenated_record = records[0]#+Seq('z')
    junctions = []
    for r in records[1:]:
        #concatenated_record += Seq('z')
        concatenated_record += r
        junctions += [len(concatenated_record.seq)]
    #seqrecord.seq.alphabet = generic_dna
    concatenated_record.seq = concatenated_record.seq.lower()
    if out_directory == None:
        with open(rf'{gbk_directory}\{filepath_gbk}'.replace('/','\\'), "w") as output_handle:
            SeqIO.write(concatenated_record, output_handle, "genbank")
    else:
        with open(rf'{out_directory}\{filepath_gbk}'.replace('/','\\'), "w") as output_handle:
            SeqIO.write(concatenated_record, output_handle, "genbank")
    return {'concatenated_record' : concatenated_record, 
            'junctions' : junctions}

# =============================================================================
# READ FROM HERE
# =============================================================================

#looks like you need a pam to activate trans cleavage and that protospacer mismatches are very poorly tolerated 
#- so dont worry too much if there is off target sites that are not identical
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6938695/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6628903/    
 
#put path to your input gbk file
results_dir = 'hyalina_31_3_23'
accepted_suffixes = ['.gbk', 'gb']
dir_with_input_gbks = r'C:/Users/u03132tk/.spyder-py3/hyalina_assembly_bif/ncbi-genomes-2023-03-31/GCF_003967355/gbks'
concatenated_gbk_filename = 'NZ_BIFH01000000.gbk'

#which receiver is making up your backbone with pbe45 (which is compulsory and hardcoded in)
forward_receiver = 'pbe48'#'pbe48' #or 'pbe44'.  

#probably dont want to change - this says PAM is always TT plus A, C or G
pam_constant = 'tt'
pam_variable = 'acg'

protospacer_len = 18

#up\downstream regions
upstream_search_span = 0#2500
downstream_search_span = 5000
print ()
#crude specificity check: (might integrate with BLAST, if i get time - could have issues with small seqs)
#takes letters 1 to N for protospacer of length N.  
#Checks if it has >1 copy in genome in either orientation in both strands.  
#It then repeats for bases 1 to N-1, N-2...until length of protospacer = 5.  Score is length of longest off-target hit.  
#general score is whole genome, target score is score within BGC fragment (highlighted as more important in the paper) 
max_general_score = 20
max_homology_g_stretch = 3
max_homology_gc = 72
min_homology_hairpin = -9000
max_spacer_gc = 72
max_target_score = 20

#do all bgcs for a genome on lax settings, then run individually based on settings 
#put one pair of orfs per BGC to be captured in the genome you put in 
codon_boundary_list = [#HYALINA
                       #['EMA17_RS43455', 'EMA17_RS43615']
                       #["EMA17_RS16455", "EMA17_RS16845"]
                       ["EMA17_RS34675", "EMA17_RS35025"],
                       #["EMA17_RS43435", "EMA17_RS43605"]
                       #["EMA17_RS44205", "EMA17_RS44595"]

                       #AUSTRALIENSIS
                       #["C8E97_RS15620", "C8E97_RS15970"],
                       #["C8E97_RS16320", "C8E97_RS16720"],
                       #["C8E97_RS18900", "C8E97_RS19285"],
                       #["C8E97_RS19665", "C8E97_RS35415"],
                       #["C8E97_RS22525", "C8E97_RS22865"]
                       
                       
                       #issy
                       #['ctg7_174','ctg7_185']
                       #]
                       
                       #ek
                       #['SSIG_RS06550', 'SSIG_RS06370']
                       ]
                       

#can change according to your polymerase - these are for PrimeStar
mv_conc = 0
dv_conc = 2 
dntp_conc = 0.4



#STOP READING NOW
results_filepath = rf'{os.getcwd()}\{results_dir}'
make_dir_if_absent(results_filepath)
assert concatenated_gbk_filename[concatenated_gbk_filename.rindex('.'):] in accepted_suffixes
print (f'concatenating records in {dir_with_input_gbks}...')
concatenated_data = concatenate_records(concatenated_gbk_filename, dir_with_input_gbks, results_filepath)

seqrecord = concatenated_data['concatenated_record']

list_of_codons = []
for start_codon, end_codon in codon_boundary_list:
    print (f'\nprocessing codons - {start_codon} extended upstream by {upstream_search_span}bp, {end_codon} extended downstream by {downstream_search_span}bp')
    #checks for homology region:  each DNA receiver carried 39 bp of homology to the BGC fragment. 
    #This 39 bp homology included 24 bp sequence downstream of the PAM sequence (including the 18 bp of spacer), 
    #3 bp of the TTV (V = A, C, or G) PAM sequence, and 12 bp sequence upstream of the PAM sequence. 
    #The GC-content for the 39 bp sequence was kept in 20–72% range while avoiding 4 or more consecutive G nucleotides. 
    #Examples of ssDNA templates and primers used for receiver amplification can be found in Supplementary Data 1.
    
    #def find_boundary_indexes
    found_codon_pair = False
    
    start_count = 0
    end_count = 0
    for f in seqrecord.features:
        if f.type == 'CDS':
            if isinstance(f.qualifiers['locus_tag'], list):
                assert len(f.qualifiers['locus_tag']) == 1
                check_codon = f.qualifiers['locus_tag'][0]
            elif isinstance(f.qualifiers['locus_tag'], str):
                check_codon = f.qualifiers['locus_tag']
            else:
                sys.exit('oh no')
            if check_codon == start_codon:
                start_locus = f.location.start
                #start_strand = f.strand
                start_count += 1
                if start_count >1:
                    sys.exit(f'you have a record with duplicates of {start_codon}')
            elif check_codon == end_codon:
                end_locus = f.location.end
                #end_strand = f.strand
                end_count+=1
                if end_count >1:
                    sys.exit(f'you have a record with duplicates of {end_codon}')
    if [start_count, end_count].count(0) == 1:
        sys.exit('you have only one of the target codons')

        
    for junction_locus in concatenated_data['junctions']:
        if start_locus - upstream_search_span <= junction_locus <= end_locus + downstream_search_span:
            sys.exit('''your chosen insert (including up and downstream search spans) 
                     includes an inter-contig junction at position {junction_locus} in 
                     "{results_filepath}\{concatenated_gbk_filename}", which is not allowed''')
    #else:
    #    sys.exit(f'developer logic error: there are {start_count} starts and {end_count} ends')
        
    #this works when you get a single record off of ncbi - https:\\www.ncbi.nlm.nih.gov\nuccore\NZ_RBXO01000001.  I think issue arises from snapgene potentially during import export and file conversion.  Play around with this.
    #sys.exit('success')    
    
    if not [start_count, end_count] == [1, 1]:
        sys.exit(f'could not find {start_codon} or {end_codon} in any of the records provided')
    print (f'CODONS: {start_codon} ({start_locus}), {end_codon}({end_locus})\tLENGTH: {end_locus - start_locus}')
    list_of_codons.append({'start_codon' : start_codon,
                           'start_locus' : start_locus,
                           'end_codon' : end_codon,
                           'end_locus' : end_locus})  
for dictionary in list_of_codons:
    start_codon = dictionary['start_codon']
    start_locus = dictionary['start_locus']
    end_codon = dictionary['end_codon']
    end_locus = dictionary['end_locus']
    top_5to3 = seqrecord.seq.lower()#might not work for non-biopython parsers 
    bottom_5to3 = seqrecord.seq.reverse_complement().lower()
    target_top_5to3 = seqrecord.seq[start_locus:end_locus].lower()
    target_bottom_5to3 = target_top_5to3.reverse_complement().lower()

    assert all([i>=0 for i in [start_locus, 
                               upstream_search_span, 
                               end_locus, 
                               downstream_search_span]]), '''you cannot have negative start_locus, 
               upstream_search_span, end_locus or downstream_search_span- {start_locus, 
                                                                           upstream_search_span, 
                                                                           end_locus, 
                                                                           downstream_search_span}
               '''
    
    print ('----- FORWARD -----')
    spacer_from_top_strand = False
    spacer_from_bottom_strand = True
    search_left = True
    search_right = False
    locus = start_locus
    if start_locus < upstream_search_span:
        sys.exit('start locus extended beyond sequence len')
    if forward_receiver == 'pbe48':
        constant_primer_region_5to3 = 'ACAGCGACACACTTGCATCG' 
    elif forward_receiver == 'pbe44':
        constant_primer_region_5to3 = 'AACTTATATCGTATGGGGCTG'
    terminus = 'start'
                        
    all_fw_protospacers = compute_protospacer (locus, 
                                           top_5to3, 
                                           bottom_5to3, 
                                           target_top_5to3, 
                                           target_bottom_5to3, 
                                           pam_constant, 
                                           pam_variable, 
                                           protospacer_len, 
                                           spacer_from_top_strand,
                                           spacer_from_bottom_strand,
                                           search_left, 
                                           search_right,
                                           upstream_search_span,
                                           constant_primer_region_5to3, 
                                           terminus,
                                           max_general_score,
                                           max_homology_g_stretch,
                                           max_homology_gc,
                                           min_homology_hairpin,
                                           max_spacer_gc, 
                                           max_target_score)
    print (f'you have {len(all_fw_protospacers["protospacers"])}')
    fw_filepath = make_dir_if_absent(rf'{results_filepath}\forward_oligos\{start_codon}_to_{end_codon}')
    fw_protospacers = all_fw_protospacers['protospacers']
    write_diagrams(fw_protospacers, 
                   fw_filepath)
    df = pd.DataFrame(all_fw_protospacers['failed_protospacers'], 
                      columns=['filter_param', 'count'])
    fig = px.pie(df, 
                 values='count', 
                 names='filter_param', 
                 title='Params that cause greatest loss of CAPTURE oligos')
    fig.write_html(fr'{fw_filepath}\filter_pie.html')
    
    with open (fr'{fw_filepath}\fw_protospacers_{start_codon}_to_{end_codon}_{upstream_search_span}.txt', 'w') as fw_file:
        for index, protospacer in enumerate(fw_protospacers):
            fw_file.write(f'\n\nProtospacer {index}:\n')
            for key, val in protospacer.items():
                fw_file.write(f'{key}\t{val}\n')
    df = pd.DataFrame(all_fw_protospacers['protospacers'])
    df.to_csv(rf"{fw_filepath}\{start_codon}_fw.csv", index=False)
    print ('----- REVERSE -----')
    spacer_from_top_strand = True
    spacer_from_bottom_strand = False
    search_left = False
    search_right = True
    locus = end_locus
    #does it extend past concatenated sequwnce
    if end_locus + downstream_search_span > len(seqrecord):
        sys.exit('end locus extends beyond sequence len')

    constant_primer_region_5to3 = 'GACGCTCAGTGGAACGAAAAC' #45
    terminus = 'end'
    all_rev_protospacers = compute_protospacer(locus, 
                                           top_5to3, 
                                           bottom_5to3, 
                                           target_top_5to3, 
                                           target_bottom_5to3, 
                                           pam_constant, 
                                           pam_variable, 
                                           protospacer_len, 
                                           spacer_from_top_strand,
                                           spacer_from_bottom_strand,
                                           search_left, 
                                           search_right,
                                           downstream_search_span,
                                           constant_primer_region_5to3, 
                                           terminus,
                                           max_general_score,
                                           max_homology_g_stretch,
                                           max_homology_gc,
                                           min_homology_hairpin,
                                           max_spacer_gc, 
                                           max_target_score)
    print (f'you have {len(all_rev_protospacers["protospacers"])}')
    rev_filepath = make_dir_if_absent(rf'{results_filepath}\reverse_oligos\{start_codon}_to_{end_codon}')
    rev_protospacers = all_rev_protospacers['protospacers']
    write_diagrams(rev_protospacers, 
                   rev_filepath)
    

    df = pd.DataFrame(all_rev_protospacers['failed_protospacers'], 
                      columns=['filter_param', 'count'])
    fig = px.pie(df, 
                 values='count', 
                 names='filter_param', 
                 title='Params that cause greatest loss of CAPTURE oligos')
    fig.write_html(fr'{rev_filepath}\filter_pie.html')
    with open (fr'{rev_filepath}\rev_protospacers_{start_codon}_to_{end_codon}_{downstream_search_span}.txt', 'w') as rev_file:
        for index, protospacer in enumerate(rev_protospacers):
            rev_file.write(f'\n\nProtospacer {index}:\n')
            for key, val in protospacer.items():
                rev_file.write(f'{key}\t{val}\n')
    
    df = pd.DataFrame(all_rev_protospacers['protospacers'])
    df.to_csv(rf"{rev_filepath}\{end_codon}_rev.csv", index=False)