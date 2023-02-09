
# -*- coding: utf-8 -*-
"""
FUNCTIONS TO ALIGN MAUS WORDS TO ORIGINAL ELAN WORDS

- read in ELAN file
    - should have already imported MAUS TextGrid into ELAN
- extract original word tier, MAUS orthographic word tier (W_O_), and MAUS SAMPA word tier (W_S_)
    - create new SAMPA tier from original word tier
        - use 'g2p_sample_CV.xlsx' and parse_orthography() function
    - perform cleaning
- get matrix of NW distance between every pair of original and MAUS SAMPA words
    - create substitution dictionary: key = (MAUS word, original word), value = distance
- do NW alignment of list of MAUS and original SAMPA words
    - use substitution dictionary (requires pairwise2.align.globalds/globaldx/etc.)
- write aligned MAUS and original sequences to CSV for manual checking
"""

import numpy as np
import os,re,pandas
from Bio import pairwise2
from string import punctuation
    # Global variable
punct = ''.join([x for x in punctuation if x != "'" and x != "`" and x!= ':'])
punct = punct + ' '
weird,d_weird = False,{}

    # Saddest function ever
def fillWDict(I):
    """Fills the dictionary 'd_weird' for WEIRD languages."""
    path = os.path.join(I.d['in_dir'],I.d['eaf'].name+".csv")
    with open(path,'r',encoding="utf-8") as f:
        for line in f:
            wd,ws = line.split(";")
            wd = strip_text(wd,punct)
            if not wd in d_weird:
                d_weird[wd] = ws

    # Still not sure which functions call this but it cleans stuff
def strip_text(text, punct):
    """ may need to specify for language (e.g. glottal stop apostrophe in Arapaho)"""
    text = str(text)
    for char in punct:
        text = text.replace(char, '')
    return text.lower()
def strip_maus(text, g2p_dict={}, lang="", add_space=True):
    """ extract string content from MAUS tag (return original text if non-matching)"""
    match = re.search("<<[^>]+>([^>]+)>", text)
    if match:
        if add_space:
            return parse_orthography(match.group(1), g2p_dict.get(lang), sep_char=' ')[1]
        else:
            return match.group(1) # add option to separate by spaces if possible
    else:
        return text
    # Turns tier content to SAMPA
def parse_orthography(text, dff, V=False, sep_char = '-'):
    """ INPUT: text string, g2p df, list of punctuation to strip
        - loop through characters in text string, searching for longest possible
        character substrings first, and identifying it as C or V
        - also create list of phoneme substrings, stored as dash-separated string
        OUTPUT: dash-separated phoneme string and CV string
        
        TODO: deal with orthographic strings that correspond to two MAUS phones"""
    text = strip_text(text, punct)
    phonemes, maus_phones, cvs = [], [], []
    text_len = len(text)
    if weird:
        result = d_weird.get(text.strip(),"").strip()
        return "",result,""
    max_len = dff['text_len'].max()
    grouped = dff.groupby('text_len')
    i = 0; flag = False
    len_range = [x for x in range(max_len, 0, -1) if x in dff['text_len'].values] # in case missing lengths
    while i < len(text):
        flag = False
        for l in len_range:
            sub_str = text[i:i + l]
            if V: print(f"- searching for {sub_str} ({i})")
            if i + l > text_len:
                if V: print(f"--- {i + l} longer than {text_len + 1}")
                continue
            group = grouped.get_group(l)
            for ind, row in group.iterrows():
                text_char = row['text']
                if sub_str == text_char:
                    row_cvs = list(row['CV']); row_mauses = row['MAUS'].split(' ')
                    flag = (text_char, [], [])
                    for row_maus, row_cv in zip(row_mauses, row_cvs):
                        flag[1].append(row_maus); flag[2].append(row_cv)
                    break
            if flag != False:
                phonemes.append(flag[0])
                for row_maus, row_cv in zip(flag[1], flag[2]):
                    maus_phones.append(row_maus); cvs.append(row_cv)
                if V: print(f"-- found string: {flag[0]}, {flag[1]}; incrementing i to {i + l}")
                i += l
                break
        if flag == False:
            phonemes.append('%'); maus_phones.append('%'); cvs.append('X')
            #print(f"--- couldn't find character: {text[i:i+1]}, {i}, {text}, {text.encode('unicode_escape')}")
            i += 1
    phoneme_str = sep_char.join(phonemes); maus_str = sep_char.join(maus_phones); cv_str = ''.join(cvs)
    return phoneme_str, maus_str, cv_str
    # Core pairwise functions
def nw_matrix(matrix, maus, orig, kind='list', pol=-1):
    """ INPUT: zero matrix, maus and original sequences, list/str flag, pos/neg score arg"""
    new_mat = matrix.copy()
    for i, mi in enumerate(maus):
        for j, oi in enumerate(orig):
            if kind == 'list':
                score = pol * pairwise2.align.globalms(list(mi), list(oi), 2, -1, -1, -1, gap_char=['€'], score_only=True)
            elif kind == 'str':
                score = pol * pairwise2.align.globalms(mi, oi, 2, -1, -1, -1, score_only=True)
            new_mat[i][j] = score
    return new_mat
def create_subst_dict(x_str_set, y_str_set, matrix):
    """ create pairwise string substitution dict from pairwise distance matrix"""
    subst_dict = {}
    for i, x_str in enumerate(x_str_set):
        for j, y_str in enumerate(y_str_set):
            subst_dict[(x_str, y_str)] = matrix[i][j]
    return subst_dict
    # 'tier_to_list()' seems deprecated
    # (I removed 'get_maus_tiers')
def tier_to_list(tier):
    return [x.content for x in tier.elem]
    # Those are for cleaning tier content
def clean_str(str_):
    try:
        new_str = str_.translate(str.maketrans(dict.fromkeys(punct)))
        return new_str
    except AttributeError as e:
        print(str_)
        raise e
    return str(str_).translate(str.maketrans(dict.fromkeys(punct)))
def clean_maus(str_, g2p_dict={}, g2p=False):
    """ clean MAUS tags; optionally convert tags to SAMPA (g2p=[language name])"""
    if str_.startswith('<<') and str_.endswith('>'):
        new_str = str_.split('>')[-2]
        if g2p and new_str: # g2p it (don't apply on empty strings)
            new_str = parse_orthography(new_str, g2p_dict.get(g2p), sep_char=' ')[1] # maus_str
        return new_str
    else:
        return str_
def clean_orig_sampa(list_, lang, g2p_dict={}, trans_dict={},
                     clean_maus=True, morphs=False):
    """ clean punctuation and <<>> from original/maus words and convert
        original words to sampa"""
    # clean original words and convert to SAMPA
    cl_orig = [clean_str(x).lower() for x in list_[0]]
    cl_orig_samp = [parse_orthography(x, g2p_dict.get(lang), sep_char=' ')[1] for x in cl_orig]
    #cl_orig_samp = [x if not all(y in ['%', ' '] for y in x) else '-' for x in cl_orig_samp] # double '-' was causing errors
    # clean MAUS orthographic words
    cl_maus = [x.translate(trans_dict) for x in list_[1]]
    if clean_maus:
        cl_maus = [clean_maus(x) for x in cl_maus]
    # clean MAUS SAMPA words
    cl_maus_samp = [x.translate(trans_dict) for x in list_[2]]
    if clean_maus:
        cl_maus_samp = [clean_maus(x, g2p_dict,g2p=lang) for x in cl_maus_samp]
    ret_list = [cl_orig, cl_orig_samp, cl_maus, cl_maus_samp]
    # clean original morphemes and convert to SAMPA (optional)
    if morphs:
        cl_morph = [[clean_str(morph).lower() for morph in morphs] for morphs in list_[3]]
        cl_morph_samp = [[parse_orthography(x, g2p_dict.get(lang), sep_char=' ')[1] for x in morphs] for morphs in cl_morph]
        orig_morph = [[morph.replace('~', '=') for morph in morphs] for morphs in list_[3]] # new!
        ret_list += [cl_morph, cl_morph_samp, orig_morph] # orig_morph is new!
    return ret_list
    # Shouldn't that one be all the way to the top?
def get_g2p_dict(g2p,lang):
    """ read in g2p dictionary (and clean it up a little)"""
    import unicodedata
    g2p_dict = pandas.read_excel(g2p, sheet_name=None, engine='openpyxl')
    g2p_dict = {k:v.dropna(subset=['text']) for k, v in g2p_dict.items()}
    for k, v in g2p_dict.items():
        if not k == lang:
            continue
        v['MAUS'] = v['MAUS'].astype(str)
        if k in ['Hoocak', 'Northern Alta']: # stupid Excel won't let me put a "'" alone in a cell
            v.loc[len(v)] = [len(v), "'", '?', 'ʔ', 0, 'C']
        if k in ['Kakabe', 'Hoocak', 'Kamas']:
            v['text'] = v['text'].apply(lambda x: unicodedata.normalize('NFC', x))
        v['text_len'] = v['text'].apply(len)
    return g2p_dict
    # So that one is basically 'fillList', but we can't pick 'gap_char'/'add'
def realign_with_gaps(aligned_list, aligned_list2, unaligned_list, gap_char='€'):
    """ INPUT: list aligned by pairwise2 and unaligned list of the same length
        as the original aligned list
        OUTPUT: aligned second list, with gaps to match the aligned first list"""
    new_list = []; align_i = 0; new_i = 0
    la = len(aligned_list)
    
    while len(new_list) < la:
        #check = unaligned_list[new_i] if new_i < len(unaligned_list) else ""
        #print(align_i,aligned_list[align_i],aligned_list2[align_i],
        #      "|",new_i,check,la,len(unaligned_list))
        if aligned_list[align_i] == gap_char:
            if aligned_list2[align_i] == gap_char:
                new_list.append(unaligned_list[new_i])
            else:
                new_list.append(gap_char)
        else:
            new_list.append(unaligned_list[new_i])
            new_i += 1
        align_i += 1
    return new_list
def realign_with_gaps_OLD(aligned_list, unaligned_list, gap_char='€'):
    """ INPUT: list aligned by pairwise2 and unaligned list of the same length
        as the original aligned list
        OUTPUT: aligned second list, with gaps to match the aligned first list"""
    new_list = []; align_i = 0; new_i = 0
    while len(new_list) < len(aligned_list):
        if aligned_list[align_i] == gap_char:
            new_list.append(gap_char)
        else:
            new_list.append(unaligned_list[new_i])
            new_i += 1
        align_i += 1
    return new_list
    # Get word/morph groups
def get_word_morph_pairs(trans, w_tier, mb_tier, content=True):
    """ INPUT: Transcription object, word tier object, morph tier object
        OUTPUT: pairs of word segments and their child morph segments"""
    pair_list = []
    for wd_seg in w_tier:
        #print('get_word_morph_pairs():', wd_seg.content)
        d_child = wd_seg.childDict()
        if not mb_tier in d_child:
            pair_list.append((wd_seg.content,[])); continue
        mb_segs = d_child[mb_tier]
        if content:
            wd_seg = wd_seg.content; mb_segs = [x.content for x in mb_segs]
        pair_list.append((wd_seg, mb_segs))
    return pair_list
    # You really need a dict here
def lang_specific_cleaning(orig, morph, lang, trans_dict):
    """ INPUT: original word and morph lists, language, and translation dict
        - perform language-specific cleaning
        OUTPUT: cleaned word and morph lists"""
    orig = [x.translate(trans_dict) for x in orig]; morph = [[x.translate(trans_dict) for x in y] for y in morph]
    if lang in ['Bora']:
        orig = [x.replace(':', '') for x in orig] # about 20 of these in the whole corpus - could add V: to g2p?
    if lang in ['Mojeno']:
        orig = [x.replace('ʔ', "'") for x in orig]
    if lang in ['Beja']: # get rid of "//" and "/" words and morphs in Beja
        '''orig2, morph2 = [], []
        for o, m in zip(orig, morph):
            if any(x.isalpha() for x in o): # must contain some text characters?
                if not re.search('\w+_\d+', o): # exclude 'BI_572' type indices
                    mini_morph = [x for x in m if any(y.isalpha() for y in x)]
                    mini_morph = [x for x in mini_morph if len(x) > 0]
                    #mini_morph = [x.replace('ː', ':') for x in mini_morph] # wait on this till morphAlign works
                    #o = o.replace('ː', ':')
                    orig2.append(o); morph2.append(mini_morph)
        orig, morph = orig2, morph2'''
        #orig = [x for x in orig if any(y.isalpha() for y in x)]
        #orig = [x for x in orig if not re.search('\w+_\d+', x)] # meant to get rid of 'BI_572', which has no corresponding morphs
        #morph = [[x for x in m if any(y.isalpha() for y in x)] for m in morph]
        #morph = [x for x in morph if len(x) > 0]
    if lang in ['Movima']:
        #orig = [x.translate(trans_dict) for x in orig]; morph = [[x.translate(trans_dict) for x in y] for y in morph]
        '''orig2, morph2 = [], []
        for o, m in zip(orig, morph):
            if not (all(x in punct+'1234567890' for x in o)):
                if not len(o) == 0:
                    orig2.append(o); morph2.append(m)
        orig, morph = orig2, morph2
        for i, m in enumerate(morph):
            morph[i] = ['-' if (x is np.nan or x == '') else x for x in m]'''
    return orig, morph
def drop_nonwords(orig, morph):
    """ INPUT: lists of original words and morph-sequences
        OUTPUT: lists of same, minus any in which the word has no alpha chars"""
    new_orig, new_morph = [], []
    for o, m in zip(orig, morph):
        if any(y.isalpha() for y in o):
            new_orig.append(o); new_morph.append(m)
    return new_orig, new_morph

    #### MAIN FUNCTION ####
def wordAlign(I,spk,g2p_dict):
    # PREPARE OBJECTS
    # language info
    global weird
    weird = I.d['weird']
    d_styp = I.d['d_typ'][spk]
    lang,w_tier,mb_tier = I.d['lang'],d_styp['wd'],d_styp['mb']
    m_tier,ms_tier = d_styp['maus'],d_styp['ms']
        # g2p dict' pre-read (due to speaker loop)
    trans_dict = I.d['tr_dict']
    
    if not w_tier or not m_tier or not ms_tier:
        return None
    if weird:
        fillWDict(I)
    # get tier content
    orig, morph = zip(*get_word_morph_pairs(I.d['eaf'],w_tier,mb_tier))
    maus_orth = [seg.content for seg in m_tier if seg.content != '']
    maus_samp = [seg.content for seg in ms_tier if seg.content != '']
    # CLEAN TIERS
    # language-specific cleaning of original word tier (and morpheme tier?)
    # very important here to exclude word-segments that have no morphological counterparts!
    orig, morph = lang_specific_cleaning(orig, morph, lang, trans_dict)
    orig, morph = drop_nonwords(orig, morph)
    # clean and convert original words/morphemes to SAMPA
    clean_tup = clean_orig_sampa((orig, maus_orth, maus_samp, morph),
                                 lang, g2p_dict,trans_dict,
                                 clean_maus=False, morphs=True)
    cl_orig, cl_orig_samp, cl_maus, cl_maus_samp, cl_morph, cl_morph_samp, orig_morph = clean_tup
    if any(x == '-' for x in cl_orig_samp):
        print('--- wordAlign: "-" present in cl_orig_samp 3')
    # PREPARE INPUT
    # create lists of tuples of unique SAMPA phones, for pairwise word-aligning in nw_matrix()
    if not I.d['weird']:
        cl_maus_samp = [strip_maus(x, g2p_dict, lang) for x in cl_maus_samp]
    orig_tup_set = list(set(tuple(x.split(' ')) for x in cl_orig_samp))
    maus_tup_set = list(set(tuple(x.split(' ')) for x in cl_maus_samp))
    maus_str_set = [' '.join(x) for x in maus_tup_set]
    orig_str_set = [' '.join(x) for x in orig_tup_set]
    # create distance matrix and dictionary (as input for pairwise2.align, using SAMPA phones)
    zeroes = np.zeros((len(maus_tup_set), len(orig_tup_set)))
    dist_matrix = nw_matrix(zeroes, maus_tup_set, orig_tup_set, kind='list', pol=1)
    subst_d = create_subst_dict(maus_str_set, orig_str_set, dist_matrix)

    #RUN ALIGNMENT
    # get alignment and write to file
    if any(x == '-' for x in cl_orig_samp):
        print('--- wordAlign: "-" present in cl_orig_samp 4')
    if any(x == '-' for x in cl_maus_samp):
        print('--- wordAlign: "-" present in cl_maus_samp')
    alignment = pairwise2.align.globalds(cl_maus_samp, cl_orig_samp, subst_d, -1, -1, gap_char=['€'], one_alignment_only=True)
    # take first alignment and realign with original word/morph lists, accounting for alignment gaps
    al_maus, al_orig = alignment[0][0], alignment[0][1]
    new_orig_wd = realign_with_gaps(al_orig, al_maus, orig)
    new_maus_wd = realign_with_gaps(al_maus, al_orig, maus_orth)
    #new_morph = realign_with_gaps(al_orig, al_maus, cl_morph)
    #new_morph_samp = realign_with_gaps(al_orig, al_maus, cl_morph_samp)
    #new_morph_orig = realign_with_gaps(al_orig, al_maus, orig_morph)
    # trying to align everything to al_maus instead of al_orig, because al_maus has gaps that al_orig doesn't
    '''new_orig_wd = realign_with_gaps(al_maus, orig)
    new_maus_wd = realign_with_gaps(al_maus, maus_orth)
    new_morph = realign_with_gaps(al_maus, cl_morph)
    new_morph_samp = realign_with_gaps(al_maus, cl_morph_samp)
    new_morph_orig = realign_with_gaps(al_maus, orig_morph)'''
    return (al_orig,al_maus,new_orig_wd,new_maus_wd)


    ## Turns tier content into tier pointers
def getPntr(I,spk):
    """Let's have actual pointers/ids.
    Returns lists of same length as those in 'tl_lists'."""

    if not I.d['output'][spk]:
        return
        # Variables
    d_styp = I.d['d_typ'][spk]
    w_tier,m_tier,mb_tier = d_styp['wd'],d_styp['maus'],d_styp['mb']
    tl_lists,lang,trans_dict = I.d['output'][spk],I.d['lang'],I.d['tr_dict']
    add = I.d['add']

        # Moar variables
    l_orig,l_maus,l_oOrth,l_mOrth = tl_lists
    l_pmaus = []; l_porig = []
    w_pos = 0; w_max = len(w_tier); m_pos = 0; m_max = len(m_tier)
        # re-run lang_specific_cleaning() MS
    orig, morph = zip(*get_word_morph_pairs(I.d['eaf'],w_tier,mb_tier))
    w_cont,dummy = lang_specific_cleaning(orig, morph, lang, trans_dict)
    for a in range(len(w_cont)):
        w_cont[a] = w_cont[a].strip()
    ls = len(l_maus)
    
        # Main loop (on 'l_maus/etc')
    for a in range(0,ls):
        morth,oorth = l_mOrth[a].strip(), l_oOrth[a].strip()
            # Get MAUS segment reference
        if not morth == add:
            for b in range(m_pos,m_max):
                mseg = m_tier.elem[b]
                if mseg.content.strip() == morth:
                    l_pmaus.append((mseg.name,b,mseg))
                    m_pos = b+1; break
        else:
            l_pmaus.append(("",-1,None))
            # Get ORIG segment reference
        try:
            if not oorth == add:
                for b in range(w_pos,w_max):
                    wseg = w_tier.elem[b]
                    if ((not oorth and not w_cont[b]) or
                        (w_cont[b] == oorth)):
                        l_porig.append((wseg.name,b,wseg))
                        w_pos = b+1; break
            else:
                l_porig.append(("",-1,None))
        except TypeError as te:
            print('TypeError in getPntr:', a-1, l_mOrth[a-1])
            print(a, oorth)
            raise(te)
    return (l_porig,l_pmaus)

