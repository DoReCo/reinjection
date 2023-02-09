# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import os,re,pandas
from string import punctuation
from Bio import pairwise2

punct = ''.join([x for x in punctuation if x != "'" and x != "`" and x!= ':'])
punct = punct + ' '
weird = False; d_weird = {}

    # Saddest function ever
def fillWDict(I):
    """Fills the dictionary 'd_weird' for WEIRD languages."""
    path = os.path.join(I.d['in_dir'],I.d['eaf'].name+".csv")
    with open(path,'r',encoding="utf-8") as f:
        for line in f:
            wd,ws = line.split(";")
            wd = strip_text(wd,punct)
            if not wd in d_weird:
                d_weird[wd] = ws.strip()


    #### ALIGNING ####
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

def strip_text(text, punct):
    """ may need to specify for language (e.g. glottal stop apostrophe in Arapaho)"""
    text = str(text)
    for char in punct:
        text = text.replace(char, '')
    return text.lower()

def parse_orthography(text, dff, ch_punct=False, V=False, sep_char = '-'):
    """ INPUT: text string, g2p df, list of punctuation to strip
        - loop through characters in text string, searching for longest possible
        character substrings first, and identifying it as C or V
        - also create list of phoneme substrings, stored as dash-separated string
        OUTPUT: dash-separated phoneme string and CV string
        
        TODO: deal with orthographic strings that correspond to two MAUS phones"""
    if ch_punct != False:
        text = strip_text(text, punct)
    phonemes, maus_phones, cvs = [], [], []
    text_len = len(text)
    if weird:
        if not "****" in text:
            result = d_weird.get(text.strip().lower()).strip()
        else:
            result = ""
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

def clean_str(str_, punct):
    try:
        new_str = str_.translate(str.maketrans(dict.fromkeys(punct)))
        return new_str
    except AttributeError as e:
        print(str_)
        raise e
    return str(str_).translate(str.maketrans(dict.fromkeys(punct)))

def create_morph_phone_index_df(morph_list, kind='list', sep=' '):
    """ INPUT: list of morphs of a word, as lists of SAMPA phones
        - create df of morphs, phones, indices of morphs (m_ind) and phones
        (p_ind), and indices of phones within morphs (pm_ind)
        OUTPUT: df"""
    morph_list = [phone for phone in [morph.split(' ') if kind == 'string' else morph for morph in morph_list]]
    m, m_ind, p, pm_ind = [], [], [], []
    for i, morph in enumerate(morph_list):
        for j, phone in enumerate(morph):
            m.append(morph)
            m_ind.append(i)
            p.append(phone)
            pm_ind.append(j)
    ret_df = pandas.DataFrame({'morph':m, 'morph_ind':m_ind, 'phone':p,
                               'morph_phone_ind':pm_ind})
    ret_df['phone_ind'] = range(len(ret_df))
    return ret_df

def check_alignment_errors(seq1, seq2, gap_char='€'):
    """ INPUT: aligned sequences from pairwise2 (MAUS first, original second)
        OUTPUT: counts of gaps in seq1 and seq2, and of replacements"""
    gap1, gap2, repl = 0, 0, 0
    for char1, char2 in zip(seq1, seq2):
        if char1 == gap_char and char2 == gap_char: raise # shouldn't be any '-' matches
        elif char1 == gap_char: gap1 += 1
        elif char2 == gap_char: gap2 += 1
        elif char1 != char2: repl += 1
    return gap1, gap2, repl

def join_morph_df_and_maus(mdf, maus_align, orig_align, gap_char='€'):
    """ INPUT: morph_df (from create_morph_phone_index_df()) and the MAUS and
        original alignment sequences from pairwise2
        OUTPUT: morph_df integrated with MAUS phones, with gaps in appropriate
        places"""
    mdf = mdf[:]
    errors = []; uncertains = {k:[] for k in mdf['morph_ind'].unique()}; col_num = len(mdf.columns)
    c = 0 # counter to adjust for multiple row insertions
    for i, (m, o) in enumerate(zip(maus_align, orig_align)):
        if o == gap_char:
            mdf.loc[i - c - (.5/(1+c))] = [np.nan] * col_num # insert new row; weird math to handle adjacent NaNs
            errors.append('OGAP')
            c += 1
        elif m == gap_char:
            errors.append('MGAP')
        elif m != o:
            errors.append('REPL')
        else:
            errors.append('')
    mdf = mdf.sort_index().reset_index(drop=True)
    # check for OGAP NaN rows and fill 'morph_ind' col appropriately
    if not all(x == gap_char for x in maus_align): # skip full NaN words
        first_ind, last_ind = 0, len(mdf) - 1
        for i, row in mdf.iterrows():
            if last_ind == 0: # skip single-phone words
                break
            if pandas.isnull(row['morph_ind']): # check NaN rows
                if i == first_ind:
                    next_ind = mdf.loc[i + 1, 'morph_ind']
                    if pandas.isnull(next_ind): # if following morph_ind is ALSO NaN, assign morph_ind = 0
                        mdf.loc[i, 'morph_ind'] = 0
                    else: # otherwise assign morph_ind = next_ind
                        mdf.loc[i, 'morph_ind'] = next_ind
                elif i == last_ind: # shouldn't be possible for previous morph_ind to be NaN (already been processed)
                    prev_ind = mdf.loc[i - 1, 'morph_ind']
                    mdf.loc[i, 'morph_ind'] = prev_ind
                else: # logically, the remaining words should have at least three phones
                    prev_ind, next_ind = mdf.loc[i - 1, 'morph_ind'], mdf.loc[i + 1, 'morph_ind']
                    mdf.loc[i, 'morph_ind'] = prev_ind
                    # flag adjacent morpheme ids as uncertain
                    if not prev_ind == next_ind:
                        if not pandas.isnull(next_ind):
                            uncertains[prev_ind].append('end')
                            uncertains[next_ind].append('start')
                        #else:
                        #    uncertains[prev_ind].append('end')
    mdf['phone_orig'] = orig_align; mdf['phone_maus'] = maus_align
    mdf['error'] = errors
    # clean up uncertains
    ret_uncertain = []
    for k, v in uncertains.items():
        v = list(set(v))
        if len(v) > 1:
            ret_uncertain.append((k, 'start/end')) # both boundaries uncertain
        elif len(v) == 1:
            ret_uncertain.append((k, v[0])) # just one boundary uncertain
    return mdf, ret_uncertain

def morph_format(mdf, col, inner_sep=' ', outer_sep=' '):
    """ INPUT: morph_df and column to format (phone_orig or phone_maus)
        OUTPUT: string of SAMPA chars with brackets around morphemes"""
    strs = []
    for name, group in mdf.groupby('morph_ind'):
        inner = inner_sep.join(group[col].tolist())
        strs.append(f"[{inner}]")
    return outer_sep.join(strs)

def create_WtP_aligned_df(dff):
    """ INPUT: word-aligned pandas df
        - loop through all words in word-aligned df
        - align phones to morph-level SAMPA
        OUTPUT: same pandas df with new columns for phone-to-morph alignments,
        problems, and uncertainty info"""
    # initialize new columns
    cols = ['maus_align', 'orig_align', 'num_morph', 'num_acc', 'morph_issue',
            'maus_gap', 'orig_gap', 'repl', 'num_align', 'uncertain']
    col_dict = {x:[] for x in cols}
    # loop through all words in df
    for i, row in dff.iterrows():
        # skip rows with missing data
        ## fill skipped rows with NaNs
        if (row['maus_samp'] == '-' or row['orig_samp'] == '-' or row['morph_samp'] == '-'):
            for col in cols:
                col_dict[col].append(np.nan)
            continue
        # create indexed morph-phone df
        lmdf = create_morph_phone_index_df(row['morph_samp'].split(','), kind='string')
        maus_input, orig_input = row['maus_samp'].split(' '), lmdf['phone'].tolist()
        # perform alignment and join aligned MAUS string back into morph-phone df
        al = pairwise2.align.globalms(maus_input, orig_input, 2, -1, -1, -1, gap_char=['€'])
        maus_al, orig_al = al[0][0], al[0][1]
        lmdf2, uncertains = join_morph_df_and_maus(lmdf, maus_al, orig_al)
        # get problems
        maus_gap, orig_gap, repl = check_alignment_errors(maus_al, orig_al)
        problem_morphs = [str(int(i)) for i, group in lmdf2.groupby('morph_ind') if any(group['error'].tolist())]
        # update columns
        col_dict['maus_align'].append(morph_format(lmdf2, 'phone_maus'))
        col_dict['orig_align'].append(morph_format(lmdf2, 'phone_orig'))
        col_dict['num_morph'].append(lmdf2['morph_ind'].max() + 1)
        col_dict['num_acc'].append(col_dict['num_morph'][-1] - len(problem_morphs))
        col_dict['morph_issue'].append(','.join(problem_morphs))
        col_dict['maus_gap'].append(maus_gap); col_dict['orig_gap'].append(orig_gap)
        col_dict['repl'].append(repl)
        col_dict['num_align'].append(len(al))
        col_dict['uncertain'].append(uncertains)
    # make new df and concatenate with old df
    new_col_df = pandas.DataFrame(col_dict)
    new_df = pandas.concat([dff, new_col_df], axis=1)
    return new_df

def create_morph_phone_id_df(dff):
    """ INPUT phone-to-morph aligned df
        - loop through words, morphs, and phones
        - assign content, id for all three, as well as uncertainty strings
        OUTPUT: phone-level df with id/content columns for phones/morphemes
        and uncertainty content"""
    def parse_maus_align(string):
        return [x.split(' ') for x in re.findall("\[(.*?)\]", string)]
    # initialize new columns
    col_dict = {k:[] for k in ['phone', 'phone_id', 'morph',# 'morph_samp',
                               'morph_id', 'word', 'word_id', 'maus_align',
                               'orig_align', 'uncertainty', 'scope']}
    wd_id, mo_id, ph_id = 0, 0, 0
    # loop through all words
    for i, row in dff.iterrows():
        # if no alignment, append NaNs
        '''if pandas.isnull(row['m_format']):
            col_dict['phone'].append(''); col_dict['morph'].append('')
            col_dict['word'].append(row['maus_samp']); col_dict['phone_id'].append(ph_id)
            col_dict['morph_id'].append(mo_id); col_dict['word_id'].append(wd_id)
            col_dict['uncertainty'].append('');
            col_dict['scope'].append(row['scope'])
            wd_id += 1; mo_id += 1; ph_id += 1
        else:'''
            # get lists of aligned morphs in word
        ma_string = row['m_format']
        morph_samp_lists = parse_maus_align(ma_string)
        morph_orth_lists = row['orig_mb']
        # create uncertainty dict (k=index, v='start' or 'end' content)
        if len(row['uncertain']) > 0:
            unc_d = {k:v for k, v in row['uncertain']}
        else:
            unc_d = {}
        # loop through morphs
        for j, (m_samp, m_orth) in enumerate(zip(morph_samp_lists, morph_orth_lists)):
            content = ''
            if j in unc_d.keys(): # add uncertainty string, if any for this morph
                content = unc_d[j]
            # loop through phones and append new row data for each column
            for k, phone in enumerate(m_samp):
                col_dict['phone'].append(phone)
                col_dict['morph'].append(m_orth)
                #col_dict['morph_samp'].append('[' + ' '.join(m_samp) + ']')
                col_dict['word'].append(row['maus_wd'])
                col_dict['phone_id'].append(ph_id)
                col_dict['morph_id'].append(mo_id)
                col_dict['word_id'].append(wd_id)
                col_dict['maus_align'].append(row['maus_align'])
                col_dict['orig_align'].append(row['orig_align'])
                col_dict['uncertainty'].append(content)
                col_dict['scope'].append(row['scope'])
                ph_id += 1
            mo_id += 1
        wd_id += 1
    ret_df = pandas.DataFrame(col_dict)
    return ret_df

def align_maus_to_morph(maus_wd, orig_wd, maus_mb, orig_mb, scope):
    """"""
    lmdf = create_morph_phone_index_df(maus_mb, kind='string')
    maus_input, morph_input = maus_wd.split(' '), lmdf['phone'].tolist()
    # perform alignment and join aligned MAUS string back into morph-phone df
    al = pairwise2.align.globalms(maus_input, morph_input, 2, -1, -1, -1, gap_char=['€'])
    maus_al, orig_al = al[0][0], al[0][1]
    lmdf2, uncertains = join_morph_df_and_maus(lmdf, maus_al, orig_al)
    # get problems
    maus_gap, orig_gap, repl = check_alignment_errors(maus_al, orig_al)
    problem_morphs = [str(int(i)) for i, group in lmdf2.groupby('morph_ind') if any(group['error'].tolist())]
    ret_d = {'m_format':morph_format(lmdf2, 'phone_maus'),
             'o_format':morph_format(lmdf2, 'phone_orig'),
             'maus_align':' '.join(maus_al), 'orig_align':' '.join(orig_al),
             'maus_wd':maus_wd, 'orig_wd':orig_wd, 'orig_mb':orig_mb,
             'scope':scope,
             'num_morph':lmdf2['morph_ind'].max()+1,
             'morph_issue':','.join(problem_morphs),
             'maus_gap':maus_gap, 'orig_gap':orig_gap, 'repl':repl,
             'num_align':len(al), 'uncertain':uncertains}
    return ret_d

def fill_skipped_morphs(maus_wd, orig_wd, maus_mb, orig_mb, scope):
    """ fill morph_d when there is no phone-to-morph alignment to be done"""
    if len(maus_mb) == 0:
        ret_d = {'m_format':maus_wd, 'o_format':maus_wd, # using maus_wd for both - no sampification applied!
                 'maus_align':maus_wd, 'orig_align':maus_wd,
                 'maus_wd':maus_wd, 'orig_wd':orig_wd, 'orig_mb':orig_mb,
                 'scope':scope, 'num_morph':'', 'morph_issue':'', 'maus_gap':'',
                 'orig_gap':'', 'repl':'', 'num_align':'', 'uncertain':''}
    elif len(maus_mb) == 1:
        maus_align = '[' + maus_wd + ']'
        ret_d = {'m_format':maus_align, 'o_format':maus_align, # using maus_wd for both - no sampification applied!
                 'maus_align':maus_wd, 'orig_align':maus_wd,
                 'maus_wd':maus_wd, 'orig_wd':orig_wd, 'orig_mb':orig_mb,
                 'scope':scope, 'num_morph':'', 'morph_issue':'', 'maus_gap':'',
                 'orig_gap':'', 'repl':'', 'num_align':'', 'uncertain':''}
    return ret_d

def update_morph_to_phone_dict(orig_d, new_d):
    """ update old dict of morph-to-phone lists with info from new word"""
    for k, v in new_d.items():
        orig_d[k].append(v)
    return orig_d


    #### MATCHING ####
def fillLists(tl_lists,df,value):
    """Fixes the damn mess (apparently deprecated?)."""
    
        # We fill new lists from old ones
    l_orth,l_nph,l_nmb,l_nid,l_ncert = [],[],[],[],[]
    l_phone, l_phid, l_morph,l_mbid, l_word,l_woid, l_ocert = tl_lists
        # Fill l_orth ('morph_orth')
    for i,row in df.iterrows():
        orth = row['morph_orth']
        if "," in orth:
            orth = orth.split(",")
        else:
            orth = [orth]
        l_orth = l_orth+orth
        # Fill the new lists
    nid = 0; oid = -1; pos = 0; lm = len(l_orth); ocert = ""
    for a in range(len(l_morph)):
        ph,id,mb,cert = l_phone[a],l_mbid[a],l_morph[a],l_ocert[a]
        if not mb or mb == value:
            continue
            # Add missing morphemes (hopefully)
        if id > oid:
            while not l_orth[pos] == mb:
                if not l_orth[pos] == value:
                    l_nph.append(""); l_nmb.append(l_orth[pos])
                    l_nid.append(nid); nid += 1; l_ncert.append(ocert)
                pos += 1
            oid = id; pos += 1
        l_nph.append(ph); l_nmb.append(mb); l_nid.append(nid); nid += 1
        l_ncert.append(cert); ocert = cert
    return (l_nph,l_nmb,l_nid,l_ncert)
def getMorphPntr(I,spk):
    """Turns lists of indexes into lists of references.
    Also generates a list 'l_cert' with tuples (index,content)."""
    
        # Variables
    d_styp = I.d['d_typ'][spk]
    mb_tier,p_tier = d_styp['mb'],d_styp['ph']
        # No tiers? Move on
    if not mb_tier or not p_tier:
        return ([],[],[])
    l_ph,l_mbid,l_ocert = I.d['output'][spk]
    sym,psym = I.d['sym'],I.d['psym']
    
        # Just associate phonemes with morphemes
        #### OKAY LISTEN
        ## 1. You have a 1-to-1 phoneme correspondance, so each row of your
        #  dataframe is a phone Segment.
        #  2. The 'morph_id' column, in turn, turns out to match morpheme
        #  Segment indexes. 
        #  > Long as this stands true, the code below will work.
    l_pmaus = []; l_porig = []; l_cert = []; i_ph = 0; omb = None
    ob_ind = -1

    for a,ph in enumerate(l_ph):
        if a-i_ph >= len(p_tier):    # Out of range / end of phones
            print("getMorphPntr() out of range:",p_tier.name,len(p_tier),
                  a,i_ph,len(l_ph),len(l_porig),len(l_pmaus),len(l_cert))
            break
        phseg = p_tier.elem[a-i_ph]
        mb_ind = int(l_mbid[a])
        mbseg = mb_tier.elem[mb_ind]
        if ph == I.d['add']:
            i_ph += 1; continue
        l_pmaus.append((phseg.name,a-i_ph,phseg))
        l_porig.append((mbseg.name,mb_ind,mbseg))
        if not omb == mbseg:
            omb = mbseg
            l_cert.append((mbseg,l_ocert[a]))
        #print(l_porig[-1][0], l_porig[-1][1], l_porig[-1][2].content,
        #      l_pmaus[-1][0], l_pmaus[-1][1], l_pmaus[-1][2].content)
    return l_porig,l_pmaus,l_cert

def print_gmp_log(l_porig, l_pmaus, l_cert):
    import pandas
    cols = ['l_porig_id', 'l_porig_i', 'l_porig_c', 'l_pmaus_id', 'l_pmaus_i',
            'l_pmaus_c']
    data = [x + y + z for x, y, z in zip(*(l_porig, l_pmaus))]
    df = pandas.DataFrame.from_records(data, columns=cols)
    for c in df:
        if c in ['l_porig_c', 'l_pmaus_c']:
            df[c] = df[c].apply(lambda x: x.content)
    df.to_excel('getMorphPntr_log.xlsx', index=False)


    #### Main function ####
def getPhones(wseg,a,ph_tier):
    """Returns a string of phones content (space separated)."""
    l_phsegs = []; maus_wd = ""
    for b in range(a,len(ph_tier)): # Get all phone segments
        l_phsegs.append(ph_tier.elem[b].content)
        if ph_tier.elem[b].end >= wseg.end:
            break
    for cont in l_phsegs:           # Return space-separated content
        maus_wd = maus_wd+" "+cont
    return maus_wd.strip()
def morphAlign(I,spk,g2p_dict):
    """Aligns morphemes and phonemes."""

    def iterTier(otier):
        for a in range(len(otier)):
            oseg = otier.elem[a]
            yield (oseg.name,a,oseg)
    
        # Variables
    trans,lang = I.d['eaf'],I.d['lang']; d_styp = I.d['d_typ'][spk]
    w_tier,ms_tier,mb_tier = d_styp['wd'],d_styp['ms'],d_styp['mb']
    ph_tier = d_styp['ph']
    
        # Lack tiers? Move on
    if not w_tier or not ms_tier or not mb_tier:
        return ([],[],[])
    global weird
    weird = I.d['weird']
    if weird:
        fillWDict(I)
    morph_d = {k:[] for k in ['m_format', 'o_format','maus_align',
                              'orig_align','maus_wd', 'orig_wd',
                              'orig_mb', 'scope', 'num_morph', 'morph_issue',
                              'maus_gap', 'orig_gap', 'repl', 'num_align',
                              'uncertain']}
    value = '-';
    punct = ''.join([x for x in punctuation if not x in ["'","`",':']]) # should be lg-specific dict
    punct = punct + ' '
    # Iterate over words
    for wid,wind,wseg in iterTier(w_tier):
        # Assume always at least one child
        l_mbchild = wseg.childDict()[mb_tier]
        scope = (wind,wind+1) # (MAUS word ind, MAUS word ind + 1)
        l_mb = []; l_mbsamp = [] # list of mb content strings and sampified strings
        for mb_seg in l_mbchild:
            l_mb.append(mb_seg.content)
            # don't try to parse MAUS tags
            if len(mb_seg.content) > 0 and not mb_seg.content[0] == '<':
                cl_morph = clean_str(mb_seg.content.lower(), punct)
                cl_morph = parse_orthography(cl_morph, g2p_dict.get(lang), 
                                                 sep_char=' ')[1]
            else:
                cl_morph = ''
            l_mbsamp.append(cl_morph)
        # MAUS sampa (words)
        ph_n,ph_ind,ph_seg = wseg.getTime(wseg.start,ph_tier,det=True)
        maus_wd = ""
        if ph_seg:
            maus_wd = getPhones(wseg,ph_ind,ph_tier)
        if len(l_mbchild) < 2:
                # Deal with other cases
            #### POST_PROCESS CASES (<<wip>>, etc.) ####
            ############################################
                # Else ignore and continue
            # fill and update morph_d (no alignment needed)
            new_d = fill_skipped_morphs(maus_wd, wseg.content, l_mbsamp, l_mb, scope)
            morph_d = update_morph_to_phone_dict(morph_d, new_d)
            scope = False
        if not scope:
            continue
            # Obtain lists of strings (for morphemes 'l_mb' 
            #                          and MAUS sampa 'l_ms')

            ## Variables at disposal at that point
            # lang
            # Inj.trans, w_tier (MAUS words),ms_tier (MAUS SAMPA), 
            #            mb_tier (morphemes)      
            #           (the Transcription and 3 tier objects)
            # 'l_mbchild','l_mschild'                           (list of Segments (morphemes for 'mbchild' and MAUS SAMPA words for 'mschild')
            # wseg,scope                                        (technicalities)
            # l_mb,l_ms                                         (lists of strings, all morpheme contents and SAMPA content for that sequence)

        # in the output, I would need 'scope','l_mb' and 4 lists 
        # (align morph/phones,morph_indexes and uncertainty)
        ### INSERT FUNCTION HERE ####
        #############################
        # phone-to-morph alignment
        # FIX THIS INPUT (like monomorphs)
        new_d = align_maus_to_morph(maus_wd, wseg.content, l_mbsamp, l_mb, scope) # are these Segments? might need to extract content
        # update morph_d
        morph_d = update_morph_to_phone_dict(morph_d, new_d)
    morph_align_df = pandas.DataFrame(morph_d)
    phone_df = create_morph_phone_id_df(morph_align_df)
    # convert phone-level df to list of tuples
    tl_lists = (phone_df['phone'].tolist(),phone_df['morph_id'].tolist(),
                phone_df['uncertainty'].tolist())
    return tl_lists
