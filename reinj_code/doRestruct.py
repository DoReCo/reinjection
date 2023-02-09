"""24.01.2022 Script for Transcription restructuration.

Handles:
    - adding missing cor tiers (tx,ref,wd)
    - restructuring (and sorting) the Transcription"""

from corflow import Transcription

    # Generic parenting functions
def time_find_children(c_tier, start, end):
    """ find children between start and end
    Note: not used for 'lexify_morph_list'."""
    segs = []; lc = len(c_tier)
    cid,cind,cseg = c_tier.searchByFloat(c_tier,start)
    if not cseg:
        return segs
    while cseg.end <= end:
        segs.append(cseg); cind += 1
        if cind >= lc:
            break
        cseg = c_tier.elem[cind]
    return segs
    
    # Create wd_tier sub-functions
def get_type(text, seps):
    """ input text and separator; return type of morpheme"""
    if len(text) == 0:
        return 'pre' # if empty string, treat as prefix (and join to preceding prefix, if necessary)
    elif text[0] in seps and not text[-1] in seps:
        return 'suf'
    elif text[0] in seps and text[-1] in seps:
        return 'inf'
    elif not text[0] in seps and text[-1] in seps:
        return 'pre'
    else:
        return 'rt'
def start_new_wd(prev, cur, i, seg):
    """ input previous and current morph types; return False if it's not clear
        whether the previous morph is word-final, or True if it is (or '?')"""
    if cur == 'inf':
        ret = False
    elif cur == 'suf':
        ret = False # allowing suffixes to follow prefixes
    elif cur == 'rt':
        ret = False if prev in ['pre', 'inf'] else True
    elif cur == 'pre':
        ret = False if prev == 'pre' else True
    if ret == '?':
        print("--- PREFIX FOLLOWED BY SUFFIX AT INDEX {}: {}!".format(i, seg.content))
        raise
    return ret
def lexify_morph_list(wd_tier,mb_tier):
    """ input list of Segments and list of separators
        iterate and group them into words; return list of lists of Segments"""
        # Variables
    seps = ["-","=","~"]
    start = -1.; end = -1.; text = ""
    id = "wd"; incr = 0
    prev_type = 'suf'; cur_type = None
        # Iterate over morpheme segments
    for i, seg in enumerate(mb_tier.elem):
        cur_type = get_type(seg.content, seps) # pre, suf, inf, rt
        cur_state = start_new_wd(prev_type, cur_type, i, seg) # True if prev morph was word-final, else True
        #print(i, seg.content, prev_type, cur_type, cur_state) # uncomment for testing validity of word merges
        if cur_state == True: # complete word and start new one
            if start >= 0.:
                nseg = wd_tier.create(-1,id+str(incr),start,end,text); incr += 1
            start = seg.start; end = seg.end
            text = seg.content
        else: # keep appending those morphs
            end = seg.end
            text = text+seg.content
        prev_type = cur_type
        # Last word segment
    if start >= 0:
        wd_tier.create(-1,id+str(incr),start,end,text); incr += 1
def tokenize(wd_tier,tx_tier,syms=" "):
    """Split utterance tier into a word tier"""
    id = "wd"; incr = 0
    for seg in tx_tier:
        l_txt = seg.content.split(syms)
        lt = len(l_txt); dur = (seg.end-seg.start)/lt
        for a,txt in enumerate(l_txt):
            s = seg.start+(a*dur)
            e = seg.start+((a+1)*dur)
            nseg = wd_tier.create(-1,id+str(incr),s,e,txt); incr += 1
            nseg.setParent(seg)
    # Create tiers (utterance/words)
def insert_wd_tier(trans,mb_tier,tx_tier,spk,d_core):
    """ input Transcription object, mb tier, tx tier, separator list, speaker
        - create list of word-level lists of mb segs with lexify_morph_list()
        - merge these into new wd Segments
        - create wd Tier and populate with wd Segments
        - replace tx tier's mb child metadata with wd metadata
        - update mb tier's tx parent metadata with wd metadata"""

        # Create new word tier
    wd_tier = trans.create(-1,'wd@'+spk,trans.start,trans.end)
    wd_tier.setMeta('speaker',spk); wd_tier.setMeta('type','wd')
    if mb_tier:
            # Add segments
        lexify_morph_list(wd_tier,mb_tier)
    elif tx_tier:
        tokenize(wd_tier,tx_tier)
        # Add new tier to d_core
    d_core[spk][wd_tier] = {'level':'word',
                       'type':'wd',
                       'speaker':spk,'old_spk':"",
                       'parent':None}
def insert_tx_tier(trans,ctier,ptier,spk,d_core):
    """Create a new utterance tier (based on 'ft' (ptier) and 'wd' (ctier))."""
    
        # Variables
    id = "tx"; incr = 0; pos = 0; lc = len(ctier)
        # Create new utterance tier
    tier = trans.create(-1,'tx@'+spk,ptier.start,ptier.end)
    tier.setMeta('speaker',spk); tier.setMeta('type','tx')
        # Add segments
    pos = 0; lc = len(ctier)
    for pseg in ptier:
        text = ""; cseg = None
        for a in range(pos,lc):                         # Get segments
            cseg = ctier.elem[a]
            if cseg.start >= pseg.end:
                pos = a; break
            elif cseg.start >= pseg.start:
                text = text + " " + cseg.content
        text = text.strip()
        nseg = tier.create(-1,id+str(incr),pseg.start,pseg.end,text); incr += 1
        # Add new tier to d_core
    d_core[spk][tier] = {'level':'utterance',
                    'type':'tx',
                    'speaker':spk,'old_spk':"",
                    'parent':None}
def insert_rf_tier(trans,ptier,spk,d_core):
    """Creates a new reference tier (and adds it to d_core)."""
    
        # Variables
    id = "ref_"; incr = 0; pos = 0; lc = len(ptier)
        # Create new reference tier
    tier = trans.create(-1,'ref@'+spk,trans.start,trans.end)
    tier.setMeta('speaker',spk); tier.setMeta('type','ref')
        # Add segments
    for pseg in ptier.elem:
        nseg = tier.create(-1,id+str(incr),pseg.start,pseg.end,
                           id+"{:04d}".format(incr)); incr += 1
        # Add new tier to d_core
    d_core[spk][tier] = {'level':'utterance',
                    'type':'ref',
                    'speaker':spk,'old_spk':"",
                    'parent':None}
    return tier
def addTiers(trans,d_core,d_else,d_typ,weird=False):
    """Checks if tiers must be created (utterance/words)."""
    
        # For each speaker
        # We create new core tiers and add them to 'd_core'
    for spk, d_tiers in d_typ.items():
        if not weird and not d_tiers['ft']:         # no 'ft', invalid spk
            if not spk in d_else:
                d_else[spk] = {}
            for typ,tier in d_tiers.items():        # move tiers to 'd_else'
                if not tier:
                    continue
                d_else[spk][tier] = d_core[spk][tier].copy()
                d_core[spk].pop(tier)
                d_tiers[typ] = None
            continue
            # Check utterance case
        if (not d_tiers['tx']) and d_tiers['wd']:
            tier = insert_tx_tier(trans,d_tiers['wd'],d_tiers['ft'],spk,d_core)
            d_tiers['tx'] = tier
            # Check word case
        elif (not d_tiers['wd']) and (d_tiers['mb'] or d_tiers['tx']):
            tier = insert_wd_tier(trans,d_tiers['mb'],d_tiers['tx'],spk,d_core)
            d_tiers['wd'] = tier
            # Check reference case
        if not d_tiers['ref']:
            ptier = d_tiers['tx']
            if not ptier: # Case?
                ptier = d_tiers['ft']
            tier = insert_rf_tier(trans,ptier,spk,d_core)
            d_tiers['ref'] = tier


    # Rebuilds the Transcription's hierarchy and ordering
def restruct(I):
    """Rebuilds the Transcription's hierarchy and ordering"""

    def cleanTi(tier,typ,spk,rename=True):
        tier.setParent(None,False,False)
        tier.clearChildren(False,False)
        tier.setMeta('speaker',spk); tier.setMeta('type',typ)
        if rename:
            tier.name = typ+"@"+spk
    def ch_par(trans,tier,key,d_ch):
        if key in d_ch:
            par = d_ch[key]
            if par:
                tier.setParent(par)
            return True
        else:
            return False
    
        # Variables
    trans,d_core,d_else = I.d['eaf'],I.d['d_core'],I.d['d_else']
    d_styp = I.d['d_typ']; d_par = {}
        # For each core tier
    l_tmp = []; debug = 0
    for spk,d_tiers in d_core.items():
        if not d_tiers:     # Actual case
            continue
            # We need two dicts
        d_par = {'ref':None, # initialize parent relation dict (tier type:tier)
                 'wd':d_styp[spk]['ref'],
                 'mb':d_styp[spk]['wd']}
        d_typ = {'utterance':d_styp[spk]['ref'], # initialize level-tier dict
                 'word':d_styp[spk]['wd'],
                 'morpheme':d_styp[spk]['mb']}
            # We clean the tiers
        for tier,d_vals in d_tiers.items():
            cleanTi(tier,d_vals['type'],spk)
            # Then we restructure
        for tier,d_vals in d_tiers.items():
                # And we give them parents
            if not ch_par(trans,tier,d_vals['type'],d_par):
                ch_par(trans,tier,d_vals['level'],d_typ)
            # Then we reorder
        top_tier = d_styp[spk]['ref']; l_tmp.append(top_tier)
        l_par = top_tier.children(); w_tier = d_styp[spk]['wd']
        if w_tier in l_par:
            l_par.remove(w_tier); l_par.append(w_tier)
        while l_par:
            l_cop = []
            for child in l_par:
                l_tmp.append(child)
                l_cop = l_cop+child.children()
            l_par = l_cop
            # Make sure MAUS tiers appear
        for tier,d_vals in d_tiers.items():
            if d_vals['level'] == "phoneme":
                l_tmp.append(tier)
    
        ## This is where we re-add non-core tiers
        # Get a dict of core tiers by name
    d_names = {}; incr = 0
    for spk,d_tiers in d_core.items():
        d_names[spk] = {}
        for tier in d_tiers:
            d_names[spk][tier.name] = tier
        # Loop through non-core tiers
    for spk,d_tiers in d_else.items():
            # Clean old structure
        for tier,d_vals in d_tiers.items():
            typ = d_vals['type']
            if (not typ or typ == 'nan') and tier.checkMeta('type'):
                typ = tier.meta('type')
            cleanTi(tier,typ,spk,False)
        for tier,d_vals in d_tiers.items():
            # Re-set the structure and add to 'l_tmp'
                # We look for an old speaker
            name = tier.name
            if 'old_spk' in d_vals:
                old_spk = d_vals['old_spk']
                    # We try and "safely" replace the old speaker
                if "@" in name:
                    typ,spk = name.rsplit("@",1)
                    spk = spk.replace(old_spk,spk)
                    name = typ+"@"+spk
                # We prevent duplicate names
            if spk in d_names and name in d_names[spk]:
                if "@" in name:
                    typ,spk = name.rsplit("@",1)
                else:
                    typ,spk = name,"@"
                typ = typ+"_old"+str(incr); incr += 1
                name = typ+"@"+spk
            if not name == tier.name:
                tier.name = name
                # And we give them their old parent
            par = d_vals['parent']
            pn = par.name if par else "None"
            if par:
                tier.setParent(par)
            l_tmp.append(tier)
    trans.elem = l_tmp; del l_tmp   # replace 'elem' completely
        # Last restructuration for the segments
    id = "a"; incr = 0
    for a,tier in enumerate(trans):         # reset indexes
        trans.d_elem[tier][0] = a
        tier.d_elem = {}                    # reset segment structure
        for a,seg in enumerate(tier):       # fill with default (+rename)
            seg.name = id+str(incr); incr += 1
            tier.d_elem[seg] = [a,None]
    l_par = trans.getTop(); l_tmp = []
    while l_par:                            # set parent segments
        for ptier in l_par:
            l_child = ptier.children(); l_tmp = l_tmp+l_child
            for ctier in l_child:
                for a in range(len(ctier)-1,-1,-1):
                    cseg = ctier.elem[a]
                    mid = cseg.start+((cseg.end-cseg.start)/2)
                    pseg = cseg.getTime(mid,ptier)
                    if not pseg:
                        ctier.pop(a)
                    else:
                        cseg.setParent(cseg.getTime(mid,ptier))
        l_par = l_tmp.copy(); l_tmp = []
    for tier in trans:
        tier.sortByTime()
        if tier.name.startswith("wd@"):
            tier.setChildTime(ch=False)