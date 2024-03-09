"""24.01.2022 - Master script for reinjection.

Handles:
    - merging MAUS tiers.
    - adding missing cor tiers (tx,ref,wd)
    - restructuring (and sorting) the Transcription
    - word-alignement
    - morph-alignment
Relies on:
    - doRestruct
    - wordAlign,restructWord
    - morphAlign,restructMorph
Uses an 'Inj' class for all operations."""

    # More modules are imported 'locally', on demand
from corflow import fromElan,toElan,fromPraat,toPraat
import sys,os,re,time,pandas,math
pandas.options.mode.chained_assignment = None

    # Variable container
class Inj:
    """A container class for all operations."""
    
    def __init__(self,d_vals={},set=True):
        self.home = os.path.abspath(os.path.dirname(__file__))
        self.d = \
        {   # paths
         'in_dir':"",           # (path) input folder
         'out_dir':"",          # (path) output folder
         'g2p':"",              # (path grapheme-to-phoneme CSV
         'tiers_csv':"",        # (path) tier_metadata CSV
         'log_path':"",
         'spk_csv':"",          # (path) file output for doRestruct
         'wd_csv':"",           # (path) file output for wordAlign
         'mb_csv':"",           # (path) file output for morphAlign
            # transcription
         'eaf':None,            # (pntr) ELAN Transcription object
         'tgd':None,            # (pntr) TextGrid Transcription object
         'weird':False,         # (bool) if WEIRD language (no 'ft')
         'd_tiers':{},          # (dict) all core tiers by speaker by file
         'd_core':{},           # (dict) core tiers by speaker
         'd_else':{},           # (dict) non-core tiers by speaker
         'd_typ':{},            # (dict) core tiers by speaker by type
            # info
         'l_ch':[],             # (list) list of all labels (getTYP)
         'd_tgd':{},            # (dict) labels for the MAUS tiers
         'lang':"",             # (str) the language name
         'd_lang':{},           # (dict) dictionary for languages
         'tr_dict':{},          # (dict) a dict' for wordAlign
         'output':None,         # (tpl<lst<...>>) output from word/morphAlign
            # symbols
         'add':"",              # (str) symbol for filling lists
         'sym':"",              # (str) symbol for non-pause filler segments
         'psym':""}             # (str) symbol for pause filler segments
            # We set default values
        if set:
            setInj(self)
            # We set custom values
        if d_vals:
            setInj(self,d_vals)

def log(log_path="",text=""):                       # logging
    """Log function. This one requires a 'path' argument."""
    if not log_path:
        if text.endswith("\n"):
            text = text[:-1]
        print(text)
        return
    elif not os.path.isfile(log_path):
        f = open(log_path,'w',encoding="utf-8"); f.close()
    with open(log_path,'a',encoding="utf-8") as f:
        f.write(text)
    # Get data
def setInj(I=None,d_vals={}):
    """Filling whatever variables we want in 'Inj()'."""
    if not I:       # We need an object to fill
        I = Inj()
    if not d_vals:  # We can set default values here
        d_vals = {'in_dir':"input",
                  'out_dir':"output",
                  'g2p':"g2p_sample_CV.xlsx",
                  'tiers_csv':"tier_metadata.xlsx",
                  'log_path':'',
                  'weird':False,
                  'wd_csv':"",
                  'mb_csv':"",
                  'l_ch':['ref','tx','ft','wd','mb','gl','ps',
                          'maus','ms','ph'],
                  'd_tgd':{"W_O_":'maus',
                           "W_S_":'ms',
                           "P_S_":'ph'},
                  'add':"€",
                  'sym':"****",
                  'psym':"<p:>"}
    for key,val in d_vals.items():
        if not key in I.d:
            continue
        I.d[key] = val
        # Sets paths (dealing with relative paths)
    l_paths = ['g2p','tiers_csv']
        # Folders (relying on 'Inj' script location)
    if not os.path.isdir(I.d['in_dir']):
        I.d['in_dir'] = os.path.join(I.home,I.d['in_dir'])
    if not os.path.isdir(I.d['out_dir']):
        I.d['out_dir'] = os.path.join(I.home,I.d['out_dir'])
        # Files (according to 'l_paths')
    for path in l_paths:
        nval = val = I.d[path]
        if not os.path.isfile(nval):            # try from script location
            nval = os.path.join(I.home,val)
            if not os.path.isfile(nval):        # try from input folder
                nval = os.path.join(I.d['in_dir'],val)
        I.d[path] = nval
    return I
def resetInj(I):
    """In case of a loop, we want to reset some transcription info."""
    I.d['eaf'] = None; I.d['tgd'] = None
    I.d['d_core'] = {}; I.d['d_else'] = {}
    I.d['output'] = None
def readCSV(path):
    """Creates a pandas DataFrame."""
    if not os.path.isfile(path):
        return None
    _,ext = os.path.splitext(path); df = None
    if ext.lower() == ".csv":
        df = pandas.read_csv(path,keep_default_na=False)
    else:
        df = pandas.read_excel(path,keep_default_na=False)
    return df
def writeCSV(data,path,sheet_name=""):
    """Writes a(n) CSV/Excel file from a dict."""
    
    _,ext = os.path.splitext(path)
    out = pandas.DataFrame(data=data)
    if ext.lower() == ".csv":
        out.to_csv(path)
    elif sheet_name:
        out.to_excel(path,sheet_name=sheet_name)
    else:
        out.to_excel(path)
def getTGD(dir,name,compl="_Algn",l_path=""):
    """Gets a TextGrid file.
    'name' is expected to have no extension."""
    
    if name.lower().endswith(".textgrid"):
        name,ext = os.path.splitext(name)
    else:
        ext = ".TextGrid"
    tgd = os.path.join(dir,name+compl+ext)
    if not os.path.isfile(tgd):
        if l_path:
            log(l_path,"TGD path invalid: {}".format(tgd))
        return None
    elif l_path:
        log(l_path,"TGD found.")
    return fromPraat.fromPraat(tgd)
def getEAF(dir,name):
    """Gets your ELAN file."""
    
    if not name.lower().endswith(".eaf"):
        name = name+".eaf"
    eaf = os.path.join(dir,name)
    if not os.path.isfile:
        return None
    trans = fromElan.fromElan(eaf)
    for tier in trans:
        tier.sortByTime()
    return trans
def exportEAF(trans,dir,name):
    """Writes an ELAN file."""
    
        # We get rid of empty tiers
    a = len(trans)-1
    while True:
        if a < 0:
            break
        elif a > len(trans)-1:
            a = a-1; continue
        tier = trans.elem[a]
        if "doreco-mb-algn@" in tier.name:
            a = a-1; continue
        elif "wd_old@" in tier.name:
            trans.pop(a)
        elif not tier.elem:
            trans.pop(a)
        a = a-1
    """for a in range(len(trans)-1,-1,-1):
        tier = trans.elem[a]
        if "doreco-mb-algn@" in tier.name:
            continue
        elif not tier.elem:
            trans.pop(a)"""
    
        # We also hate the HEADER
    for key,l_val in trans.iterMeta('elan',ch_list=True):
        if key == "LINKED_FILE_DESCRIPTOR":
            trans.metadata['elan'].remove("LINKED_FILE_DESCRIPTOR")
        elif key == "MEDIA_DESCRIPTOR":
            ch_header = True
            for a in range(len(l_val)-1,-1,-1):
                val = l_val[a]
                if  ch_header and trans.name+".wav" in val:
                    ch_header = False; continue
                l_val.pop(a)

        # We check the output
    if not name.lower().endswith(".eaf"):
        name = name+".eaf"
    path = os.path.join(dir,name)
    trans.renameSegs()
    toElan.toElan(path,trans)
    name = name.split(".eaf",1)[0]+".TextGrid"
    path = os.path.join(dir,name)
    toPraat.toPraat(path,trans)
def fillSPK(file,fn="file",ft="name",fs="speaker",fl="tier_label",
                 fc="tier_core",fv="tier_level",fy="tier_type",
                 fa="language",l_path=""):
    """Gets all rows from the 'tier_metadata' excel file."""
    
    start = time.time()
    d_tiers = {}
    d_file = pandas.read_excel(file,dtype=str)
    for index, row in d_file.iterrows():
        tr = str(row[fn]); ti = str(row[ft])
        if tr.lower().endswith(".eaf"):
            tr,_ = os.path.splitext(tr)
        value = {'spk':str(row[fs]),'label':str(row[fl]),'core':str(row[fc]),
                 'level':str(row[fv]),'type':str(row[fy]),'lang':str(row[fa])}
        if tr not in d_tiers:
            d_tiers[tr] = {}
        d_tiers[tr][ti] = value.copy()
    end = time.time()-start
    if l_path:
        log(l_path,"Tier metadata loaded: {:04f}s\n".format(end))
    return d_tiers
def _getOldSpk(tier):
    """Retrieve speaker from tier name."""
    
    old_spk = "#"
        # @ case
    if "@" in tier.name:
        old_spk = tier.name.rsplit("@",1)[1]
        # _ case
    elif "_" in tier.name:
        old_spk = tier.name.split("_",1)[0]
    return old_spk
def getSPK(d_tiers,trans,l_ch):
    """We want the actual tiers in an actual Transcription.
    Note: we want to distinguish between core/non-core tiers.
    Note: we still take non-cored 'ref' tiers into account."""
    
    def addSPK(d_spk,tier,spk,d_vals):
        if not spk in d_spk:
            d_spk[spk] = {tier:d_vals.copy()}
        else:
            d_spk[spk][tier] = d_vals.copy()
    d_nan = {'tx':'utterance','ft':'utterance','wd':'word',
             'mb':'morpheme'}
    
    d_core = {}; d_else = {}; lang = ""; max = len(trans)
        # For each tier
    for a,tier in enumerate(trans):
        old_spk = tier.meta('speaker')
        if not old_spk:
            old_spk = _getOldSpk(tier)
            # Not in 'tier_metadata'
        if ((not d_tiers) or
            (not trans.name in d_tiers) or
            (not tier.name in d_tiers[trans.name])):
                # MAUS case
            """if ("@" in tier.name and tier.name.split("@",1)[0] in l_ch
                and a >= max-3):
                typ,spk = tier.name.split("@",1); core = "yes"
                level = "unknown"
                if typ and not typ == "nan":
                    tier.setMeta('type',typ,i=0)
                    level = d_nan.get(typ,"unknown")
                    if level == "unknown":
                        core = "no"
                d_vals = {'level':level,'type':typ,'old_spk':spk,
                          'parent':tier.parent()}"""
                # Others
            spk = old_spk; core = "no"
            tier.setMeta('type',"",i=0)
            d_vals = {'level':"",'type':"",'old_spk':old_spk,
                      'parent':tier.parent()}
            # In 'tier_metadata'
        else:
            d_vals = d_tiers[trans.name][tier.name]
            lang = d_vals['lang']
            spk = d_vals['spk']; tier.setMeta('speaker',spk,i=0)
            core = d_vals['core']; typ = d_vals['label']
            if typ and not typ == "nan":    # Pandas 'numpy.nan'
                if (not core == 'yes') and typ in l_ch:
                    typ = ""
                tier.setMeta('type',typ,i=0)
            old_typ = d_vals['type']
            if (not core == 'yes') and old_typ == 'ref':
                core = 'yes'; typ = old_typ
                if typ and not typ == "nan":
                    tier.setMeta('type',typ,i=0)
            d_vals = {'level':d_vals['level'],'type':typ,'old_spk':old_spk,
                      'parent':tier.parent()}
        if core == 'yes':
            addSPK(d_core,tier,spk,d_vals)
        else:
            typ = tier.meta('type')
            if typ and typ in l_ch:
                tier.setMeta('type',"o"+typ)
            addSPK(d_else,tier,spk,d_vals)
    return d_core,d_else,lang
def getTYP(d_ti,l_ch=['tx','ft','ref','wd','mb','maus','mb','ph']):
    """Rework 'd_spk' by speaker-then-type."""
    
    d_spk = {}
    for spk,d_tiers in d_ti.items():
        d_spk[spk] = {}
        for tier,d_vals in d_tiers.items():
            typ = tier.meta('type',empty=d_vals['type'])
            d_spk[spk][typ] = tier
            # We need to ensure all labels are in there
    for spk, d_tiers in d_spk.items():
        for typ in l_ch:
            if typ not in d_tiers:
                d_spk[spk][typ] = None
    return d_spk
def getNewTYP(I):
    """Same as 'getTYP()' but from 'trans'."""
    d_spk = {}
    for tier in I.d['eaf']:
        spk = tier.meta('speaker',empty="#")
        if not spk in d_spk:
            d_spk[spk] = {}
        typ = tier.meta('type')
        if typ and typ in I.d['l_ch']:
            d_spk[spk][typ] = tier
    for spk,d_tiers in d_spk.items():
        for typ in I.d['l_ch']:
            if typ not in d_tiers:
                d_spk[spk][typ] = None
    return d_spk
    # MAUS tiers merging
def corrInput(I):
    """Language-specific corrections."""

    def tierYongningNa(O):  
        """Restructure the YongningNa file."""
        d_typ = getTYP(O.d['d_core'],O.d['l_ch'])
        d_lvl = {'morpheme':'mb','word':'wd','utterance':'tx'}
        for spk,d_tiers in O.d['d_else'].items():
            d_styp = d_typ[spk]
            for tier,d_vals in d_tiers.items():
                typ = d_lvl.get(d_vals['level'])
                if typ:
                    O.d['d_else'][spk][tier]['parent'] = d_styp[typ]
    def spkGoemai(O):
        if not O.d['tgd']:
            return
        for tier in O.d['tgd']:
            one,two,name = tier.name.split("_",2)
            name,spk = name.split("@",1)
            tier.name = one+"_"+two+"_ft@"+spk
    def spkNuu(O):
        if not O.d['tgd']:
            return
        for tier in O.d['tgd']:
            if tier.name.startswith("MAUS_"):
                tier.name = tier.name.split("MAUS_",1)[1]
                l_s = tier.name.split("_")
                if len(l_s) == 3:           # 'maus'
                    tier.name = "W_O_"+l_s[1]
                elif l_s[0] == 'W':         # 'ms'
                    tier.name = "W_S_"+l_s[1]
                elif l_s[0] == 'S':         # 'ph'
                    tier.name = "P_S_"+l_s[1]
    def spkSvan(O):
        """Remove '_ALGN' tier suffixes from Svan TextGrid."""
        if not O.d['tgd']:
            return
        for tier in O.d['tgd']:
            if tier.name.endswith("_ALGN"):
                tier.name = tier.name.rsplit("_ALGN",1)[0]+"-"
    def corr_ref(O,sym,nsym):
        """Replaces a tier 'ref' type to something else."""
        d_tmp = {}
        for spk,d_tiers in O.d['d_core'].items():
            d_tmp[spk] = {}
            for tier,d_vals in d_tiers.items():
                d_tmp[spk][tier] = d_vals.copy()
        for spk,d_tiers in d_tmp.items():
            for tier,d_vals in d_tiers.items():
                if not tier.name.startswith(sym):
                    continue
                O.d['d_core'][spk][tier]['type'] = nsym
                if not spk in O.d['d_else']:
                    O.d['d_else'][spk] = {}
                O.d['d_else'][spk][tier] = O.d['d_core'][spk][tier].copy()
                O.d['d_core'][spk].pop(tier)
    
        # Universal boundaries alignment
    at = -1.
    for tier in I.d['eaf']:
        for seg in tier:
            if at < 0. or seg.end-seg.start < at:
                at = seg.end-seg.start
    I.d['eaf'].fixBounds(abs_tol=(at-0.001))
    
        # A dictionary per correction, with a key per language...
    lang = I.d['lang'].lower()
    d_spk = {'goemai':spkGoemai,
             'nuu':spkNuu,
             'svan':spkSvan,
             'yongning na':tierYongningNa}
    d_ref = {'jejuan':('interlinear','oref'),
             'tabaq':('sound@','oref'),
             'svan':('c@','oref')}
        # Check old speakers
    f = d_spk.get(lang)
    if f:
        f(I)
        # Check ref tiers
    tpl = d_ref.get(lang)
    if tpl:
        corr_ref(I,tpl[0],tpl[1])
    return I    
def mergeTGD(I,edit=True):
    """Merging two files together.
    Also getting rid of abusive pauses."""
    
    def findRef(O,tier,d_ref,rn):
        ref = None
        for otier in O.d['eaf']:
            if otier.name == rn:
                ref = otier; break
        if not ref in d_ref:
            d_ref[ref] = [tier]
        else:
            d_ref[ref].append(tier)
    def fuse(seg1,seg2,mend):
        """Actual fusion."""
        seg1.content = seg1.content+seg2.content
        seg1.end = mend
        seg2.struct.remove(seg2)
    
    if not I.d['tgd']:
        log(I.d['log_path'],"\t\tNo TextGrid found.\n")
        return I
        # Find the original tier
    d_ref = {}; l_tiers = []
    for tier in I.d['tgd']:
        if not "_" in tier.name:
            continue
        findRef(I,tier,d_ref,tier.name.split("_",2)[2])
        # Attribute MAUS tiers to speaker in 'd_core'
    le = len(I.d['eaf'])
    for spk,d_tiers in I.d['d_core'].items():   # for all core tiers
        for rtier in d_ref:
            if not rtier in d_tiers:            # check if tier is there
                continue
            for mtier in d_ref[rtier]:          # add MAUS tiers to speaker
                typ = I.d['d_tgd'].get(mtier.name[:4])
                if not typ:
                    continue
                ntier = I.d['eaf'].add(-1,mtier)
                ntier.setMeta("old_tier",mtier,"tech")
                ntier.name = typ+"@"+spk; l_tiers.append(ntier)
                ntier.setMeta('speaker',spk); ntier.setMeta('type',typ)
                for a in range(len(ntier)-1,-1,-1):
                    seg = ntier.elem[a]
                    if seg.start == seg.end:
                        ntier.elem.pop(a)
                I.d['d_core'][spk][ntier] = \
                   {'type':typ,'parent':None,'level':"phoneme",
                    'old_spk':I.d['d_core'][spk][rtier]['old_spk']}
            break
        # Fix boundaries
    I.d['eaf'].setBounds()
    for tier in l_tiers:
        if not tier.elem:
            continue
        seg = tier.elem[-1]
        if seg.end < tier.end:
            if (not seg.content) or ("<p:>" in seg.content):
                seg.end = tier.end
            else:
                tier.create(-1,"add"+str(len(tier)),seg.end,tier.end,
                            content=I.d['psym'])
        seg = tier.elem[0]
        if seg.start > tier.start:
            if (not seg.content) or ("<p:>" in seg.content):
                seg.start = tier.start
            else:
                tier.create(0,"add0",tier.start,seg.start,content=I.d['psym'])
    if len(I.d['eaf']) == le:
        log(I.d['log_path'],"\t\tNo MAUS tier added.\n")
        # Fix multiple pauses
        # Also merge <<wip>> plz
    if not edit:
        return I
    for tier in l_tiers:
        if len(tier) < 2:
            continue
        for a in range(len(tier)-2,-1,-1):  # multiple pauses
            seg1,seg2 = tier.elem[a],tier.elem[a+1]
            if ((not seg1.content or seg1.content == "<<p:>>") and
                (seg1.content == seg2.content)):
                fuse(seg1,seg2,seg2.end)
        #if len(tier) < 3:
        #    continue
        #for a in range(len(tier)-3,-1,-1):  # WIP
        #    seg1,seg2,seg3 = tier.elem[a],tier.elem[a+1],tier.elem[a+2]
        #    if seg2.content == "<<wip>>" and not tier.name.startswith("ph@"):
        #        fuse(seg2,seg3,seg3.end)
        #        fuse(seg1,seg2,seg2.end)
    
    return I
def addTiers(I):
    """Adds missing core tiers ('doRestruct')."""
    import reinj_code.doRestruct as dR
    
    I.d['d_typ'] = getTYP(I.d['d_core'],I.d['l_ch'])
    dR.addTiers(I.d['eaf'],I.d['d_core'],I.d['d_else'],I.d['d_typ'],
                I.d['weird'])
    return I
def _writeRestruct(I):
    """Sub-function to write CSV/Excel output from 'restruct()'."""
    
        # Get path
    path = os.path.join(I.d['out_dir'],I.d['eaf'].name+"_"+I.d['spk_csv'])
        # We need a dict' for pandas' DataFrame
    d_out = {'name':[],'speaker':[],'tier_label':[],'parent':[]}
        # We loop over 'd_core' to fill 'd_out'
    for spk,d_tiers in I.d['d_core'].items():
        for tier,d_vals in d_tiers.items():
            d_out['name'].append(tier.name)
            d_out['speaker'].append(tier.meta('speaker'))
            d_out['tier_label'].append(tier.meta('type'))
            if tier.parent():
                d_out['parent'].append(tier.parent().name)
            else:
                d_out['parent'].append("")
        # Pandas export according to 'spk_csv' extension
    writeCSV(d_out,path,"d_core")
def restruct(I):
    """Rebuilds the tier hierarchy and ordering ('doRestruct')."""
    import reinj_code.doRestruct as dR
    
    I.d['d_typ'] = getTYP(I.d['d_core'],I.d['l_ch'])
    dR.restruct(I)

        # Write 'd_core' to 'spk_csv'
    if I.d['spk_csv']:
        _writeRestruct(I)
    return I
def _readWord(I):
    """Sub-function to read CSV/Excel input for 'aWord()'."""
    
        # Get path
    path = os.path.join(I.d['in_dir'],I.d['eaf'].name+"_"+I.d['wd_csv'])
    if not os.path.isfile(path):
        I.d['output'] = None
        return
    df = readCSV(path)     # pre-aligned words
        # We need a dict' for pandas' DataFrame
    d_spk = {}; l_rows = ['orig_samp','maus_samp','orig_ortho','maus_ortho']
    for a,row in df.iterrows():
        spk = str(row['speaker'])
        if not spk in d_spk:
            d_spk[spk] = ([],[],[],[])
        for b,k in enumerate(l_rows):
            val = row[k] if not pandas.isnull(row[k]) else ""
            d_spk[spk][b].append(str(val))
    I.d['output'] = d_spk
def _writeWord(I):
    """Sub-function to write CSV/Excel output from 'aWord()'."""
    
        # Get path
    path = os.path.join(I.d['out_dir'],I.d['eaf'].name+"_"+I.d['wd_csv'])
        # We need a dict' for pandas' DataFrame
    d_out = {'speaker':[],'orig_samp':[],'maus_samp':[],
             'orig_ortho':[],'maus_ortho':[]}
        # We loop over 'output' to fill 'd_out'
    l_out = ['orig_samp','maus_samp','orig_ortho','maus_ortho']
    lo = len(l_out)
    for spk,tpl in I.d['output'].items():
        if not tpl:
            continue
        for a in range(len(tpl[0])):
            d_out['speaker'].append(spk)
            for b in range(lo):
                d_out[l_out[b]].append(tpl[b][a])
        # Pandas export according to 'spk_csv' extension
    writeCSV(d_out,path,"wordAlign")
def aWord(I):
    """Aligns and restructures the Word tier ('wordAlign/restructWord')."""
    from reinj_code import wordAlign,restructWord

        # Get core tiers sorted by speaker and type
    I.d['d_typ'] = getNewTYP(I)

        # Get word-to-MAUS alignment
        ## Either from CSV/Excel file
    path = os.path.join(I.d['in_dir'],I.d['eaf'].name+"_"+I.d['wd_csv'])
    log(I.d['log_path'],"\t\t... aligning\n")
    if os.path.isfile(path):
        _readWord(I)
        ## Or from calling 'wordAlign()'
    else:
        I.d['output'] = {}
        g2p = wordAlign.get_g2p_dict(I.d['g2p'],I.d['lang'])
        for spk in I.d['d_core']:
            I.d['output'][spk] = wordAlign.wordAlign(I,spk,g2p)
            # Write 'output' to 'wd_csv'
        if I.d['wd_csv']:
            _writeWord(I)
        # Turn it into pointers
    log(I.d['log_path'],"\t\t... to pointers\n")
    for spk,tpl in I.d['output'].items():
        I.d['output'][spk] = wordAlign.getPntr(I,spk)
        # Restructuring
    log(I.d['log_path'],"\t\t... restructuring\n")
    restructWord.restructWord(I)
    return I
def _readMorph(I):
    """Sub-function to read CSV/Excel input for 'aMorph()'."""
    from reinj_code import morphAlign
    
        # Get path
    path = os.path.join(I.d['in_dir'],I.d['eaf'].name+"_"+I.d['mb_csv'])
    if not os.path.isfile(path):
        I.d['output'] = None
        return
    df = readCSV(path)     # pre-aligned words
        # We need a dict' for pandas' DataFrame
    d_spk = {}; l_rows = ['phone','morph_id','uncertainty']
    for a,row in df.iterrows():
        spk = str(row['speaker'])
        if not spk in d_spk:
            d_spk[spk] = ([],[],[])
        for b,k in enumerate(l_rows):
            val = row[k] if not pandas.isnull(row[k]) else ""
            d_spk[spk][b].append(str(val))
    I.d['output'] = d_spk
def _writeMorph(I):
    """Sub-function to write CSV/Excel output from 'aMorph()'."""
    
        # Get path
    path = os.path.join(I.d['out_dir'],I.d['eaf'].name+"_"+I.d['mb_csv'])
        # We need a dict' for pandas' DataFrame
    d_out = {'speaker':[],'phone':[],'morph_id':[],'uncertainty':[]}
        # We loop over 'output' to fill 'd_out'
    for spk,tpl in I.d['output'].items():
        if not tpl:
            continue
        for a in range(len(tpl[0])):
            d_out['speaker'].append(spk)
            d_out['phone'].append(tpl[0][a])
            d_out['morph_id'].append(tpl[1][a])
            d_out['uncertainty'].append(tpl[2][a])
        # Pandas export according to 'spk_csv' extension
    writeCSV(d_out,path,"morphAlign")
def aMorph(I):
    """Aligns and restructures the Morph tier ('morphAlign/restructMorph')."""
    from reinj_code import morphAlign,restructMorph
    
        # Get core tiers sorted by speaker and type
    I.d['d_typ'] = getNewTYP(I)  
        # Get morph-to-phone alignment
        ## Either from CSV/Excel file
    path = os.path.join(I.d['in_dir'],I.d['eaf'].name+"_"+I.d['mb_csv'])
    log(I.d['log_path'],"\t\t... aligning\n")
    if os.path.isfile(path):
        _readMorph(I)
    else:
        I.d['output'] = {}
        g2p = morphAlign.get_g2p_dict(I.d['g2p'],I.d['lang'])
        for spk in I.d['d_core']:
            I.d['output'][spk] = morphAlign.morphAlign(I,spk,g2p)
            # Write 'output' to 'wd_csv'
        if I.d['mb_csv']:
            _writeMorph(I)
        # Turn it into pointers
    log(I.d['log_path'],"\t\t... to pointers\n")
    for spk in I.d['d_core']:
        I.d['output'][spk] = morphAlign.getMorphPntr(I,spk)
        # Restructuring
    log(I.d['log_path'],"\t\t... restructuring\n")
    restructMorph.restructMorph(I)
    
    return I

    # Main function
def reinject(l_ops=['corr','merge','add','restruct','word','morph'],
             d_vals = {}):
    """Does everything?"""
    
    d_ops = {'corr':corrInput,
             'merge':mergeTGD,
             'add':addTiers,
             'restruct':restruct,
             'word':aWord,
             'morph':aMorph}
    def readOp(I):
        """Iterate over 'l_ops'."""
        for a,op in enumerate(l_ops):
            if op in d_ops:
                log(I.d['log_path'],"\t "+op+"...\n")
                I = d_ops[op](I)
    
        # Create and fill 'Inj'
    I = Inj(d_vals)
        # Import core tiers
        # We need that if at least for the 'lang' variable...
    if 'merge' in l_ops or 'add' in l_ops or 'restruct' in l_ops:
        I.d['d_tiers'] = fillSPK(I.d['tiers_csv'],l_path=I.d['log_path'])
        # Move to loop
    for file in os.listdir(I.d['in_dir']):
            # Skip non-ELAN files
        name,ext = os.path.splitext(file)
        if not ext.lower() == ".eaf":
            continue
            # Reset 'Inj'
        resetInj(I)
        log(I.d['log_path'],"Processing '"+file+"':\n")
        #try:
            # Imports files
        I.d['eaf'] = getEAF(I.d['in_dir'],file) # ELAN file
        for a in range(len(I.d['eaf'])-1,-1,-1):
            tier = I.d['eaf'].elem[a]
            if tier.name == "Pangloss":
                I.d['eaf'].pop(a)
        I.d['tgd'] = getTGD(I.d['in_dir'],name,
                            l_path=I.d['log_path']) # TextGrid file
            # Get 'd_core'
        I.d['d_core'],I.d['d_else'],I.d['lang'] = getSPK(I.d['d_tiers'],
                                                         I.d['eaf'],
                                                         I.d['l_ch'])
        l_lang = ['english','french']
        if I.d['lang']:
            if I.d['lang'].lower() in l_lang:
                I.d['weird'] = True
            I.d['eaf'].setMeta('lang',I.d['lang'],i=0)
        elif I.d['eaf'].checkMeta('lang'):
            I.d['lang'] = I.d['eaf'].meta('lang')
            # Operates
        readOp(I)
            # Exports an ELAN file
        log(I.d['log_path'],"\tExport.\n")
        exportEAF(I.d['eaf'],I.d['out_dir'],file)
        #except:
        #    log(I.d['log_path'],"\t-- FAILURE --\n")

if __name__ == '__main__':
    reinject()
    