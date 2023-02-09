

from corflow.Transcription import Corpus,Transcription
import os

mc = "mc-zero"
un = "doreco-mb-algn"
d_old = {}

    # Technical functions
def _escape(data):
    """Double-quotes for CSVs."""
    return "\""+data.replace("\"","\"\"")+"\""

def getSpk(l_trans):
    """Returns a dictionary of speaker-tiers.
    Note: naming conventions specific to DoReCo."""
    
    def fillDict(d_d,k,kk,v):
        if k in d_d:
            d_d[k][kk] = v
        else:
            d_d[k] = {kk:v}
    def getSingle(trans,d_trans):
        d_trans[trans.name] = {}
        for tier in trans:
            typ,spk = "",""
            if "@" in tier.name:
                typ,spk = tier.name.rsplit("@",1)
            elif "_" in tier.name:
                spk,typ = tier.name.split("_",1)
            if spk:
                tier.setMeta("typ",typ,"tech")
                fillDict(d_trans[trans.name],spk,typ,tier)
                for ntier in tier.allTree():
                    if ntier.name not in d_trans[trans.name][spk]:
                        d_trans[trans.name][spk][ntier.name] = tier
    
    d_trans = {}; ch_mc,ch_un = False,False     # fill 'd_spk'
    for trans in l_trans:
        getSingle(trans,d_trans)
    for name,d_spk in d_trans.items():          # Check for 'mc' and 'un'
        for spk,d_tiers in d_spk.items():
            wd_tier = d_tiers.get('wd'); mc_tier = d_tiers.get(mc)
            if wd_tier and mc_tier:
                ch_mc = True; wd_tier.setMeta("cmpl",mc_tier,"tech")
            elif wd_tier:
                wd_tier.setMeta("cmpl",None,"tech")
            mb_tier = d_tiers.get('mb'); un_tier = d_tiers.get(un)
            if mb_tier and un_tier:
                ch_un = True; mb_tier.setMeta("cmpl",un_tier,"tech")
            elif mb_tier:
                mb_tier.setMeta("cmpl",None,"tech")
        if ch_mc == True and ch_un == True:
            break
    return d_trans,ch_mc,ch_un
def checkTyp(ntier,ch_mc=False,ch_un=False):
    """Returns the complement for a given tier type ('ntier')."""
    if ntier == "wd" and ch_mc:
        return "id",mc
    elif ntier == "wd":
        return "id",""
    elif ntier == "mb" and ch_un:
        return "id",un
    elif ntier == "mb" or ntier == "ph":
        return "id",""
    return "",""
def writeHeader(f,l_tiers,ch_mc,ch_un,sep):
    """Writes the columns' names."""

    r_typ = l_tiers[0]
    id,cmpl = checkTyp(r_typ,ch_mc,ch_un)
    txt = ("{1}{0}{2}{0}{3}{0}{4}{0}{5}{0}{6}{0}{7}{0}{8}"
            .format(sep,"\"lang\"","\"file\"","\"core_extended\"",
                    "\"speaker\"",r_typ+"_ID",r_typ,"\"start\"","\"end\""))
    if cmpl:
        txt = txt+("{0}{1}".format(sep,cmpl))
    if len(l_tiers) < 2:
        f.write(txt+"\n"); return
    for typ in l_tiers[1:]:
        if not ch_mc and (typ == "refind" or typ == "isnref"): # non-mc corpus
            continue
        id,cmpl = checkTyp(typ,ch_mc,ch_un); id_typ = typ+"_ID"
        if (typ == "wd" and not ch_mc) or (typ == "mb" and not ch_un):
            txt = txt+("{0}{1}{0}{2}".format(sep,_escape(id_typ),_escape(typ)))
        elif (typ == "wd" and ch_mc) or (typ == "mb" and ch_un):
            txt = txt+("{0}{1}{0}{2}{0}{3}"
                       .format(sep,_escape(id_typ),_escape(typ),_escape(cmpl)))
        elif typ == "ph":
            txt = txt+("{0}{1}{0}{2}".format(sep,_escape(id_typ),
                       _escape(typ)))
        else:
            txt = txt+"{0}{1}".format(sep,_escape(typ))
    f.write(txt+"\n")
    return
def writeLines(f,spk,d_tiers,trans,l_tiers,ch_mc,ch_un,sep):
    """Writes the lines."""

    def writeLine(otxt,name,cont,nid,cmpl,ncmpl):
        if nid and not cmpl:
            otxt = otxt+("{0}{1}{0}{2}".format(sep,_escape(name),
                                               _escape(cont)))
        elif cmpl:
            otxt = otxt+("{0}{1}{0}{2}{0}{3}"
                         .format(sep,_escape(name),
                         _escape(cont),_escape(ncmpl)))
        else:
            otxt = otxt+("{0}{1}".format(sep,_escape(cont)))
        return otxt
    
    l_parse = []
    r_tier = d_tiers.get(l_tiers[0])
    if not r_tier:                  # No reference tier for that speaker
        return
    r_id,r_co = checkTyp(l_tiers[0],ch_mc,ch_un)
    if r_co:
        rc_tier = d_tiers.get(r_co,None)
    lang = trans.meta("lang")
    for rseg in r_tier:
        txt = ("{1}{0}{2}{0}{3}{0}{4}{0}{5}{0}{6}{0}{7:.03f}{0}{8:.03f}"
               .format(sep,"\""+lang+"\"", _escape(trans.name),
              trans.meta('core','tech'),"\""+spk+"\"",rseg.name.strip(),
              _escape(rseg.content),rseg.start,rseg.end))
        if r_co:
            cmpl = ""
            if rc_tier:
                cmpl = r_tier.getTime(rseg.start,rc_tier)
            cmpl = cmpl.content if cmpl else ""
            txt = txt+("{0}{1}"
                  .format(sep,_escape(cmpl)))
        if len(l_tiers) < 2:
            txt = txt+"\n"; continue
        for ntier in l_tiers[1:]:           # 'getTime()', no gap in DoReCo
            tier = d_tiers.get(ntier)
            nid,cmpl = checkTyp(ntier,ch_mc,ch_un)
            if not tier:                    # No such tier for that speaker
                txt = writeLine(txt,"","",nid,cmpl,""); continue
            cmpl_tier = d_tiers.get(cmpl,None); ncmpl = ""
            sn,si,seg = trans.getTime(rseg.start,tier,det=True)
            if ((tier.meta('type') in l_parse) and
                ((not tier in d_old) or (not d_old[tier] == seg))):
                d_old[tier] = seg
            elif tier.meta('type') in l_parse:
                txt = writeLine(txt,"","",nid,cmpl,""); continue
            cont = ""; name=""
            for a in range(si,len(tier)):   # List of segments to 'rseg'
                seg = tier.elem[a]
                if cmpl_tier:
                    cseg = cmpl_tier.getTime(seg.start,cmpl_tier)
                if seg.start >= rseg.end:
                    break
                elif seg.end < rseg.start:
                    continue
                cont = cont+" "+seg.content
                name = name+" "+seg.name
                if cmpl_tier and cseg:
                    ncmpl = ncmpl+" "+cseg.content
            name = name.strip()
            cont = cont.replace("\n",""); cont = cont.strip()
            ncmpl = ncmpl.replace("\n",""); ncmpl = ncmpl.strip()
            txt = writeLine(txt,name,cont,nid,cmpl,ncmpl)
        f.write(txt+"\n")
def saveTables(path,l_trans,sep,l_tiers,encoding):
    """Exports a set of Transcriptions into a tabular file.
    ARGUMENTS:
    - path          : (str) Full path to a directory or file.
    - trans         : (pntr) A Transcription instance.
    - sep           : (str) A CSV separator
    - encoding      : (str) The Elan file encoding.
    RETURNS:
    - Creates a tabular file at 'path' from 'l_trans'."""
    
        # Path
    if os.path.isdir(path):                     # If it's a directory
        path = os.path.join(path,trans.name+".csv") # Use 'trans.name'
    
    ch_mc,ch_un = False,False
    with open(path,'w',encoding=encoding) as f: # Open file
        d_trans,ch_mc,ch_un = getSpk(l_trans)
        writeHeader(f,l_tiers,ch_mc,ch_un,sep)
        for trans in l_trans:
            d_old = {}
            ntrans = trans.copy()               # We use a copy from there
            if not ntrans.meta('core','tech'):
                ntrans.setMeta('core','extended','tech')
            #for tier in ntrans:
            #    tier.sortByTime()
            d_spk = d_trans[ntrans.name]
            for spk,d_tiers in d_spk.items():
                txt = writeLines(f,spk,d_tiers,ntrans,l_tiers,ch_mc,ch_un,sep)
def _saveOne(path,trans,sep,l_tiers,encoding):
    """Exports a single Transcription into a tabular file."""
    saveTables(path,[trans],sep,l_tiers,encoding)

    # Main function
def toDTabular(path,trans,**args):
    """Exports one or more Tables.
    ARGUMENTS:
    - path          : (str) A full path to either a directory or a file.
    - trans         : (overloaded) A Transcription, Corpus or list of
                                   Transcriptions.
    - encoding      : (str) The file encoding.
    - sep           : (str) A CSV separator.
    - tiers         : (lst) A list of tiers for the crosstable.
    RETURNS:
    - Creates the Tables at 'path' from 'trans'.
    Note: Creates a copy for each Transcription while exporting.
    Note: limited to restructured DoReCo files."""
    
        # Args
    encoding = args.get('encoding','utf-8')# file encoding (for all files)
    sep = args.get('sep',",")           # A CSV separator
    l_tiers = args.get('tiers')         # The tiers to process
    if not l_tiers:
        raise("Need a 'tiers' parameter as list of tiers for crosstable.\n")
        # Overload
    f = d_load.get(type(trans))
    if f:
        f(path,trans,sep,l_tiers,encoding)
    else:
        raise KeyError("First argument must be of type 'Transcription/"+
                       "/Corpus/list'.")
d_load = {Transcription:_saveOne,Corpus:saveTables,
          list:saveTables}