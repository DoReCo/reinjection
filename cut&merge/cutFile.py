""""""

from corflow import fromPraat,fromElan,toPraat,toElan
import sys,os
home = os.path.abspath(os.path.dirname(__file__))

def iterBoth(indir,s,e=""):
    """Main loop iterator. Assumes a certain folder structure."""
    if not indir:
        indir = home
    if not e and s == "00_source":
        e = "10_cut"
    for d in os.listdir(indir):
        d_files = {}
        if not os.path.isdir(d):
            continue
        sdir = os.path.join(indir,d,s)
        if e:
            sout = os.path.join(indir,d,e)
            if not os.path.isdir(sout):
                os.mkdir(sout)
            for file in os.listdir(sdir):
                fi,ext = os.path.splitext(file)
                ftgd = os.path.join(sdir,fi+"_Algn"+".TextGrid")
                if (ext.lower() == ".eaf" and 
                   os.path.isfile(ftgd)):
                    d_files[fi] = {"EAF":os.path.join(sdir,file),
                                   "TGD":os.path.join(ftgd)}
        else:
            sout = os.path.join(indir,d)
            for file in os.listdir(sdir):
                fi,ext = os.path.splitext(file)
                if ext.lower() != ".eaf":
                    continue
                fin,n = fi.rsplit("_",1)
                if fin in d_files:
                    d_files[fin].insert(int(n)-1,os.path.join(sdir,file))
                else:
                    d_files[fin] = [os.path.join(sdir,file)]
        yield d,sout,d_files

def remSeg(seg):
    """Removes a segment and all of its children."""
    l_csegs = seg.allChildren()
    for cseg in l_csegs:
        cseg.struct.remove(cseg)
    seg.struct.remove(seg)
def clean(l_files,ch_t=300.):
    """Before merging, removes segments above 'ch_t' duration."""
    for path in l_files:
        trans = fromElan.fromElan(path)
        for tier in trans:
            if not tier.name.startswith("ref@") or not tier.elem:
                continue
            seg = tier.elem[-1]
            if seg.end-seg.start > ch_t:
                remSeg(seg)
            m = tier.elem[-1].end
        for tier in trans:
            for a in range(len(tier)-1,-1,-1):
                seg = tier.elem[a]
                if seg.end > m:
                    if seg.start >= m:
                        tier.remove(seg)
                    else:
                        seg.end = m
        trans.setBounds()
        toElan.toElan(path,trans)

def cutCut(tr,t1,t2,ch):
    """Sub-function of 'cut()' to remove segments before/after 't'.
    Note: 'ch' controls which part of the transcription to remove."""
    if t2 < 0.:
        t2 = tr.end
    for tier in tr:
        for a in range(len(tier)-1,-1,-1):
            seg = tier.elem[a]
            if seg.start >= t2:
                tier.pop(a)
            elif ch == ">" and seg.start >= t1:
                tier.pop(a)
            elif ch == "<" and seg.start < t1:
                tier.pop(a)
def cutCop(tr1,tr2,t1,t2):
    """Sub-function of 'cut()' to add segments to tr2 between 't1' and 't2'."""
    d_par = {}
    if t2 < 0:
        t2 = tr1.end
    for tier in tr1:
        pt = tier.parent().name if tier.parent() else ""
        ntier = tr2.add(-1,tier.copy(empty=True))
        d_par[ntier] = pt
        for seg in tier:
            if seg.start >= t1 and seg.start < t2:
                ntier.add(-1,seg)
        # Restructure
    s,e = -1.,-1.
    for tier in tr2:
        if tier.elem and (s < 0 or s > tier.elem[0].start):
            s = tier.elem[0].start
        if tier.elem and (e < 0 or e < tier.elem[-1].end):
            e = tier.elem[-1].end
        ptier = tr2.getName(d_par[tier])
        if ptier:
            tier.setParent(ptier)
            for seg in tier:
                mid = seg.start + ((seg.end-seg.start)/2)
                pseg = ptier.getTime(mid,ptier)
                seg.setParent(pseg)
    tr2.start = s; tr2.end = e
    for tier in tr2:
        tier.start = tr2.start; tier.end = tr2.end
    return tr2
def cut(path,outdir,l_t,fc1,fc2):
    """Splits a transcription into 2+ transcriptions by time code.
    ARGUMENTS:
    - 'path'        :   (str) Path to the input file.
    - 'outdir'      :   (str) Path to the output folder.
    - 'l_t'         :   (list) List of time codes for cuts.
    - 'fc1/fc2'     :   (pntr) Corflow import/export functions."""

    trans = fc1(path)
    p,file = os.path.split(path)
    fi,ext = os.path.splitext(file)
    cmpl = ""
    if ext.lower() == ".textgrid":
        fi = fi.rsplit("_",1)[0]; cmpl = "_Algn"
    incr = 2
    for a,t1 in enumerate(l_t):
        t2 = l_t[a+1] if a < len(l_t)-1 else -1.
        ntrans = trans.copy(empty=True)
        ntrans = cutCop(trans,ntrans,t1,t2)
        npath = os.path.join(outdir,fi+"_"+str(incr)+cmpl+ext); incr += 1
        fc2(npath,ntrans)
    if l_t:
        ntrans = trans.copy(empty=True)
        ntrans = cutCop(trans,ntrans,0.,l_t[0])
    else:
        ntrans = trans
    npath = os.path.join(outdir,fi+"_1"+cmpl+ext)
    fc2(npath,ntrans)
def cutEAF(path,outdir,l_t):
    """Calls 'cut()' with Elan functions."""
    cut(path,outdir,l_t,fromElan.fromElan,toElan.toElan)
def cutTGD(path,outdir,l_t):
    """Calls 'cut()' with Praat functions."""
    cut(path,outdir,l_t,fromPraat.fromPraat,toPraat.toPraat)
def mergeEAF(l_files,outdir,fi=""):
    """Merges 2+ EAF files, outputs EAF and TGD.
    ARGUMENTS:
    - 'l_files'     :   (list) List of EAF files.
    - 'outdir'      :   (str) Path to output folder.
    - 'fi'          :   (str) Name (no extension) of the output file."""
    
    def parse(tr1,tr2):
        d_tiers = {}                        # getName, basically
        m = -1.
        for tier in tr1:                    # tier access and max_time (m)
            d_tiers[tier.name] = tier
            if not tier.elem:
                continue
            seg = tier.elem[-1]
            if seg.end-seg.start > 300.:
                tier.remove(seg)
            if tier.elem[-1].end > m:
                m = tier.elem[-1].end
        n = -1.
        for tier in tr2:                 # Remove bunk segments
            if not tier.elem:
                continue
            seg = tier.elem[0]
            if seg.end-seg.start > 300.:
                tier.remove(seg)
            if not tier.elem:
                continue
            if tier.elem and (n < 0. or tier.elem[0].start < n):
                n = tier.elem[0].start
        return (d_tiers,m,n)
    def avoidOverlaps(tr1,tr2):
        d_tiers,m,n = parse(tr1,tr2)
        diff = m-n
        #if diff > 0.:
        for tier in tr2:            # increment time codes
            for seg in tier:
                seg.start = seg.start+diff
                seg.end = seg.end+diff
        for tier in tr2:            # remove previous segments
            for a in range(len(tier)-1,-1,-1):
                seg = tier.elem[a]
                if seg.start < m:
                    tier.remove(seg)
        return d_tiers
    def merge(tr1,tr2):
        d_tiers = avoidOverlaps(tr1,tr2)
        
        for tier in tr2:                    # Merge
            ntier = d_tiers[tier.name]
            for seg in tier:
                ntier.add(-1,seg)
        for tier in tr1:                    # Reparent
            ptier = tier.parent()
            if not ptier:
                for seg in tier:
                    seg.setParent(None)
                continue
            for seg in tier:
                mid = seg.start+((seg.end-seg.start)/2)
                pseg = ptier.getTime(mid,ptier)
                if not seg.parent() == pseg:
                    seg.setParent(pseg)
        rtier = tr1.findName("ref@")
        if rtier:
            incr = 1
            for seg in rtier:
                if not "_" in seg.content:
                    continue
                seg.content = ("{:04d}_"
                               .format(incr)+seg.content.split("_",1)[1])
                incr += 1
        tr1.setBounds()
        tr1.renameSegs()
    
    if not l_files:
        return
    if not fi:
        path = l_files[0]
        p,file = os.path.split(path)
        fi,ext = os.path.splitext(file)
        fi = fi.rsplit("_",1)[0]
    trans = fromElan.fromElan(l_files[0])
    for a in range(1,len(l_files)):
        ntrans = fromElan.fromElan(l_files[a])
        merge(trans,ntrans)
    toElan.toElan(os.path.join(outdir,fi+".eaf"),trans)
    toPraat.toPraat(os.path.join(outdir,fi+".TextGrid"),trans)

def main(mdir,indir,outdir,d_lang):
    """Either merges or cuts.
    ARGUMENTS:
    - 'mdir'        :   (str) Path to the main folder
    - 'indir'       :   (str) Path to the input sub-folder.
    - 'outdir'      :   (str) Path to the output sub-folder.
    - 'd_lang'      :   (dict) For each language
    Note: 'indir' must start with '00_' for cut, '20_' for merge
    Note: if 'mdir' is empty, defaults to script folder.
    Note: '00_source' indir generates a default '10_cut' outdir."""
    l_skip = ['Cashinahua','Popoluca','Sanzhi']
    for ndir,outdir,d_files in iterBoth(home,indir,outdir):
        #if ndir in l_skip:
        #    continue
        print(ndir)
        if indir.startswith("00_"):                   # cut
            for fi,d_typ in d_files.items():
                l_t1,l_t2 = d_lang.get(ndir)
                cutEAF(d_typ["EAF"],outdir,l_t1)      # cut EAF
                cutTGD(d_typ["TGD"],outdir,l_t2)      # cut TGD
        elif indir.startswith("20_"):
            for fi,l_files in d_files.items():
                print("\t",fi)
                #clean(l_files)                      # remove last pause...
                mergeEAF(l_files,outdir,fi)         # merge EAFs

if __name__ == "__main__":
    indir = "20_reinj" # "00_source" or "20_reinj"
    outdir = ""   # "10_cut" or ""
        # get TimeCodes (for EAF and TGD)
    d_lang = {'Cashinahua':([1773.266],[1776.256]),
              'Popoluca':([1983.172],[1983.312]),
              'Sanzhi':([1216.659],[]),
              'Yali':([1033.417,2121.164],[1034.869])}
    main(home,indir,outdir,d_lang)