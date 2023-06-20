""" 02/08/2022
After reinjection, code to export in different format (plus post-process).

The main functions are 'exp_1()' and 'exp_2()'.
- Call 'exp_1()' first on reinjected files.
- Then generate TEIs from <ct3.ortolang.fr/teiconvert>.
- Call 'exp_2()' on all resulting files.
Both functions take 'indir/outdir' (folder paths) as parameters.
'exp_1()' also takes 'glot' (file path) as a parameter.

'indir' is the input folder; 'outdir' is the output folder. 'glot' is the
'glottocodes.txt' file to generate a 'lang-to-glottocode' dictionary.

The main functions simply call the other functions in order.
'exp_1()' will do all corrections and generate the tables.
'ex_2()' will only rename files (TEI extensions, file prefixes).

"""

from corflow import fromPraat,fromElan,toPraat,toElan,toTEI
import toDTabular
import sys,os,re,time,shutil
import xml.etree.cElementTree as ETree

home = os.path.abspath(os.path.dirname(__file__))
def iter(indir,stop=".eaf"):
    for file in os.listdir(indir):
        fi,ext = os.path.splitext(file)
        path = os.path.join(indir,file)
        if stop and (not ext.lower() == stop):
            continue
        yield fi,ext,file,path
def iterAll(indir):
    for root,dirs,files in os.walk(indir):
        if root == indir:
            continue
        for file in files:
            fi,ext = os.path.splitext(file)
            path = os.path.join(root,file)
            yield fi,ext,file,path
def iterSingle(root,files):
    for file in files:
        fi,ext = os.path.splitext(file)
        path = os.path.join(root,file)
        yield (fi,ext,file,path)

def getGlot(p=""):
    if not p:
        p = os.path.join(home,"glottocodes.txt")
    d_glot = {}
    with open(p,'r',encoding="utf-8") as f:
        for line in f:
            line = line[:-1]
            name,code = line.split(",")
            d_glot[name] = code
    return d_glot

def cleanNegs(trans):
    for tier in trans:
        for a in range(len(tier)-1,-1,-1):
            if tier.elem[a].start < 0. or tier.elem[a].end < 0.:
                tier.allPop(a)

def chTr(trans):
    for tier in trans:
        for seg in tier:
            if seg.end <= seg.start:
                print("\t",trans.name,tier.name,seg.name,
                       seg.start,seg.end,seg.content)
def fixMC(trans):
    for tier in trans:
        if tier.name.startswith("mc-zero") and not "@" in tier.name:
            ptier = tier.parent()
            if not ptier:
                continue
            if not "@" in ptier.name:
                continue
            typ,spk = ptier.name.split("@",1)
            tier.name = tier.name+"@"+spk
def fixSort(trans):
    for tier in trans:
        tier.sortByTime()   
def checkLims(ptier,l_allChild,seg):
    for cseg in seg.allChildren():
        if cseg.end > seg.end:
            cseg.end = seg.end
def goshDarnPauses(trans):
    for tier in trans.getTop():
        for seg in tier:
            checkLims(tier,tier.allChildren(),seg)
        l_allChild = tier.allChildren()
        for a in range(len(tier)-1,0,-1):
            seg1,seg2 = tier.elem[a-1],tier.elem[a]
            nseg = None
            if seg1.end < seg2.start:
                nseg = tier.create(seg2.index(),"",
                                   seg1.end,seg2.start,"<p:>")
            if not nseg:
                continue
            for ctier in l_allChild:
                cseg = ctier.findTime(nseg.start)
                if not cseg:
                    cnseg = ctier.create(0,"",nseg.start,nseg.end,"<p:>")
                    cnseg.setParent(ctier); continue
                elif "<p:>" in cseg.content:
                    continue
                cnseg = ctier.create(cseg.index()+1,"",
                                     nseg.start,nseg.end,"<p:>")
                cnseg.setParent(nseg)
                if cseg.end > cnseg.start:
                    cseg.end = cnseg.start
    trans.renameSegs()
def alsoGaps(trans):
    
        def gapPar(tier,l_segs=[],start=-1.,end=-1.):
            if tier.name.startswith("ph@"):
                return []
            if not l_segs:
                l_segs = tier.elem
            l_res,s0,s1 = [],l_segs[0],l_segs[-1]
            if s1.end < end:  # end
                ncseg = tier.create(s1.index()+1,"",s1.end,end,"<p:>")
                l_res.append(ncseg)
            for a in range(len(l_segs)-1,0,-1):   # middle
                seg1,seg2 = l_segs[a-1],l_segs[a]
                if seg1.end < seg2.start:
                    ncseg = tier.create(seg2.index(),"",seg1.end,seg2.start,
                                        "<p:>")
                    l_res.append(ncseg)
            if s0.start > start:  # start
                ncseg = tier.create(s0.index(),"",start,s0.start,"<p:>")
                l_res.append(ncseg)
            #l_res.reverse()
            return l_res
        def findContext(tier,ind,pseg):
            prev = (-1,None); l_segs = []; next = (-1,None)
            for a in range(ind,-1,-1):
                seg = tier.elem[a]
                if seg.end-0.0001 <= pseg.start:
                    prev = (a,seg); ind = a; break
                elif seg.start-0.0001 >= pseg.end:
                    next = (a,seg)
                else:
                    l_segs.append((a,seg))
            return (ind,prev,l_segs,next)
        def gapChild(ctier,i,pseg,prev,l_segs,nex):
            if ctier.name.startswith("ph@"):
                return
            ccont = "****" if not "<p:>" in pseg.content else pseg.content
            if not l_segs:
                cseg = ctier.create(i+1,"",pseg.start,pseg.end,ccont)
                cseg.setParent(pseg)
            else:
                l_res = gapPar(ctier,l_segs,pseg.start,pseg.end)
                for cseg in l_res:
                    cseg.content = ccont
                    cseg.setParent(pseg)

        l_par,l_tmp = [],[]
        for tier in trans:      # parent tiers
            if (tier.elem and (not tier.parent()) and not 
                tier.name.startswith("doreco-mb-align")):
                gapPar(tier,tier.elem,trans.start,trans.end)
                l_par.append(tier)
            elif tier.name.startswith("doreco-mb-align"):
                for a in range(len(tier)-1,-1,-1):
                    seg = tier.elem[a]
                    if "****" in seg.content or "<p:>" in seg.content:
                        tier.pop(a)
            else:
                tier.setMeta("i",len(tier)-1,"tech")
        while l_par:            # child tiers
            for ptier in l_par:
                l_child = ptier.children()
                if not l_child or ptier.name.startswith("doreco-mb-align"):
                    continue
                for a in range(len(ptier)-1,-1,-1):
                    pseg = ptier.elem[a]
                    for ctier in l_child:
                        i = ctier.meta("i","tech")
                        i,prev,l_segs,nex = findContext(ctier,i,pseg)
                        for a in range(len(l_segs)):
                            l_segs[a] = l_segs[a][1]
                        l_segs.reverse()
                        ctier.setMeta("i",i,"tech")
                        gapChild(ctier,i,pseg,prev,l_segs,nex)
                    # append
                for ctier in l_child:
                    ctier.sortByTime()
                    for cseg in ctier:
                        mid = cseg.start+((cseg.end-cseg.start)/2)
                        cseg.setParent(ptier.getTime(mid,ptier))
                    if ctier.elem and ctier.children():
                        l_tmp.append(ctier)
            l_par = l_tmp.copy(); l_tmp.clear()
        trans.renameSegs()
  
def setupWipAlign(g2p,tol):
    from reinj_code import morphAlign as mbA
    from string import punctuation
    import pandas
    pandas.options.mode.chained_assignment = None
    
    def fillDict(kk,k,v,l_check,d_d):
        if not k in l_check:
            return
        if kk in d_d:
            d_d[kk][k] = v
        else:
            d_d[kk] = {k:v}
    def setChTime(pseg,l_l):
        ln,s,e,ch_seg = len(l_l),-1.,-1.,True
        if isinstance(pseg,tuple):
            s,e = pseg[0],pseg[1]; ch_seg = False
        else:
            s,e = pseg.start,pseg.end
        dur = e-s
        for a,cseg in enumerate(l_l):
            if ch_seg:
                cseg.setParent(pseg)
            cseg.start = s+((a*dur)/ln)
            cseg.end = s+(((a+1)*dur)/ln)
    def setMoarTime(segm,ch_moar=True,stop=[]):
        l_tmp,l_ctiers = [],segm.struct.children()
        d_csegs = segm.allChildDict()
        while l_ctiers:
            for ctier in l_ctiers:
                if ctier in stop: # stop part
                    continue
                l_csegs = d_csegs.get(ctier,[])
                if ch_moar and l_csegs and l_csegs[0].start >= 0.: # ch part
                    continue
                segm._split(l_csegs,segm)
                l_tmp = l_tmp+ctier.children()
            l_ctiers = l_tmp.copy(); l_tmp.clear()
    def remCh(oseg,stop=[]):
        d_csegs = oseg.allChildDict()
        l_tmp,l_ctiers = [],oseg.struct.children()
        while l_ctiers:
            for ctier in l_ctiers:
                if ctier in stop:
                    continue
                l_csegs = d_csegs.get(ctier,[])
                for a in range(len(l_csegs)-1,-1,-1):
                    ctier.remove(l_csegs[a])
                l_tmp = l_tmp+ctier.children()
            l_ctiers = l_tmp.copy(); l_tmp.clear()
        oseg.struct.remove(oseg)
    def cleanWD(wd_tier,mb_tier,tol=0.0):
        for a in range(len(wd_tier)-1,-1,-1):
            wseg = wd_tier.elem[a]
            if wseg.start == wseg.end:
                wd_tier.pop(a)
        for mseg in mb_tier:
            pseg = wd_tier.getTime(mseg.start+tol,wd_tier)
            mseg.setParent(pseg)
    def copSeg(oseg,s,e,stop=[]):
        nseg = oseg.struct.create(oseg.index(),"",s,e,oseg.content)
        for ctier,l_csegs in oseg.allChildDict().items():
            if not l_csegs or ctier in stop:
                continue
            ptier,l_ncsegs = ctier.parent(),[]
            i = l_csegs[0].index()
            for a in range(len(l_csegs)-1,-1,-1):
                cseg = l_csegs[a]
                l_ncsegs.append(ctier.create(i,"",-1.,-1.,cseg.content))
            l_ncsegs.reverse()
            setChTime((s,e),l_ncsegs)
            for ncseg in l_ncsegs:
                ncseg.setParent(ptier.getTime(ncseg.start,ptier))
        ptier = oseg.struct.parent()
        nseg.setParent(ptier.getTime(nseg.start,ptier))
    def addSeg(atier,s,e):
        l_seg = []
        aseg = atier.getTime(s+tol,atier)
        if not aseg:
            return []
        la = len(atier)
        while aseg.start < (e-tol):
            if not ("<<wip>>" in aseg.content):
                l_seg.append(aseg)
            i = aseg.index()+1
            if i >= la:
                break
            aseg = atier.elem[i]
        return l_seg   
    def fixTx(seg,ph_tier):
        pseg = seg.struct.getTime(seg.start,seg.struct.parent())
        d_csegs = pseg.childDict()
        check = True
        for ctier,l_csegs in d_csegs.items():
            if ctier == seg.struct:
                continue
            elif len(l_csegs) > 1:
                check = False; break
            elif not "****" in l_csegs[0].content:
                check = False; break
        if check:
            aseg = None
            if pseg.end > seg.end:
                aseg = pseg.struct.elem[pseg.index()-1]
                if aseg:
                    aseg.end = pseg.end
            elif pseg.start < seg.start:
                aseg = pseg.struct.elem[pseg.index()+1]
                if aseg:
                    aseg.start = pseg.start
            if aseg:
                remCh(pseg,stop=[seg.struct,ph_tier])
                setMoarTime(aseg,False,stop=[seg.struct,ph_tier])
    def processWIP(seg,g2p,punct,wd_tier,mb_tier,ph_tier):
        """Realigns for a given 'wip' segment."""
        i = seg.index()
        prev,nex = wd_tier.elem[i-1],wd_tier.elem[i+1]
        while prev.content.startswith("<"):
            prev = wd_tier.elem[prev.index()-1]
        while nex.content.startswith("<"):
            nex = wd_tier.elem[nex.index()+1]
            # Get me the 'mb' and 'ph' for those words
            ## /!\ Doesn't include <<wip>>
        l_phchild = addSeg(ph_tier,prev.start,nex.end)
        l_mbchild = addSeg(mb_tier,prev.start,nex.end)
            # Get the actual input for alignment
        morph_d = {k:[] for k in ['m_format', 'o_format','maus_align',
                              'orig_align','maus_wd', 'orig_wd',
                              'orig_mb', 'scope', 'num_morph', 'morph_issue',
                              'maus_gap', 'orig_gap', 'repl', 'num_align',
                              'uncertain']}
        wcont = prev.content+" "+nex.content
        phcont = ""
        for ph_seg in l_phchild:
            phcont = phcont+" "+ph_seg.content
        phcont = phcont.strip()
        l_mb,l_mbsamp = [],[]
        for a in range(len(l_mbchild)-1,-1,-1):
            mb_seg = l_mbchild[a]
            if mb_seg.content == "****":
                remCh(mb_seg,[ph_tier])
                l_mbchild.pop(a); continue
            if (mb_seg.content and not mb_seg.content[0] == '<'):
                cl_morph = mbA.clean_str(mb_seg.content.lower(), punct)
                cl_morph = mbA.parse_orthography(cl_morph, g2p, 
                                                 sep_char=' ')[1]
            else:
                cl_morph = ''
            l_mbsamp.insert(0,cl_morph)
            l_mb.insert(0,mb_seg.content)
            ## morphAlign
            # single morpheme (?)
        if len(l_mbchild) < 2:
            new_d = mbA.fill_skipped_morphs(phcont,wcont,l_mbsamp,l_mb,None)
            morph_d = mbA.update_morph_to_phone_dict(morph_d, new_d)
        else:
            new_d = mbA.align_maus_to_morph(phcont,wcont,l_mbsamp,l_mb,None)
            morph_d = mbA.update_morph_to_phone_dict(morph_d, new_d)
        morph_align_df = pandas.DataFrame(morph_d)
        phone_df = mbA.create_morph_phone_id_df(morph_align_df)
        l_ph = phone_df['phone'].tolist()
        l_in = phone_df['morph_id'].tolist()
        a = 0; omb = None; l_elim = []
        for b in range(len(l_ph)): # realign
            ph = l_ph[b]
            if ph == "€":
                continue
            elif (ph.startswith("<") and not "<<wip>>" in ph):
                a += 1; continue
            ph_seg = l_phchild[a]
            mb_seg = l_mbchild[l_in[b]]
            if not omb == mb_seg:
                mb_seg.start = ph_seg.start
                omb = mb_seg
                l_elim.append(omb)
            mb_seg.end = ph_seg.end
            a += 1
        cut_mb = (None,-1.,-1.) # fix intra-morph pause (1)
        for mb_seg in l_mbchild:
            if mb_seg not in l_elim:
                remCh(mb_seg,[ph_tier])
                mb_tier.remove(mb_seg); continue
            if mb_seg.start < seg.start and mb_seg.end > seg.end:
                cut_mb = (mb_seg,mb_seg.start,prev.end)
                mb_seg.start = seg.end
            for ctier,l_csegs in mb_seg.allChildDict().items():
                if ctier == ph_tier:
                    continue
                setChTime(mb_seg,l_csegs)
            mb_seg.setParent(wd_tier.getTime(mb_seg.start,wd_tier))
        if cut_mb[0]: # fix intra-morph pause (2)
            copSeg(*cut_mb,stop=[ph_tier])
        new = prev
        for ph_seg in l_phchild: # re-parent phones
            mid = ph_seg.start+((ph_seg.end-ph_seg.start)/2)
            ph_seg.setParent(mb_tier.getTime(mid,mb_tier))
        fixTx(seg,ph_tier) # merge utterances
    
    punct = ''.join([x for x in punctuation if (x != "'" and 
                                                x != "`" and x!= ':')])+' '
    return (mbA,punct,pandas,fillDict,setChTime,setMoarTime,remCh,
            cleanWD,copSeg,addSeg,fixTx,processWIP)
def wipSingleAlign(root,files,g2p="g2p_sample_CV.xlsx",tol=0.0):
    mbA,punct,pandas,fillDict,setChTime,setMoarTime,remCh, \
    cleanWD,copSeg,addSeg,fixTx,processWIP = setupWipAlign(g2p,tol)
    
    lang = os.path.basename(root)
    if lang.lower() == "english" or root.lower() == "french":
        return
    for fi,ext,file,path in iterSingle(root,files):
        if not ext.lower() == ".eaf":
            continue
        trans = fromElan.fromElan(path)
        d_spk,d_g2p = {},{}
        for tier in trans: # only wd/mb/ph tiers (per speaker)
            if not "@" in tier.name:
                continue
            typ,spk = tier.name.split("@",1)
            fillDict(spk,typ,tier,['wd','mb','ph'],d_spk)
        ch = False
        for spk,d_tiers in d_spk.items(): # find wip contexts
            wd_tier,mb_tier = d_tiers.get('wd'),d_tiers.get('mb')
            ph_tier = d_tiers.get('ph'); ch = False
            if not (wd_tier and ph_tier):
                continue
            if mb_tier:
                cleanWD(wd_tier,mb_tier,tol)
            for seg in wd_tier:
                if "<<wip>>" in seg.content: # found a context
                    if not d_g2p:
                        d_g2p = mbA.get_g2p_dict(g2p,lang)
                    if not mb_tier:
                        fixTx(seg,ph_tier); continue
                    processWIP(seg,d_g2p.get(lang),punct,
                               wd_tier,mb_tier,ph_tier) # process
                    ch = True
        for tier in trans:
            tier.sortByTime()
        goshDarnPauses(trans)
        alsoGaps(trans)
        for tier in trans:
            ptier = tier.parent()
            if ptier:
                for seg in tier:
                    mid = seg.start+((seg.end-seg.start)/2)
                    seg.setParent(ptier.getTime(mid,ptier))
        npath = os.path.join(root,file)
        cleanNegs(trans)
        toElan.toElan(npath,trans)
        npath = os.path.join(root,fi+".TextGrid")
        toPraat.toPraat(npath,trans)
def wipAlign(indir="",outdir="",g2p="g2p_sample_CV.xlsx",tol=0.0):
    """Realigns WIP contexts using biopairwise."""
    
    mbA,punct,pandas,fillDict,setChTime,setMoarTime,remCh, \
    cleanWD,copSeg,addSeg,fixTx,processWIP = setupWipAlign(g2p,tol)
    
    
    for root,dirs,files in os.walk(indir): # languages
        if root == indir:
            continue
        lang = os.path.basename(root)
        if (root == indir or lang.lower() == "english" or
            root.lower() == "french"):
            continue
        ndir = os.path.join(outdir,lang)
        if not os.path.isdir(ndir):
                os.mkdir(ndir)
        print(lang)
        for file in files: # file per language
            fi,ext = os.path.splitext(file)
            if not ext.lower() == ".eaf":
                continue
            path = os.path.join(root,file)
            trans = fromElan.fromElan(path)
            d_spk,d_g2p = {},{}
            for tier in trans: # only wd/mb/ph tiers (per speaker)
                if not "@" in tier.name:
                    continue
                typ,spk = tier.name.split("@",1)
                fillDict(spk,typ,tier,['wd','mb','ph'],d_spk)
            ch = False
            for spk,d_tiers in d_spk.items(): # find wip contexts
                wd_tier,mb_tier = d_tiers.get('wd'),d_tiers.get('mb')
                ph_tier = d_tiers.get('ph'); ch = False
                if not (wd_tier and ph_tier):
                    continue
                if mb_tier:
                    cleanWD(wd_tier,mb_tier,tol)
                for seg in wd_tier:
                    if "<<wip>>" in seg.content: # found a context
                        if not d_g2p:
                            d_g2p = mbA.get_g2p_dict(g2p,lang)
                        if not mb_tier:
                            fixTx(seg,ph_tier); continue
                        processWIP(seg,d_g2p.get(lang),punct,
                                   wd_tier,mb_tier,ph_tier) # process
                        ch = True
            for tier in trans:
                tier.sortByTime()
            alsoGaps(trans)
            for tier in trans:
                ptier = tier.parent()
                if ptier:
                    for seg in tier:
                        mid = seg.start+((seg.end-seg.start)/2)
                        seg.setParent(ptier.getTime(mid,ptier))
            #npath = os.path.join(ndir,file)
            #toElan.toElan(npath,trans)
            npath = os.path.join(ndir,fi+".TextGrid")
            toPraat.toPraat(npath,trans)
def renameSingleGlot(glot,root,files):
    l_ext = [".eaf",".xml",".csv",".textgrid"]
    for fi,ext,file,path in iterSingle(root,files):
        if fi.startswith("doreco_"):
            continue
        if not ext.lower() in l_ext:
            os.remove(path); continue
        elif ext.lower() == ".csv":
            continue
        npath = os.path.join(root,"doreco_{}_{}".format(glot,file))
        if os.path.isfile(npath):
            os.remove(npath)
        os.rename(path,npath)
        if os.path.isfile(path):
            os.remove(path)
def renameGlot(indir="",d_glot={}):
    if not indir:
        indir = os.path.join(home,"input")
    l_ext = [".eaf",".xml",".csv",".textgrid"]
    if not d_glot:
        d_glot = getGlot()
    
    for root,dirs,files in os.walk(indir):
        if root == indir:
            continue
        lang = os.path.basename(root)
        glot = d_glot.get(lang)
        if not glot:
            print(lang,"\tuhm...")
            continue
        sys.stdout.write("\rRenameGlot: "+lang+","+glot)
        for file in files:
            fi,ext = os.path.splitext(file)
            if fi.startswith("doreco_"):
                continue
            path = os.path.join(root,file)
            if not ext.lower() in l_ext:
                os.remove(path); continue
            elif ext.lower() == ".csv":
                continue
            npath = os.path.join(root,"doreco_{}_{}".format(glot,file))
            if os.path.isfile(npath):
                os.remove(npath)
            os.rename(path,npath)
def renameSingleAll(root,files):
    for fi,ext,file,path in iterSingle(root,files):
        if fi.startswith("doreco_") and not ext.lower() == ".csv":
            _1,_2,nfi = fi.split("_",2)
        elif (not fi.startswith("doreco_")) and (ext.lower() == ".csv"):
            nfi = "doreco_"+fi
        else:
            nfi = fi
        npath = os.path.join(root,nfi+ext)
        if path == npath:
            continue
        elif os.path.isfile(npath):
            os.remove(npath)
            os.rename(path,npath)
        else:
            shutil.copyfile(path,npath)
        if os.path.isfile(path):
            os.remove(path)
def renameAll(indir="",outdir=""):
    if not indir:
        indir = os.path.join(home,"input")
    if not outdir:
        outdir = os.path.join(home,"output")
    
    for root,dirs,files in os.walk(indir):
        if root == indir:
            continue
        d = os.path.basename(root)
        nd = os.path.join(outdir,d)
        if not os.path.isdir(nd):
            os.mkdir(nd)
        sys.stdout.write("\rRenameAll: "+d)
        for file in files:
            fi,ext = os.path.splitext(file); nfi = fi
            path = os.path.join(root,file)
            if fi.startswith("doreco_") and not ext.lower() == ".csv":
                _1,_2,nfi = fi.split("_",2)
            npath = os.path.join(nd,nfi+ext)
            if path == npath:
                continue
            elif nd == root:
                if os.path.isfile(npath):
                    os.remove(npath)
                os.rename(path,npath)
            else:
                shutil.copyfile(path,npath)
            if os.path.isfile(path):
                os.remove(path)
def renameTEI(indir=""):
        # Paths
    if not indir:
        indir = os.path.join(home,"input")
        #indir = os.path.join("C:\\\\Users\\delafontainef\\Downloads")
        
        # Main loop
    for root,dirs,files in os.walk(indir):
        for file in files:
            fi,ext = os.path.splitext(file)
            if not ext.lower().endswith(".xml"):
                continue
            path = os.path.join(root,file)
            if ".eaf.tei_corpo" in fi:
                fi = fi.split(".eaf",1)[0]
            npath = os.path.join(root,fi+".xml")
            os.rename(path,npath)
def rn(tr):
    incr = 1
    for tier in tr:
        if tier.name.startswith("ref@"):
            for seg in tier:
                if "<" in seg.content:
                    continue
                seg.content = "{:04d}_{}".format(incr,tr.name)
                incr += 1
    if 'elan' in tr.metadata:
        if 'MEDIA_DESCRIPTOR' in tr.metadata['elan']:
            tr.metadata['elan'].pop('MEDIA_DESCRIPTOR')
        tr.metadata['elan']['MEDIA_DESCRIPTOR'] = ["MEDIA_URL"
            "=\"{}.wav\" MIME_TYPE=\"audio/x-wav\"".format(tr.name)]
        if 'PROPERTY' in tr.metadata['elan']:
            l_prop = tr.metadata['elan']['PROPERTY']
            for a in range(len(l_prop)-1,-1,-1):
                prop = l_prop[a]
                prop = prop.replace("\t","").replace("\n","")
                prop = prop.strip()
                if not prop:
                    l_prop.pop(a)
            if l_prop:
                tr.metadata['elan']['PROPERTY'] = l_prop
            else:
                tr.metadata['elan'].pop('PROPERTY')
def reSingleAll(root,files):
    for fi,ext,file,path in iterSingle(root,files):
        if ext.lower() == ".textgrid":
            trans = fromPraat.fromPraat(path)
            rn(trans)
            toPraat.toPraat(path,trans)
        elif ext.lower() == ".eaf":
            trans = fromElan.fromElan(path)
            rn(trans)
            toElan.toElan(path,trans)
def reAll(indir=""):
    if not indir:
        indir = os.path.join(home,"input")
    for fi,ext,file,path in iterAll(indir):
        if ext.lower() == ".textgrid":
            sys.stdout.write("\r"+"ReAll: "+file)
            trans = fromPraat.fromPraat(path)
            rn(trans)
            toPraat.toPraat(path,trans)
        elif ext.lower() == ".eaf":
            sys.stdout.write("\r"+"ReAll: "+file)
            trans = fromElan.fromElan(path)
            rn(trans)
            toElan.toElan(path,trans)
def setupFixWIP():
    def cleanProps(trans):
        """Tries to keep metadata as clean as possible."""
        if 'elan' in trans.metadata and 'PROPERTY' in trans.metadata['elan']:
            l_props = trans.metadata['elan']['PROPERTY']
            for a in range(len(l_props)-1,-1,-1):
                prop = l_props[a]
                prop = prop.replace("\\t","").replace("\\n","")
                prop = prop.strip()
                if not prop:
                    l_props.pop(a)
            if not l_props:
                trans.metadata['elan'].pop('PROPERTY')
            else:
                trans.metadata['elan']['PROPERTY'] = l_props
    def noOldWD(trans):
        """Double-checks to remove 'wd_old'."""
        for a in range(len(trans)-1,-1,-1):
            tier = trans.elem[a]
            if tier.name.startswith("wd_old@"):
                trans.elem.pop(a)
    def fixZero(trans):
        """Removes zero-duration segments."""
        for tier in trans:
            ptier = tier.parent()
            if not ptier:
                continue
            for seg in tier:
                mid = seg.start+((seg.end-seg.start)/2)
                pseg = ptier.getTime(mid,ptier)
                seg.setParent(pseg)
        for tier in trans:
            for a in range(len(tier)-1,-1,-1):
                seg = tier.elem[a]
                if seg.end-seg.start < 0.0001:
                    l_csegs = seg.allChildren()
                    for cseg in l_csegs:
                        cseg.struct.remove(cseg)
                    tier.remove(seg)
    def renameMC(trans):
        """Handles refind and insref tier names."""
        for tier in trans:
            if (tier.name.startswith("refind") or 
                tier.name.startswith("isnref")):
                ptier = tier.parent()
                if not ptier:
                    print(trans.name,tier.name,"No parent?")
                elif "@" in ptier.name:
                    typ,spk = ptier.name.split("@")
                    tier.name = tier.name+"@"+spk
    def fixRefind(trans):
        """Should add separators to 'refind' segments."""
        refind = trans.getName("refind")
        if not refind:
            return
        d_len = {}
        r = 4
        for seg in refind:
            if not seg.content or "<p:>" in seg.content:
                continue
            i,txt = 0,""
            for a,char in enumerate(seg.content):
                if a > 0 and a%r == 0:
                    txt = txt+";"
                txt = txt+char
            seg.content = txt.strip("; ")
    return (cleanProps,noOldWD,fixZero,renameMC,fixRefind)
def fixSingleWIP(root,files):
    cleanProps,noOldWD,fixZero,renameMC,fixRefind = setupFixWIP()
    def fillEmpty(trans):
        for tier in trans:
            for seg in tier:
                if not seg.content:
                    seg.content = "****"
    for fi,ext,file,path in iterSingle(root,files):
        if not ext.lower() == ".eaf":
            continue
        trans = fromElan.fromElan(path)
        cleanProps(trans)                   # Get rid of empty 'PROPERTY'
        noOldWD(trans)                      # Double-check for no 'wd_old'
        fixZero(trans)                      # remove 0-words
        alsoGaps(trans)                     # Fill gaps everywhere
        fixRefind(trans)                    # more refind
        renameMC(trans)                     # refind and insref
        fillEmpty(trans)                    # put fillers in empty segs
        npath = os.path.join(root,file)
        toElan.toElan(npath,trans)
        npath = os.path.join(root,fi+".TextGrid")
        toPraat.toPraat(npath,trans)
def fixWIP(indir="",outdir=""):
    """Tries and fills the gaps in 'ps' tiers."""
    cleanProps,noOldWD,fixZero,renameMC,fixRefind = setupFixWIP()
    
    if not indir:
        indir = os.path.join(home,"input")
    if not outdir:
        outdir = os.path.join(home,"input")
    
    for root,dirs,files in os.walk(indir):
        if root == indir:
            continue
        d = os.path.basename(root)
        ndir = os.path.join(outdir,d)
        print(d)
        if not os.path.isdir(ndir):
            os.mkdir(ndir)
        for file in files:
            fi,ext = os.path.splitext(file)
            path = os.path.join(root,file)
            if not ext.lower() == ".eaf":
                continue
            sys.stdout.write("\r\t"+fi)
            #print("\t"+fi)
            trans = fromElan.fromElan(path)
            cleanProps(trans)                   # Get rid of empty 'PROPERTY'
            noOldWD(trans)                      # Double-check for no 'wd_old'
            fixZero(trans)                      # remove 0-words
            alsoGaps(trans)                     # Fill gaps everywhere
            fixRefind(trans)                    # more refind
            renameMC(trans)                     # refind and insref
            npath = os.path.join(ndir,file)
            toElan.toElan(npath,trans)
            npath = os.path.join(ndir,fi+".TextGrid")
            toPraat.toPraat(npath,trans)
def fixPauses(trans):
    
    for tier in trans.getTop():
        l_allChild = tier.allChildren()
        for a in range(len(tier)-1,-1,-1):
            seg = tier.elem[a]
            if not "<p:>" in seg.content:
                checkLims(tier,l_allChild,seg)
                continue
            for ctier in l_allChild:
                cseg = ctier.getTime(seg.start)
                if (not cseg or cseg.start == cseg.end 
                    or "<p:>" in cseg.content):
                    continue
                nseg = ctier.create(cseg.index()+1,"",seg.start,seg.end,"<p:>")
                nseg.setParent(seg)
                cseg.end = seg.start
    trans.renameSegs()
def fixSinglePauses(root,files):
    for fi,ext,file,path in iterSingle(root,files):
        if not ext.lower() == ".eaf":
            continue
        trans = fromElan.fromElan(path)
        fixPauses(trans)
        npath = os.path.join(root,file)
        toElan.toElan(npath,trans)
def notSingleTEI(glot,root,files):
    lph_trans,lwd_trans = [],[]
    d_typs = {'wd':('w',1),'mb':('m',1),'ph':('p',1)}
    for fi,ext,file,path in iterSingle(root,files):
        if not ext.lower() == ".eaf":
            continue
        trans = fromElan.fromElan(path)
        fixMC(trans)                    # quick fix for 'mc-zero' tier
        trans.setMeta('lang',glot)      # simpler glot attribution
        ch_ph = False                   # check if 'ph' tier present
            # Fix IDs for CSV tables
        for tier in trans:
            if not "@" in tier.name:
                continue
            typ,spk = tier.name.split("@",1)
            if typ == "ph":
                trans.setMeta('core','core','tech')
            tpl = d_typs.get(typ)
            if not tpl:
                continue
            n,i = tpl
            for seg in tier:
                seg.name = "{}{:06d}".format(n,i); i += 1
            d_typs[typ] = (n,i)
        if not trans.meta('core','tech'):
            trans.setMeta('core','extended','tech')
        for tier in trans:
            if not "@" in tier.name:
                continue
            typ,spk = tier.name.split("@",1)
            if typ == "ph":
                ch_ph = True; break
        if ch_ph:
            lph_trans.append(trans)
        lwd_trans.append(trans)
        fixSort(trans)
        npath = os.path.join(root,fi+".TextGrid")
        toPraat.toPraat(npath,trans)    # export to Praat
        npath = os.path.join(root,fi+".xml")
        toTEI.toTEI(npath,trans)
    if lph_trans:                       # tabular ('ph' tier)
        npath = os.path.join(root,glot+"_ph.csv")
        print("\t",glot+"_ph")
        toDTabular.toDTabular(npath,lph_trans,
                                  tiers=['ph','ref','tx','ft','wd',
                                         'mb','ps','gl',
                                         'refind','isnref'])
    if lwd_trans:                       # tabular ('wd' tier)
        npath = os.path.join(root,glot+"_wd.csv")
        print("\t",glot+"_wd")
        toDTabular.toDTabular(npath,lwd_trans,
                                  tiers=['wd','ref','tx','ft',
                                         'mb','ps','gl','ph',
                                         'refind','isnref'])
def notTEI(indir="",outdir="",d_glot={}):

    def setGlot(tr):
        gl = tr.name
        if not "_" in gl:
            return
        gl = gl.split("_",2)[1]
        tr.setMeta('lang',gl)
        
        # Paths
    if not indir:
        indir = os.path.join(home,"input")
    if not outdir:
        outdir = os.path.join(home,"input")
    if not d_glot:
        d_glot = getGlot()
        # Main loop
    for root,dirs,files in os.walk(indir):
        if root == indir:
            continue
        dir = os.path.basename(root); glot = dir
        ndir = os.path.join(outdir,dir)
        if not os.path.isdir(ndir):
            os.mkdir(ndir)
        print(dir)                          # for each directory
        lph_trans,lwd_trans = [],[]
        d_typs = {'wd':('w',1),'mb':('m',1),'ph':('p',1)}
        for file in files:
            fi,ext = os.path.splitext(file)
            if glot == dir and "doreco_" in fi:
                glot = "doreco_"+fi.split("_",2)[1]
            if not ext.lower() == ".eaf":
                continue
            path = os.path.join(root,file)
            #print("\t",fi)                  # for each (eaf) file
            trans = fromElan.fromElan(path)
            fixMC(trans)                    # quick fix for 'mc-zero' tier
            setGlot(trans)                  # 'lang' to glottocode
            ch_ph = False                   # check if 'ph' tier present
                # Fix IDs for CSV tables
            for tier in trans:
                if not "@" in tier.name:
                    continue
                typ,spk = tier.name.split("@",1)
                if typ == "ph":
                    trans.setMeta('core','core','tech')
                tpl = d_typs.get(typ)
                if not tpl:
                    continue
                n,i = tpl
                for seg in tier:
                    seg.name = "{}{:06d}".format(n,i); i += 1
                d_typs[typ] = (n,i)
            if not trans.meta('core','tech'):
                trans.setMeta('core','extended','tech')
            for tier in trans:
                if not "@" in tier.name:
                    continue
                typ,spk = tier.name.split("@",1)
                if typ == "ph":
                    ch_ph = True; break
            if ch_ph:
                lph_trans.append(trans)
            lwd_trans.append(trans)
            fixSort(trans)
            npath = os.path.join(ndir,fi+".TextGrid")
            toPraat.toPraat(npath,trans)    # export to Praat
            npath = os.path.join(ndir,fi+".xml")
            toTEI.toTEI(npath,trans)
        if lph_trans:                       # tabular ('ph' tier)
            npath = os.path.join(ndir,glot+"_ph.csv")
            print("\t",glot+"_ph")
            toDTabular.toDTabular(npath,lph_trans,
                                      tiers=['ph','ref','tx','ft','wd',
                                             'mb','ps','gl',
                                             'refind','isnref'])
        if lwd_trans:                       # tabular ('wd' tier)
            npath = os.path.join(ndir,glot+"_wd.csv")
            print("\t",glot+"_wd")
            toDTabular.toDTabular(npath,lwd_trans,
                                      tiers=['wd','ref','tx','ft',
                                             'mb','ps','gl','ph',
                                             'refind','isnref'])
def justCheck(indir=""):
        # Paths
    if not indir:
        indir = os.path.join(home,"output")
        
        # Main loop
    for root,dirs,files in os.walk(indir):
        dir = os.path.basename(root)
        print(dir)
        for file in files:
            fi,ext = os.path.splitext(file)
            if not ext.lower().endswith(".eaf"):
                continue
            print(fi)
            path = os.path.join(root,file)
            trans = fromElan.fromElan(path)
            chTr(trans)

    # Deprecated functions
def fixWO(indir="",outdir=""):
        # Paths
    if not indir:
        indir = os.path.join(home,"input")
    if not outdir:
        outdir = os.path.join(home,"output")
    
    for file in os.listdir(indir):
        fi,ext = os.path.splitext(file)
        path = os.path.join(indir,file)
        if not ext.lower() == ".textgrid":
            continue
        print(fi)
        trans = fromPraat.fromPraat(path)
        d_refs = {}
        for tier in trans:                      # Find tiers (by spk)
            if not "_" in tier.name:
                continue
            typ,ref = tier.name.rsplit("_",1)
            if not ref in d_refs:
                d_refs[ref] = {typ:tier}
            else:
                d_refs[ref][typ] = tier
        for ref,d_tiers in d_refs.items():      # Replace tiers
            wo, wf = d_tiers['W_O'],d_tiers['W_F']
            i_wo = wo.index()
            ntier = wf.copy(); ntier.name = "W_O_"+ref
            trans.add(i_wo,ntier)
            trans.remove(wo)
        npath = os.path.join(outdir,file)
        toPraat.toPraat(npath,trans)
def fixEvenki(indir="",outdir=""):
        # Paths
    if not indir:
        indir = os.path.join(home,"input")
    if not outdir:
        outdir = os.path.join(home,"output")
    
    for fi,ext,file,path in iter(indir):
        print(fi)
        trans = fromElan.fromElan(path)
        npath = os.path.join(outdir,file)
        toElan.toElan(npath,trans)
def tokThatJoke():
    """Retokenize Sanzhi."""
    path = os.path.join(home,"input","Sanzhi_09_08_2012_RMDF_001.eaf")
    trans = fromElan.fromElan(path)
    tx_tier = trans.getName("tx@KZ"); tx_ind = -1
    wd_tier = trans.getName("wd@KZ")
    ph_tier = trans.getName("ph@KZ")
    for seg in tx_tier:
        ref = seg.parent()
        d_csegs = ref._childDict(ref.allChildren())
        if ph_tier in d_csegs.keys():
            continue
        m = seg.start; tx_ind = seg.index()
        break
    for a in range(len(wd_tier)-1,-1,-1):
        wd_seg = wd_tier.elem[a]
        if wd_seg.start >= m:
            wd_tier.remove(wd_seg)
    for a in range(tx_ind,len(tx_tier)):
        seg = tx_tier.elem[a]
        l_toks = seg.content.split(" ")
        for b in range(len(l_toks)-1,-1,-1):
            if not l_toks[b]:
                l_toks.pop(b)
        lt = len(l_toks); dur = ((seg.end-seg.start)/lt)
        for b in range(lt):
            cont = l_toks[b]
            start = seg.start+(dur*b)
            end = seg.start+(dur*(b+1))
            if end > trans.end:
                end = trans.end
            wd_tier.create(-1,"",start,end,cont,seg)
    trans.renameSegs()
    npath = os.path.join(home,"output","Sanzhi_09_08_2012_RMDF_001.eaf")
    toElan.toElan(npath,trans)   
def merge(f1="",f2="",outdir=""):
    """Merging two files..."""
    if not f1:
        f1 = os.path.join(home,"input","interview_IP_Ware_1.eaf")
    if not f2:
        f2 = os.path.join(home,"input","interview_IP_Ware_3.eaf")
    if not outdir:
        outdir = os.path.join(home,"output")
        
    trans1 = fromElan.fromElan(f1)
    trans2 = fromElan.fromElan(f2)
    
    d_tiers = {}                        # getName, basically
    m = -1.
    for tier in trans1:                 # tier access and max_time (m)
        d_tiers[tier.name] = tier
        if tier.elem and tier.elem[-1].end > m:
            m = tier.elem[-1].end
    
    for tier in trans2:                 # Remove bunk segments
        if not tier.elem:
            continue
        seg = tier.elem[0]
        if seg.end-seg.start > 300.:
            tier.remove(seg)
    n = -1.
    for tier in trans2:                 # Get min time for trans2
        if tier.elem and (n < 0. or tier.elem[0].start < n):
            n = tier.elem[0].start
    diff = m-n
    if diff > 0.:
        for tier in trans2:             # avoid overlap
            for seg in tier:
                seg.start = seg.start+diff
                seg.end = seg.end+diff
        for tier in trans2:
            for a in range(len(tier)-1,-1,-1):
                seg = tier.elem[a]
                if seg.start < m:
                    tier.remove(seg)
    for tier in trans2:                 # Merge
        ntier = d_tiers[tier.name]
        for seg in tier:
            ntier.add(-1,seg)
    for tier in trans1:                 # Reparent
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
    rtier = trans1.findName("ref@")
    if rtier:
        incr = 1
        for seg in rtier:
            if not "_" in seg.content:
                continue
            seg.content = "_{:04d}".format(incr)+seg.content.split("_",1)[1] 
            incr += 1
    trans1.setBounds()
    trans1.renameSegs()
    
    file = os.path.basename(f1)         # Export
    npath = os.path.join(outdir,file)
    fi,ext = os.path.splitext(file)
    toElan.toElan(npath,trans1)
    npath = os.path.join(outdir,fi+".TextGrid")
    toPraat.toPraat(npath,trans1)
def cut(f="",outdir="",t=1033.431):
    if not f:
        f = os.path.join(home,"input","interview_IP_Ware.eaf")
    if not outdir:
        outdir = os.path.join(home,"output")
    
    def exp(p,tr,ext):
        print(p,ext)
        if ext.lower()==".eaf":
            toElan.toElan(p,tr)
        elif ext.lower()==".textgrid":
            toPraat.toPraat(p,tr)
    def split(trans,file):
        ntrans = trans.copy()
        fi,ext = os.path.splitext(file)
        for tier in trans:
            for a in range(len(tier)-1,-1,-1):
                seg = tier.elem[a]
                if seg.start > t:
                    tier.elem.pop(a)
        trans.setBounds()
        for tier in ntrans:
            for a in range(len(tier)-1,-1,-1):
                seg = tier.elem[a]
                if seg.end <= t:
                    tier.remove(seg)
        npath = os.path.join(outdir,file)
        exp(npath,trans,ext)
        npath = os.path.join(outdir,fi+"_2"+ext)
        exp(npath,ntrans,ext)
    
    file = os.path.basename(f)
    indir = os.path.dirname(f)
    fi,ext = os.path.splitext(file)
    p_tgd = os.path.join(indir,fi+"_Algn.TextGrid")
    trans = fromElan.fromElan(f)
    split(trans,file)
    if os.path.isfile(p_tgd):
        tgd = fromPraat.fromPraat(p_tgd)
        split(tgd,fi+"_Algn.TextGrid")
def fixWare(f=""):
    if not f:
        f = os.path.join(home,"input","interview_IP_Ware.eaf")
    trans = fromElan.fromElan(f)
    print(trans.name)
    for tier in trans:
        for a in range(len(tier)-2,-1,-1):
            seg1,seg2 = tier.elem[a],tier.elem[a+1]
            if seg1.end > seg2.start:
                print(tier.name,seg1.start,seg1.end,seg1.content,
                      seg2.start,seg2.end,seg2.content)
                tier.remove(seg2)
    toElan.toElan(f,trans)
    path,file = os.path.dirname(f),os.path.basename(f)
    fi,ext = os.path.splitext(file)
    npath = os.path.join(path,fi+".TextGrid")
    toPraat.toPraat(npath,trans)
def deleteStuff(indir=""):
    if not indir:
        indir = os.path.join(home,"output")
    file = os.path.join(home,"65_should_be_deleted.tsv")
    
    for root, dirs,files in os.walk(indir):
        if root == indir:
            continue
        for fi,ext,file,path in iter(root,""):
            if (not ext.lower() == ".csv") and fi.startswith("doreco_"):
                os.remove(path)
    
    """d_files = {}
    with open(file,'r',encoding="utf-8") as f:
        for line in f:
            lang,fi,typ = line.split("\t")
            typ = typ.replace("\n","")
            if not lang in d_files:
                d_files[lang] = {}
            d_files[lang][fi] = typ
    
    for root,dirs,files in os.walk(indir):
        if root == indir:
            continue
        d = os.path.basename(root)
        d_fi = d_files.get(d)
        if not d_fi:
            continue
        for fi,ext,file,path in iter(root,""):
            if fi in d_fi:
                os.remove(path)
                print(d,fi,ext)"""
  
    #### MAIN FUNCTIONS ####
def doEverything(indir="",glot=""):
    def updateFiles(root):
        files = []                        # We changed the names
        for file in os.listdir(root):
            if os.path.isfile(os.path.join(root,file)):
                files.append(file)
        return files
    
    if not indir:
        indir = os.path.join(home,"input")
    if not glot:
        glot = os.path.join(home,"glottocodes.txt")
    d_glot = getGlot()
    for root,dirs,files in os.walk(indir):
        if root == indir:
            continue
        lang = os.path.basename(root)
        glot = d_glot.get(lang)
        print(lang,glot)
        print("\tfixPauses...")
        fixSinglePauses(root,files) # I mean...
        print("\twipAlign...")
        wipSingleAlign(root,files)
        print("\trenaming...")
        renameSingleGlot(glot,root,files) # Renaming.
        print("\treAll...")
        files = updateFiles(root)
        reSingleAll(root,files) # Not sure what that's for anymore...
        print("\tfixWIP...")
        fixSingleWIP(root,files) # Bulk of the work...
        print("\tnotTEI...")
        notSingleTEI(glot,root,files) # Convert formats.
        files = updateFiles(root)
        print("\tRenameTEI...")
        renameTEI(root) # Probably deprecated...
        print("\tRenameAll...")
        renameSingleAll(root,files) # More renaming...
def exp_1(indir="",outdir="",glot=""):
        # Paths
    if not indir:
        indir = os.path.join(home,"input")
    if not outdir:
        outdir = indir
    if not glot:
        glot = os.path.join(home,"glottocodes.txt")
        # Operations
    d_glot = getGlot(glot)              # get Glottocodes as a dict
    renameGlot(indir,d_glot)            # add file prefixes
    reAll(indir)                        # fix 'ref' and 'wav' urls
    fixWIP(indir,outdir)                # fix WIP and properties...
    notTEI(indir,outdir)                # generate tables
def exp_2(indir="",outdir=""):
        # Paths
    if not indir:
        indir = os.path.join(home,"input")
    if not outdir:
        outdir = os.path.join(home,"output")
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        # Operations
    renameTEI(indir)                    # fix TEI-Corpo extension
    renameAll(indir,outdir)             # remove file prefixes
def test(indir=""):
        # Paths
    if not indir:
        indir = os.path.join(home,"input")
    for fi,ext,file,path in iterAll(indir):
        if not ext.lower() == ".eaf":
            continue
        try:
            trans = fromElan.fromElan(path)
        except:
            print("Didn't open:",fi)
        for tier in trans:
            for seg in tier:
                seg0 = None
                if seg.index() > 0:
                    seg0 = tier.elem[seg.index()-1]
                if seg.start == seg.end:
                    print("0-dur",fi,tier.name,seg.start); break
  
if __name__ == "__main__":
    indir = os.path.join(home,"03_export")
    #outdir = os.path.join(home,"export")
    doEverything(indir)



    # Code for 'fixStructure()' in 'fromElan.py'
    ## (Evenki fix.)
"""d_refs = {}
    for name,tpl in d_segs.items():
        seg,ch,ref = tpl
        if ref in d_refs:
            d_refs[ref].append(seg)
        else:
            d_refs[ref] = [seg]
    
    for tier in trans:
        if not tier.name == "fon_old":
            continue
        if not tier.elem or tier.elem[0].start < 0.:
            continue
        mid = (len(tier)//2)
        for a in range(len(tier)-1,mid-1,-1):
            seg1,seg2 = tier.elem[a],tier.elem[a-mid]
            print(a,a-mid,seg1.start,seg2.start)
            if seg1.name in d_refs:
                for cseg in d_refs[seg1.name]:
                    d_segs[cseg.name] = (cseg,d_segs[cseg.name][1],seg2.name)
            tier.elem.pop(a)"""