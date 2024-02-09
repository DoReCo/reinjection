from corflow import fromPraat,fromElan,toPraat,toElan
import reinject
from reinj_code import wordAlign
from Bio import pairwise2
import os,shutil,re,time,pandas

    # Paths
home = os.path.abspath(os.path.dirname(__file__))
idir = os.path.join(home,"input")
odir = os.path.join(home,"output")

    # Technical functions
def iter(d,stop=".eaf"):     # file iterator (for 'os')
    for file in os.listdir(d):
        fi,ext = os.path.splitext(file)
        path = os.path.join(d,file)
        if not ext.lower() == stop:
            continue
        yield fi,ext,file,path
def getG2P(g2p,lang):           # load g2p as dict[text] = MAUS
    """Get the G2P for Beja."""
    g2p_dict = pandas.read_excel(g2p, sheet_name=lang, engine='openpyxl')
    g2p_dict["text"] = g2p_dict["text"].apply(str)
    g2p_dict["MAUS"] = g2p_dict["MAUS"].apply(str)
    d_g2p = {}
    for a in range(len(g2p_dict)):
        d_g2p[g2p_dict["text"][a]] = g2p_dict["MAUS"][a]
    return d_g2p
def parentTGD(trans):           # parenting tiers to 'W_O_'
    """Parent phones to words in TGD."""
    d_spk = {}
    for tier in trans:
        if tier.name.count("_") < 2:
            continue
        typ1,typ2,spk = tier.name.split("_",2); typ = typ1+"_"+typ2+"_"
        if not spk in d_spk:
            d_spk[spk] = {}
        d_spk[spk][typ] = tier
    for spk,d_tiers in d_spk.items():
        ptier = d_tiers["W_O_"]
        for typ,tier in d_tiers.items():
            if tier == ptier:
                continue
            tier.setParent(ptier)
            for seg in tier:
                seg.setParent(seg.getTime(seg.start,ptier))
def getSpk(eaf,tgd,d_tiers):    # Sort tiers by speaker (EAF and TGD)
    """Returns a dict' with for each speaker a pair of tiers."""
    
    def fillDict(d_d,k,v):
        if k in d_d:
            d_d[k][v.name] = v
        else:
            d_d[k] = {v.name:v}
    def findStuff(eaf,d_tiers):
        d_tmp,d_rtier = {},{}
        for spk,d_stuff in d_tiers.items():
            for tier,d_vals in d_stuff.items(): # Get ELAN tiers
                #print(tier.name,d_vals)
                if d_vals['type'] == "tx":
                    fillDict(d_rtier,spk,tier)
                elif d_vals['type'] == "wd":
                    fillDict(d_rtier,spk,tier)
                    d_tmp[spk] = tier
        return d_tmp,d_rtier
    
    d_spk = {}
    d_tmp,d_rtier = findStuff(eaf,d_tiers)
    if not d_rtier:
        print("ELAN file doesn't contain expected tiers:",d_spk,d_rtier)
        return d_spk
    for spk,d_rn in d_rtier.items():      # Add TGD tiers
        tier = None
        for rn,ctier in d_rn.items():
            tier = tgd.getName("W_O_"+rn)
            if tier:
                break
        if tier:
            wd = d_tmp.get(spk)
            d_spk[rn] = [wd,tier]
    return d_spk
def getIStuff(I):
    I.d['d_core'],I.d['d_else'],I.d['lang'] = reinject.getSPK(
                                                     I.d['d_tiers'],
                                                     I.d['eaf'],
                                                     I.d['l_ch'])
    if I.d['lang']:
        I.d['eaf'].setMeta('lang',I.d['lang'],i=0)
    elif I.d['eaf'].checkMeta('lang'):
        I.d['lang'] = I.d['eaf'].meta('lang')
    I = reinject.mergeTGD(I,edit=False)
    I = reinject.addTiers(I)
    I.d['d_typ'] = reinject.getTYP(I.d['d_core'],I.d['l_ch'])
def wdAlign(I):
    I.d['d_typ'] = reinject.getNewTYP(I)
    I.d['output'] = {}
    d_g2p = wordAlign.get_g2p_dict(I.d['g2p'],I.d['lang'])
    for spk in I.d['d_core']:
        I.d['output'][spk] = wordAlign.wordAlign(I,spk,d_g2p)
    for spk,tpl in I.d['output'].items():
        I.d['output'][spk] = wordAlign.getPntr(I,spk)
    for spk,tpl in I.d['output'].items():
        if not tpl:
            continue
        l_osegs,l_msegs = tpl; mtier,omtier = None,None; incr = 0
        for mn,mi,mseg in l_msegs:        # find the TGD maus word tier
            if mseg:
                mtier = mseg.struct
                omtier = mtier.meta("old_tier","tech"); break
        if (omtier and mtier.elem[0].content == "<p:>" and 
            ((not omtier.elem[0].content == "") or
             (not omtier.elem[0].content == "<p:>"))):
            incr = -1
        for a in range(len(l_osegs)):
            oseg,mseg = l_osegs[a][2],l_msegs[a][2]
            if oseg:
                oseg.setMeta("pair",mseg,"tech")
            if mseg:
                mseg.setMeta("pair",oseg,"tech")
                if omtier:
                    ind = mseg.index()+incr
                    if ind < 0 or ind >= len(omtier):
                        continue
                    omtier.elem[mseg.index()].setMeta("pair",oseg,"tech")
                    
    # Log data
def startLog(lpath):
    f = open(lpath,'w',encoding="utf-8"); f.close()
def log(text="",p=""):
    if not p:
        print(text)
    else:
        with open(p,'a',encoding="utf-8") as f:
            f.write(text)

    # Technical functions
def _cleanString(data):
        """Remove special symbols."""
        data = wordAlign.strip_maus(data,add_space=False)
        data = wordAlign.strip_text(data,wordAlign.punct)
        return data
def toStrings(l_tiers):                 # Turn segments into a single string
    """Returns two strings and char-to-seg for each tier."""
    
    def toStringForReal(tier):
        l_res = []; l_match = []
        tier.fixGaps("")
        for seg in tier:
            if "<p:>" in seg.content or not seg.content:
                continue
            cont = _cleanString(seg.content)
            for char in cont:
                l_res = l_res+[char]; l_match.append(seg)
        return l_res,l_match
    
    ll_info = []
    for tier in l_tiers:
        l_res,l_match = toStringForReal(tier)
        ll_info.append([l_res,l_match])
    return ll_info
def pairStr(l1,l2):                 # Bio.pairwise2
    """Run biopairwise2."""
        # Variables
    alignment = pairwise2.align.globalxs(l1,l2,-1,-1,
                              gap_char=["€"],one_alignment_only=True)
    return [alignment[0][0],alignment[0][1]]
def reassociate(ll_info,l_align):       # Assign characters to segments
    """Assign everything back and find segments to split that way."""
    
    def incrP(pos,data1,ln,data2):
        if not data1 == "€":
            pos = pos+1; data2 = data2+data1
        if pos >= ln:
            pos = pos-1
        return data2,pos
    def addToSplit(d_split,seg,data):
        if seg in d_split:
            d_split[seg].append(data)
        else:
            d_split[seg] = [data]
    def getMsegs(mstart,mend):
        """Gets 'l_msegs'."""
        l_msegs = []
        mn,mi,mseg = m_tier.getTime(mstart,m_tier,det=True)
        for i in range(mi,len(m_tier)):
            mseg = m_tier.elem[i]
            if mseg.start >= mend:
                break
            l_msegs.append(mseg)
        return l_msegs
    
        # Unpack
    o_eaf,o_meaf = ll_info[0][0],ll_info[0][1]; le = len(o_meaf)
    o_tgd,o_mtgd = ll_info[1][0],ll_info[1][1]; lt = len(o_mtgd)
    n_eaf,n_tgd = l_align[0],l_align[1]
        # Reassociate
    d_split = {}
    p_eaf,p_tgd = 0,0; ch_eaf,ch_tgd = None,None; t_eaf,t_tgd = "",""
    m_tier = None
    
    for a in range(len(n_eaf)):
            # Variables
        c_eaf,c_tgd = n_eaf[a],n_tgd[a]
        eaf_seg = o_meaf[p_eaf]; tgd_seg = o_mtgd[p_tgd]
        if tgd_seg and not m_tier:
            m_tier = tgd_seg.struct
            # Find splits
        if not ch_eaf:                              # Start
            ch_eaf = eaf_seg; ch_tgd = tgd_seg
        elif (not c_eaf == "€" 
              and not ch_eaf == eaf_seg):           # New EAF segment
            if (ch_tgd == tgd_seg) and t_tgd:       # Split
                if not len(eaf_seg.content) <= 1:   # 1-char check
                    addToSplit(d_split,tgd_seg,t_tgd)
            ch_eaf = eaf_seg; ch_tgd = tgd_seg
            t_eaf,t_tgd = "",""
        if (not c_tgd == "€" and ch_eaf 
            and not ch_tgd == tgd_seg):             # New TGD segment
            ch_tgd = tgd_seg; t_tgd = ""
            # Increment
        t_eaf,p_eaf = incrP(p_eaf,c_eaf,le,t_eaf)
        t_tgd,p_tgd = incrP(p_tgd,c_tgd,lt,t_tgd)
        # Ignore 1-char strings in 'd_split'
    for mseg,l_txt in d_split.items():
        for a in range(len(l_txt)-1,-1,-1):
            txt = l_txt[a]
            if len(txt) == 1:
                if a+1 < len(l_txt):
                    l_txt[a+1] = l_txt[a]+l_txt[a+1]
                l_txt.pop(a)
        cont = _cleanString(mseg.content)
        if l_txt and len(cont.lower().rsplit(l_txt[-1],1)[1]) == 1:
            l_txt = [mseg.content]
        d_split[mseg] = l_txt
    return d_split
def toG2P(d_g2p,data):                  # Convert string with g2p
    """Turns a string into a g2p-translated version."""
    
    def recht(ch):
        if ch:
            return " "
        return ""
    def singleTransl(res,incr,cht=False):
        nc = d_g2p.get(data[incr])        # A single character
        if nc:
            res = res+nc+recht(cht)
        else:
            res = res+data[incr]+recht(cht)
        return res
    def addTransl(res,incr,cht=False):
        nc = d_g2p.get(data[incr]+data[incr+1])
        if not nc:
            if data[incr] == data[incr+1]:
                nc = d_g2p.get(data[incr])
                if nc:
                    res = res+nc+"ː"+recht(cht)
                else:
                    res = res+data[incr]+"ː"+recht(cht)
                incr += 1
            elif data[incr+1] == "ː":
                res = res+data[incr]+data[incr+1]+recht(cht); incr += 1
            else:
                res = (res+singleTransl("",incr,cht))
        else:
            res = res+nc+recht(cht); incr += 1
        return incr,res
    
    res1,res2 = "",""; a = 0; ld = len(data)
    while a < ld:
        if a+1 >= ld:
            res1 = singleTransl(res1,a)
            res2 = singleTransl(res2,a,True); break
        _,res1 = addTransl(res1,a)
        a,res2 = addTransl(res2,a,True)
        a += 1
    return res1,res2
def timeToSplit(tgd,d_g2p,spk,m_tier,d_split):  # everything else
    """Split words, split children, attribute phonemes."""
    
    def findTime(l_phsegs,a,data):
        """Find the correct time boundary from the phone tier."""
        count = len(data.split(" "))-2+a
        if count >= len(l_phsegs):
            count = len(l_phsegs)-1
        return count,l_phsegs[count].end
    def splitMaus(ms_tier,tph,s_txt):
        """Splits the MAUS sampa tier."""
        nc,ic,ms_seg = ms_tier.getTime(tph,ms_tier,det=True)
        count = len(s_txt.split(" "))-2
        l_cont = ms_seg.content.split(" ")
        cont1,cont2 = "",""
        for a in range(0,count+1):
            cont1 = cont1+l_cont[a]+" "
        for a in range(count+1,len(l_cont)):
            cont2 = cont2+l_cont[a]+" "
        cont1.strip(); cont2.strip()
        nseg = ms_tier.create(ic,"",ms_seg.start,tph,cont1)
        nseg.setMeta("pair",ms_seg.meta("pair","tech"),"tech")
        ms_seg.start = tph; ms_seg.content = cont2
    def splitElse(mseg,d_children,ph_tier,ms_tier,tph,o_txt):
        """Splits every other tier, MAUS word included."""
        m_tier = mseg.struct; im = mseg.index()
        count = len(o_txt)
        cont1,cont2 = mseg.content[:count],mseg.content[count:]
        nseg = m_tier.create(im,"",mseg.start,tph,cont1)
        nseg.setMeta("pair",mseg.meta("pair","tech"),"tech")
        mseg.start = tph; mseg.content = cont2
        for ctier in d_children:
            if ctier == ph_tier or ctier == ms_tier:
                continue
            nc,ic,cseg = ctier.getTime(tph,ctier,det=True)
            nseg = ctier.create(ic,"",cseg.start,tph,cont1)
            nseg.setMeta("pair",cseg.meta("pair","tech"),"tech")
            cseg.start = tph; cseg.content = cont2
    
    ph_tier = tgd.getName("P_S_"+spk)
    ms_tier = tgd.getName("W_S_"+spk)
    for mseg,l_txt in d_split.items():  # For each MAUS word seg' to split
        d_children = mseg._childDict(mseg.children())
        pos = 0; ntxt = ""
        for txt in l_txt:               # For each part
            ntxt = ntxt+txt
            if len(ntxt) >= len(mseg.content):
                break
            ch = mseg.content[len(ntxt):]
            if not ch:
                continue
            o_txt,s_txt = toG2P(d_g2p,ntxt)
            s_txt = s_txt.replace("  "," ") # 09.02.2024 (remove empty spaces)
            i,tph = findTime(d_children[ph_tier],pos,s_txt)
            pos = i+1
            log("\tSPLIT: ({:.02f},{:.02f}); {}; txt: {}, {}; "
                "o_txt: {}; s_txt: {}; tph: {}, ch: {}"
                .format(mseg.start,mseg.end,mseg.content,txt,
                        ntxt,o_txt,s_txt,tph,ch),glog_path)
            if mseg.end == tph:
                log(";\n",glog_path); continue
            else:
                log("; split\n",glog_path)
            splitMaus(ms_tier,tph,s_txt)
            splitElse(mseg,d_children,ph_tier,ms_tier,tph,ntxt)
def reparent(tgd):
    """Fix whatever timeToSplit() did."""
    for tier in tgd:
        ptier = tier.parent()
        if not ptier:
            continue
        for seg in tier:
            pseg = tier.getTime(seg.start,ptier)
            seg.setParent(pseg)
def p2(s1, s2):
    s1 = wordAlign.strip_text(s1,wordAlign.punct)
    s2 = wordAlign.strip_text(s2,wordAlign.punct)
    if isinstance(s2, list):
        s2 = ''.join(s2)
    al = pairwise2.align.globalms(list(s1), list(s2),
                                  1, -1, -1, -1, gap_char=['€'],
                                  one_alignment_only=True)
    return (al[0][2], al[0][2]/len(al[0][0]))
def timeToMerge(tgd,spk,m_tier):
    """Fuse stuff, look I'm... ugh."""
    
    ph_tier = tgd.getName("P_S_"+spk)
    ms_tier = tgd.getName("W_S_"+spk)
    
    def fuse(seg1,seg2,mend):
        """Actual fusion."""
        seg1.content = seg1.content+seg2.content
        seg1.end = mend
        seg2.struct.remove(seg2)
    def fuseChild(seg1,seg2,mend):
        """Fuse children."""
        d_children1 = seg1._childDict(seg1.children())
        d_children2 = seg2._childDict(seg2.children())
        for ctier,l_csegs2 in d_children2.items():  # Fuse children
            if ctier == ph_tier:                    # phonemes
                for cseg in l_csegs2:
                    cseg.setParent(seg1)
            else:
                l_csegs1 = d_children1[ctier]       # assumed always true
                if (ctier == ms_tier and not
                    l_csegs1[0].content.endswith(" ")):
                    l_csegs1[0].content = l_csegs1[0].content+" "
                fuse(l_csegs1[0],l_csegs2[0],mend)

    a = len(m_tier)-1; ch_while = True; count = 0
    while ch_while:
        ch_while = False; count += 1
        for a in range(len(m_tier)-1,-1,-1):
            mseg = m_tier.elem[a]; check = mseg.meta("pair","tech")
            if check or not mseg.content or "<<" in mseg.content:
                continue
                # Only MAUS words with no ORIG words associated left
            lmaus = m_tier.elem[a-1] if a-1 >= 0 else None
            rmaus = m_tier.elem[a+1] if a+1 < len(m_tier) else None
            lorig = lmaus.meta("pair","tech") if lmaus else None
            rorig = rmaus.meta("pair","tech") if rmaus else None
            left_m = lmaus.content+mseg.content if lmaus else ""
            right_m = mseg.content+rmaus.content if rmaus else ""
            score_dict = {'left':-1000, 'right':-1000,
                          'left-merge':-10000, 'right-merge':-10000}
                #Set scores
            if lorig and lmaus.content and not "<<" in lmaus.content:
                score_dict['left'] = p2(lmaus.content, lorig.content)[1]
                score_dict['left-merge'] = p2(left_m, lorig.content)[1]
            if rorig and rmaus.content and not "<<" in rmaus.content:
                score_dict['right'] = p2(rmaus.content, rorig.content)[1]
                score_dict['right-merge'] = p2(right_m, rorig.content)[1]
                # Compare & merge
            if (score_dict['left-merge'] > score_dict['left'] and
                score_dict['right-merge'] <= score_dict['right']):
                log("\tMERGE_{}: ({:0.2f}s,{:0.2f}s); '{}'; score: {{"
                .format(count,mseg.start,mseg.end,mseg.content,score_dict),
                    glog_path)
                for k,v in score_dict.items():
                    log("{}: {:.04f},".format(k,v),glog_path)
                log("}}; left_content: "+left_m+"\n",glog_path)
                fuseChild(lmaus,mseg,mseg.end); fuse(lmaus,mseg,mseg.end)
                ch_while = True
            elif (score_dict['right-merge'] > score_dict['right'] and
              score_dict['left-merge'] <= score_dict['left']):
                log("\tMERGE_{}: ({:0.2f}s,{:0.2f}s); '{}'; score: {{"
                .format(count,mseg.start,mseg.end,mseg.content,score_dict),
                glog_path)
                for k,v in score_dict.items():
                    log("{}: {:.04f},".format(k,v),glog_path)
                log("}}; right_content: "+right_m+"\n",glog_path)
                mseg.setMeta("pair",rmaus.meta("pair","tech"),"tech")
                fuseChild(mseg,rmaus,rmaus.end); fuse(mseg,rmaus,rmaus.end)
                ch_while = True
def splitnMerge(eaf,tgd,d_tiers,d_g2p=None,op=""):
    """Compares word tiers and splits TGD one if need be.
    Note : /!\ 'spk' is the 'tx' tier name."""
    
    parentTGD(tgd)
    d_spk = getSpk(eaf,tgd,d_tiers)                     # Sort by speaker
    for spk,l_tiers in d_spk.items():                   # For each speaker...
        log("\t"+tgd.name+", "+spk+"\n",glog_path)
        start = time.time()
        if not op or op=="split":
            ll_info = toStrings(l_tiers)                    # turn to strings
            l_align = pairStr(ll_info[0][0], ll_info[1][0]) # Bio.pairwise2
            d_split = reassociate(ll_info,l_align)          # Find segs'
            timeToSplit(tgd,d_g2p,spk,l_tiers[1],d_split)   # Split them
            reparent(tgd)
            ll_info.clear(); l_align.clear(); d_split.clear()
        if not op or op=="merge":
            timeToMerge(tgd,spk,l_tiers[1])          # Fuse them
        end = time.time()-start
        print("\t",tgd.name+", "+spk+"  |  ","{}: {:.4f}s".format(spk,end))
    return tgd
    
def main(idir=idir,odir=odir,d_vals={},log_path=""):
    """Call for master function."""
    I = reinject.Inj(); I.d['d_tiers'] = reinject.fillSPK(I.d['tiers_csv'])
    global glog_path
    glog_path = log_path
    if not glog_path:
        glog_path = os.path.abspath("log_splitMerge.txt")
    startLog(glog_path)
    for fi,ext,file,path in iter(idir):
        reinject.resetInj(I)
        log(fi+"\n",glog_path); start = time.time()
        #try:
            # Deal with 'reinject' stuff
        I.d['eaf'] = reinject.getEAF(idir,file)         # ELAN file
        I.d['tgd'] = reinject.getTGD(idir,fi,
        npath = os.path.join(odir,fi+".eaf")
        shutil.copy(path,npath)                         # copy EAF (raw)
                            l_path=glog_path)           # TextGrid file
        if not I.d['tgd']:
            log("\tNo TextGrid found.\n",glog_path); continue
        getIStuff(I)                                    # core,lang,etc.
        wdAlign(I)                                      # word alignment
        end = time.time()-start
        log("\twdAlign: {:.4f}s\n".format(end),glog_path)
        print(fi,I.d['lang'])
        d_g2p = getG2P(I.d['g2p'],I.d['lang'])
            # Back to the old script
        splitnMerge(I.d['eaf'],I.d['tgd'],I.d['d_core'],d_g2p)
        npath = os.path.join(odir,fi+"_Algn.TextGrid") # Export TextGrid
        toPraat.toPraat(npath,I.d['tgd'])
        #except:
        #    log("\t-- FAILURE --\n")

if __name__ == '__main__':
    main()