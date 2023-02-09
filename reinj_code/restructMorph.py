
    #### RESTRUCTURE ####
def assignMorphRefs(I,out,ph_tier,mb_tier,ms_tier):
    """Assigns ids/refs to the phoneme tier."""
    
        # Variables
    l_porig,l_pmaus,l_cert = out
        # Reassignment
    for a in range(len(l_porig)):               # phone-to-morph
        l_pmaus[a][2].setParent(l_porig[a][2])
        # Attach the phone tier to the morpheme tier
    ph_tier.setParent(mb_tier)
        # Pop the weasel
    if ms_tier:
        I.d['eaf'].pop(ms_tier.index())
    return I.d['eaf']
def alignMorphSegs(I,spk,ph_tier,mb_tier):
    """Realigns morphemes."""
    
    lmb = len(mb_tier)
    for a in range(lmb-1,-1,-1):
        mbseg = mb_tier.elem[a]; start = 0.; end = 0.
        d_children = mbseg.childDict()
        if not ph_tier in d_children:           # 0-duration morphemes
            #for tier,l_csegs in d_children.items():# remove
            #    for cseg in l_csegs:
            #        tier.remove(cseg)
            #mb_tier.pop(a)
            if a+1 < lmb:                    # 0-duration
               start = mb_tier.elem[a+1].start; end = start
        else:                                   # time codes
            l_phsegs = d_children[ph_tier]
            if l_phsegs:
                start = l_phsegs[0].start; end = l_phsegs[-1].end
        mbseg.start = start; mbseg.end = end    # Realign morpheme
    l_stop = [ph_tier]                          # Deal with mb_tier children
    for tier in I.d['d_core'][spk]:
        l_stop.append(tier)
        # setChildTime
    l_par = []
    for child in mb_tier.allChildren():
        par = child.parent()
        if not par in l_par:
            par.setChildTime(ch=False,stop=[ph_tier])
            l_par.append(par)
        
    return I.d['eaf']
def createCertTier(I,spk,out,ph_tier,mb_tier,c):
    """Creates the 'uncertainty' tier.
    Note: doesn't fill dummy/filler segments (pauses, non-pauses)."""
    
        # Variables
    trans = I.d['eaf']
    l_cert = out[2]; id = "aCert"
        # Uncertainty tier
    utier = trans.create(ph_tier.index()+1,"doreco-mb-algn@"+spk,
                         trans.start,trans.end)
    utier.setParent(mb_tier)
    utier.setMeta('speaker',spk); utier.setMeta('type',"cert")
    for mbseg,cont in l_cert:
        if cont and isinstance(cont,str):
            nseg = utier.create(-1,id+str(c),mbseg.start,mbseg.end,cont)
            nseg.setParent(mbseg); c += 1
    return trans,c
def fixWIP(mb_tier,ph_tier):
    """Splits everything at '<<wip>>'."""

    for a in range(len(mb_tier)-1,-1,-1):
        mb_seg = mb_tier.elem[a]
        d_child = mb_seg._childDict(mb_seg.children())
            # Find a '<<wip>>'
        l_phsegs = d_child.get(ph_tier)
        ch_wip = None
        if not l_phsegs:
            continue
        for phseg in l_phsegs:
            if phseg.content == "<<wip>>":
                ch_wip = phseg; break
        if not ch_wip:
            continue
            # Split
        if mb_tier.name.startswith("mb@"):          # normal
            wn,wi,w_seg = mb_seg.parent(det=True)
        else:                                       # _alternate()
            wn,wi,w_seg = mb_seg.name,a,mb_seg
        ptier = w_seg.parent()
        ptier = ptier.struct if ptier else None
        if "<<wip>>" in w_seg.content:          # word
            prev,nex = w_seg.content.rsplit("<<wip>>",1)
            nseg = w_seg.struct.create(wi,"",ch_wip.start,ch_wip.end,"<<wip>>")
            if ptier:
                nseg.setParent(ptier.getTime(nseg.start,ptier))
            if not prev:
                w_seg.start = ch_wip.end; w_seg.content = nex
            elif not nex:
                w_seg.end = ch_wip.start; w_seg.content = prev
            else:
                nseg = w_seg.struct.create(wi,"",w_seg.start,ch_wip.start,prev)
                if ptier:
                    nseg.setParent(ptier.getTime(nseg.start,ptier))
                w_seg.start = ch_wip.end; w_seg.content = nex
        d_child = w_seg._childDict(w_seg.allChildren())
        for ctier,l_csegs in d_child.items():   # children
            if ctier == ph_tier:
                continue
            for cseg in l_csegs:
                if cseg.end <= ch_wip.start:
                    continue
                elif cseg.start == ch_wip.start and cseg.end == ch_wip.end:
                    break
                ptier = cseg.parent()
                ptier = ptier.struct if ptier else None
                    ## Only two cases considered
                nseg = None
                if cseg.start == ch_wip.start:  # WIP at start
                    nseg = ctier.create(cseg.index(),"",cseg.start,
                                        ch_wip.end,"<<wip>>")
                    if ptier:
                        nseg.setParent(ptier.getTime(nseg.start,ptier))
                    cseg.start = ch_wip.end
                elif cseg.end == ch_wip.end:    # WIP at end
                    nseg = ctier.create(cseg.index()+1,"",ch_wip.start,
                                        cseg.end,"<<wip>>")
                    if ptier:
                        nseg.setParent(ptier.getTime(nseg.start,ptier))
                    cseg.end = ch_wip.start
                break
        # Reparent
    """for ctier in mb_tier.allChildren():
        if ctier == ph_tier:
            continue
        ptier = ctier.parent()
        for cseg in ctier:
            cseg.setParent(ptier.getTime(cseg.start,ptier))"""

    # Main
def _alternate(I,d_styp):
    """Without a morpheme tier, simple cleanup."""
    ph_tier,ms_tier,wd_tier = d_styp['ph'], d_styp['ms'],d_styp['wd']
    if not wd_tier:
        print("\tNo word ('wd') tier for aligned speaker, "
              "check file/metadata.")
        return
    if ms_tier:
        I.d['eaf'].pop(ms_tier.index())
    ph_tier.setParent(wd_tier)
    for seg in ph_tier:
        pseg = ph_tier.getTime(seg.start,wd_tier)
        seg.setParent(pseg)
    fixWIP(wd_tier,ph_tier)
def restructMorph(I):
    """Uses lists of references to reassign/realign morphs/phones."""

    for spk in I.d['d_core']:
            # Variables
        d_styp = I.d['d_typ'][spk]
        ph_tier,ms_tier,mb_tier = d_styp['ph'], d_styp['ms'],d_styp['mb']
            # No tiers? Move on
        if not ph_tier or not mb_tier:
            if ph_tier:
                _alternate(I,d_styp)
            continue
        out = I.d['output'][spk]; c = 0
            # assigns phonemes to morphemes
        assignMorphRefs(I,out,ph_tier,mb_tier,ms_tier)
            # realigns morphemes
        alignMorphSegs(I,spk,ph_tier,mb_tier)
            # Deal with WIP
        fixWIP(mb_tier,ph_tier)
            # Creates the 'uncertainty' tier
        _,c = createCertTier(I,spk,out,ph_tier,mb_tier,c)
    return I.d['eaf']