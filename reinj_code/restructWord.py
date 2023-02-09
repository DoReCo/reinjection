 # -*- coding: utf-8 -*-
"""restructWord
Replaces the original word tier ('wd') with the MAUS word tier ('W_O_tx').

A pair of EAF+Tgd files must have been loaded in an 'Inj' object, then
run through 'wordAlign()' and 'getMatches()' in 'wordAlign.py' first.

Just call 'restructWord()' with that 'Inj' object as argument. It will:
    0. Merges some MAUS segments.
       > fixMaus()
    1. Replace the original word tier ('w_tier') with the MAUs one ('m_tier').
       > assignRefs()
    2. Realign all segments reassigned to the MAUS words, delete others.
       > alignSegs()
    3. Add dummy/filler segments for (non-)pauses on all tiers.
       > addFillSegs()
Note that 'alignSegs()' and 'addFillSegs()' require structural information.
This is provided by a method in 'Inj', 'getStructure()'.

You can then either call 'Inj.export()' or use the 'Transcription' object
that 'restructWord()' returns to export a restructured EAF file.
"""

    #### REASSIGN SEGMENTS ####
def assignRefs(trans,w_tier,m_tier,l_porig,l_pmaus):
    """Reassigns ids/refs and pointers.
    ARGUMENTS:
    - 'w_tier'  :       (pntr) original word tier
    - 'm_tier'  :       (pntr) MAUS word tier
    - 'l_porig' :       (lst<tpl<str,int,pntr>>) orig' word tier pointers
    - 'l_pmaus' :       (lst<tpl<str,int,pntr>>) MAUS word tier pointers
    RETURNS:
    - 'w_tier/m_tier' with 'm_tier' having replaced 'w_tier'
      in the structure.
    Note: 'l_porig/l_pmaus' mimic the pairwise output (see 'wordAlign.py',
          'getPntr()'). The two lists must be of same length and contain
          tuples (id,index,segment) as values.
    Note: If no segment is assigned, the tuple will be ("",-1,None)."""

    def iterRefs(l_porig,l_pmaus,cho=True):
        """To iterate over the two lists of references."""
        for a in range(len(l_porig)):
            mid,mind,mseg = l_pmaus[a]; oid,oind,oseg = l_porig[a]
            if cho and not oseg:
                continue
            yield (mid,mind,mseg,oid,oind,oseg)
    def assignTiers(w_tier,m_tier):
        """Putting together everything about switching each tier."""
            # Avoid pointer error
        m_ind = trans.d_elem[m_tier][0]; w_ind = w_tier.index()
        m_tier = trans.elem[m_ind]; w_tier = trans.elem[w_ind]
            # Sort
        m_tier.move(w_ind); w_tier.move(-1)
        
        m_ind = m_tier.index(); w_ind = w_tier.index()
            # Name
        wd_spk = ""
        if not 'speaker' in w_tier.metadata:
            wd_spk = w_tier.name.rsplit("@",1)[1]
        else:
            wd_spk = w_tier.metadata['speaker']
        d_tiers = {tier.name:tier for tier in trans}; incr = 1
        m_tier.name = w_tier.name; w_tier.name = "wd_old@"+wd_spk
        while w_tier.name in d_tiers:
            w_tier.name = "wd_old"+str(incr)+"@"+wd_spk; incr+= 1
            # Type
        w_tier.setMeta('type',"orig"); m_tier.setMeta('type',"wd")
            # Assign new parent
        m_tier.setParent(w_tier.parent())
        w_tier.setParent(None)
        for child in w_tier.children():
            m_tier.addChild(child)
    
        # Deal with the tier level
    assignTiers(w_tier,m_tier)
        # Deal with the segment level
    d_csegs = {}
    for ctier in m_tier.children():
        for cseg in ctier:
            wseg = cseg.parent()
            if wseg in d_csegs:
                d_csegs[wseg].append(cseg)
            else:
                d_csegs[wseg] = [cseg]
    for mid,mind,mseg,oid,oind,oseg in iterRefs(l_porig,l_pmaus,False):
        if not mseg or not oseg:    # Ignore non-matches
            continue
        pseg = oseg.parent()
        mseg.metadata = oseg.metadata.copy()
        if pseg:                    # Parent
            mseg.setParent(pseg)
        if oseg in d_csegs:
            for cseg in d_csegs[oseg]:
                cseg.setParent(mseg)
        # Clear the rest
    for seg in w_tier:
        seg.setParent(None)
    for oseg,l_csegs in d_csegs.items():
        for cseg in l_csegs:
            if cseg.parent() == oseg:
                cseg.setParent(None)
    return (w_tier,m_tier)
    #### TIME CODES / REMOVAL ####
def alignSegs(I,spk,trans,w_tier,m_tier):
    """Removes segments without a parent.
    Also realigns parent segments.
    Exceptions: 'w_tier', 'm_tier' and independent tiers."""

        ## Parents
    ptier = tier = m_tier
    while ptier.parent(): # We loop from m_tier to the top parent
        tier = ptier; ptier = ptier.parent()
            # Get dict' <parent_seg: seg>
        d_segs = {}
        for seg in tier:
            pseg = seg.parent()
            if pseg:                            # parent seg
                if pseg in d_segs:
                    d_segs[pseg].append(seg)
                    pseg.end = seg.end          # sneak align
                else:
                    d_segs[pseg] = [seg]
                    pseg.start = seg.start      # sneak align
                    pseg.end = seg.end
            # Remove all segs not in that dict'
        for a in range(len(ptier)-1,-1,-1):
            pseg = ptier.elem[a]
            if pseg not in d_segs:
                ptier.pop(a)
        ## Children
    l_par = [ptier]
    while l_par:
        d_psegs = {}; l_child = []
        for otier in l_par:     # dict' <pseg:True> for all parent tiers
            for pseg in otier:
                d_psegs[pseg] = True
            l_child = l_child+otier.children()
        for child in l_child:   # children tiers
            if child == m_tier: # ignore m_tier plz
                continue
                # Remove
            for a in range(len(child)-1,-1,-1):
                cseg = child.elem[a]; pseg = cseg.parent()
                if (not pseg) or (not pseg in d_psegs):
                    child.pop(a)
                # Align
                ## Note: It should be impossible for 'oseg' to not exist
            oseg = None; l_segs = []
            for cseg in child:
                pseg = cseg.parent()
                if not pseg == oseg:
                    child._split(l_segs,oseg)
                    oseg = pseg; l_segs = [cseg]
                else:
                    l_segs.append(cseg)
            if oseg:
                child._split(l_segs,oseg)
        l_par = l_child

    return trans
    #### DUMMY/FILLER SEGMENTS ####
def addFillSegs(I,spk,m_tier,id,c):
    """Adds filler segments EVERYWHERE.
    'sym' for non-pause segments. 'psym' for pause segments."""

    def addSeg(tier,ind,start,end,sym,id,c,seg=None):
        nseg = tier.create(ind,id+str(c),start,end,sym,struct=tier); c += 1
        if seg and seg.parent():
            nseg.setParent(seg.parent())
        if ind-1 >= 0:
            prev = tier.elem[ind-1]
            if prev.end > start:
                prev.end = start
        if (ind+1 < len(tier)) and (not ind+1 == 0):
            next = tier.elem[ind+1]
            if next.start < end:
                next.start = end
        return (c,ind,nseg)
    
    def findContext(tier,ind,pseg):
        """Returns a context of segments.
        Yes, I know, future me, very descriptive."""
        
        prev = (-1,None); l_segs = []; next = (-1,None)
        for a in range(ind,-1,-1):
            seg = tier.elem[a]
            if (seg.end-pseg.start) < 0.0001:
                prev = (a,seg); ind = a; break
            elif (seg.start-pseg.end) > -0.0001:
                next = (a,seg)
            else:
                l_segs.append((a,seg))
        return (ind,prev,l_segs,next)
    def fillPauses(ntrans,ntier,nsym,c):
        """Adding pause fillers for 'm_tier'."""

            # End
        if not ntier.elem[-1].content:
            ntier.elem[-1].content = nsym
            ntier.elem[-1].end = trans.end
        elif ntier.elem[-1].end < ntrans.end:
            c,ind,nseg = addSeg(ntier,len(ntier),ntier.elem[-1].end,
                                ntrans.end,nsym,id,c)
            # Middle
        for a in range(len(ntier)-1,0,-1):
                # By time
            seg1 = ntier.elem[a-1]; seg2 = ntier.elem[a]
            if (seg1.end < seg2.start):
                c,ind,nseg = addSeg(ntier,a,seg1.end,seg2.start,nsym,id,c)
                # ... By content
            elif not seg2.content:
                seg2.content = nsym
            # Start
        if not ntier.elem[0].content:
            ntier.elem[0].content = nsym
            ntier.elem[0].start = trans.start
        elif ntier.elem[0].start > ntrans.start:
            c,ind,nseg = addSeg(ntier,0,ntrans.start,
                                ntier.elem[0].start,nsym,id,c)
            # Just cleaning up a bit
        if ntier.elem:
            ntier.start = ntier.elem[0].start
            ntier.end = ntier.elem[-1].end
        else:
            ntier.start = ntrans.start; ntier.end = ntrans.end
        return c
    def fillOtherPauses(tier,ptier,ind,nsym,id,c):
        """Fills the pauses of 'ptier' using 'tier' pauses."""
            # For each pause of 'tier'
        for a in range(len(tier)-1,-1,-1):
            seg = tier.elem[a]
            if not seg.content == nsym:
                continue
                # We get the previous segment from 'ptier'
            ind,prev,l_segs,next = findContext(ptier,ind,seg)
            i,prev = prev
                # Little tweak to deal with index (ind) == 0
            i = i+1 if prev else 0
                #  We do nothing if there is already something there
            if i < len(ptier) and (ptier.elem[i].start <= seg.start and
               ptier.elem[i].end >= seg.end):
                continue
                # We add a pause at that index
            c,i,nseg = addSeg(ptier,i,seg.start,seg.end,seg.content,id,c)
        if ptier.elem:
            ptier.start = ptier.elem[0].start
            ptier.end = ptier.elem[-1].end
        else:
            ptier.start = trans.start; ptier.end = trans.end
        return ind,c
    def fillMore(tier,id,c):
        """Add non-pauses filler segments."""
        
        n="DoReCo_"; d=1; cont = "{}{:04d}".format(n,d)
            # End case
        if tier.elem and tier.elem[-1].end < trans.end:
            c,i,nseg = addSeg(tier,-1,tier.elem[-1].end,
                              trans.end,cont,id,c)
            d += 1; cont = "{}{:04d}".format(n,d)
            # Middle case
        if len(tier) > 1:
            for a in range(len(tier)-1,0,-1):
                seg1 = tier.elem[a-1]; seg2 = tier.elem[a]
                if seg1.end < seg2.start:
                    c,i,nseg = addSeg(tier,a,seg1.end,seg2.start,cont,id,c)
                    d += 1; cont = "{}{:04d}".format(n,d)
            # Start case
        if tier.elem and (tier.elem[0].start > trans.start):
            c,i,nseg = addSeg(tier,0,trans.start,tier.elem[0].start,cont,
                              id,c)
            d += 1; cont = "{}{:04d}".format(n,d)
        return c
    def fillOtherMore(tier,ptier,sym,id,c):
        """Add non-pause fillers but takes 'ptier' into account."""
        
        ind = len(tier)-1
            # For each parent segment
        for a in range(len(ptier)-1,-1,-1):
            pseg = ptier.elem[a]
                # Get the 'tier' context
            ind,prev,l_segs,next = findContext(tier,ind,pseg)
                # None
            if not l_segs:
                    # We need an index for the new segment
                i = -1
                if prev[1]:
                    i = prev[0]+1
                else:
                    i = 0
                    # Add using that index
                c,i,nseg = addSeg(tier,i,pseg.start,pseg.end,sym,id,c)
                nseg.setParent(pseg)
                continue
                # End case
            if l_segs[0][1].end < pseg.end:
                c,i,nseg = addSeg(tier,l_segs[0][0]+1,l_segs[0][1].end,
                                  pseg.end,sym,id,c)
                nseg.setParent(pseg)
                # Middle case
            if len(l_segs) > 1:
                for b in range(1,len(l_segs)):
                    i1,seg1 = l_segs[b-1]; i2,seg2 = l_segs[b]
                    if seg1.end < seg2.start:
                        c,i,nseg = addSeg(tier,i2,seg1.end,seg2.start,sym,
                                          id,c)
                        nseg.setParent(pseg)
                # Start case
            if l_segs[-1][1].start > pseg.start:
                c,i,nseg = addSeg(tier,l_segs[-1][0],pseg.start,
                                  l_segs[-1][1].start,sym,id,c)
                nseg.setParent(pseg)
        return c
        # Variables
    trans,sym,psym = I.d['eaf'],I.d['sym'],I.d['psym']
        # Fill 'm_tier' with pauses
    c = fillPauses(trans,m_tier,psym,c)
    for tier in trans:          # because of morphAlign(), fill 'ms/ph'
        typ = tier.meta('type')
        if typ == "ms" or typ == "ph":
            c = fillPauses(trans,tier,psym,c)

        # Fill top parent
        # Get said top parent (ptier)
    ptier = m_tier
    while ptier.parent():
        ptier = ptier.parent()
    ind = len(ptier)-1
        # Add parent pauses first
    ind,c = fillOtherPauses(m_tier,ptier,ind,psym,id,c)
        # Fill parent
    c = fillMore(ptier,id,c)

        # Fill the children
        # (And renew the segment references)
    l_child = ptier.children()
    while l_child:
        l_tmp = []
        for ctier in l_child:
            ptier = ctier.parent(); l_tmp = l_tmp+ctier.children()
            if ctier not in I.d['d_core'][spk]:
                continue
            ind = len(ctier)-1
                # Add child pauses first
            ind,c = fillOtherPauses(ptier,ctier,ind,psym,id,c)
                # Fill child
            c = fillOtherMore(ctier,ptier,sym,id,c)
                # Renew references
            for cseg in ctier:
                pseg = cseg.getTime(cseg.start,ptier)
                cseg.setParent(pseg)
        l_child = l_tmp; l_tmp = []
    return trans,c
    #### Split phone tier ####
def addPhones(I,spk,d_styp,m_tier):
    """In preparation for wordAlign, we split (<notProcessedChunk>)."""

    mb_tier,ph_tier = d_styp['mb'],d_styp['ph'] # Moar tiers
    if not mb_tier or not ph_tier:
        return
    for wseg in m_tier:                         # For each word segment
        l_mbsegs = wseg.childDict()[mb_tier]    # morphemes
        ph_n,ph_ind,phseg = wseg.getTime(wseg.start,ph_tier,det=True)
        l_phsegs = []                           # get phonemes
        for b in range(ph_ind,len(ph_tier)):
            l_phsegs.append(ph_tier.elem[b])
            if ph_tier.elem[b].end >= wseg.end:
                break
        if not(len(l_phsegs) == 1 and l_phsegs[0].content.startswith("<not")):
            continue
        lmb = len(l_mbsegs)                     # Compare lengths
        if len(l_phsegs) < lmb:                 # And create missing segs
            i = len(l_mbsegs)-len(l_phsegs); dur = (wseg.end-wseg.start)/lmb
            phseg.start = wseg.start+(dur*i)
            for b in range(i-1,-1,-1):
                s = wseg.start+(dur*b); e = wseg.start+(dur*(b+1))
                nseg = ph_tier.create(ph_ind,start=s,end=e,
                                      content=phseg.content)
                l_phsegs.append(phseg)
    #### Rename reference segments ####
def refTier(I,spk,d_styp):
    """Replaces the content of the reference tier segments."""
    
    rtier = d_styp['ref']; incr = 0
    if not rtier:   # absurd
        return
    for seg in rtier:
        if "<p:>" in seg.content:
            continue
        txt = seg.content if not seg.content.startswith("DoReCo") else ""
        seg.content = "{:04d}_DoReCo_{}".format(incr,I.d['eaf'].name)
        incr += 1
        if txt:
            seg.content = seg.content+"_"+txt

    #### MAIN ####
def restructWord(I):
    """A function to restructure the pairwise output."""

    id="resWd";c=0
    for spk,out in I.d['output'].items():
        if not out:
            refTier(I,spk,I.d['d_typ'][spk]); continue
        l_porig,l_pmaus = out; d_styp = I.d['d_typ'][spk]
        trans,w_tier,m_tier = I.d['eaf'],d_styp['wd'],d_styp['maus']
            # Reassign parent/child tiers/segments
        assignRefs(trans,w_tier,m_tier,l_porig,l_pmaus)
        I.d['output'][spk] = None; del l_porig; del l_pmaus
            # Realign / remove parents/children segments
        I.d['eaf'] = alignSegs(I,spk,trans,w_tier,m_tier)
            # Filling segments
        I.d['eaf'],c = addFillSegs(I,spk,m_tier,id,c)
            # Splitting the phoneme tier
        addPhones(I,spk,d_styp,m_tier)
            # Rename reference segments
        refTier(I,spk,d_styp)
    
    return I.d['eaf']
    











