from corflow import fromPangloss,fromElan,toElan
import os,re

home = os.path.abspath(os.path.dirname(__file__))
indir = os.path.join(home,"tmp")
outdir = os.path.join(home,"tmp")

for file in os.listdir(indir):
    fi,ext = os.path.splitext(file)
    path = os.path.join(indir,file)
    if not ext.lower() == ".xml":
        continue
    print(fi)
    trans = fromPangloss.fromPangloss(path)
    npath = os.path.join(outdir,fi+".eaf")
    toElan.toElan(npath,trans)
    trans = fromElan.fromElan(npath)
    d_lvls = {'S':False,'W':False,'M':False}
    d_top = {'S':None,'W':'S','M':'W'}
    for tier in trans:
        tier.setParent(None)
        for seg in tier:
            seg.setParent(None)
        if "-" in tier.name:
            lvl = tier.name.split("-",1)[0]
            if lvl in d_lvls and not d_lvls[lvl]:
                d_lvls[lvl] = tier
    for tier in trans:
        if "-" in tier.name:
            lvl = tier.name.split("-",1)[0]
            ptier = d_lvls.get(lvl); nlvl = d_top.get(lvl)
            if not ptier:
                continue
            elif ptier == tier:
                if not nlvl:
                    continue
                ptier = d_lvls.get(nlvl)
                if not ptier:
                    continue
            tier.setParent(ptier)
            for a in range(len(tier)-1,-1,-1):
                seg = tier.elem[a]
                mid = seg.start+((seg.end-seg.start)/2)
                pseg = ptier.getTime(mid,ptier)
                if not pseg:
                    tier.pop(a)
                else:
                    seg.setParent(pseg)
    toElan.toElan(npath,trans)