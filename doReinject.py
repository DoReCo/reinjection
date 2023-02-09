import splitMergeTGD,reinject,doExport
import sys,os,re,time

    # Paths
home = os.path.abspath(os.path.dirname(__file__))
idir = os.path.join(home,"01_splitmerged")
odir = os.path.join(home,"02_reinjected")

def iterIn(root,files):
    for file in files:
        fi,ext = os.path.splitext(file)
        path = os.path.join(root,file)
        yield fi,ext,file,path
def suff(p):
    fi,ext = os.path.splitext(os.path.basename(p))
    return int(fi.rsplit("_",1)[1])

def main(idir,odir):
    for root,dirs,files in os.walk(idir):
        if root == idir:
            continue
        d = os.path.basename(root)
        print(d)
        nd = os.path.join(odir,d)
        if not os.path.isdir(nd):
            os.mkdir(nd)
        #splitMergeTGD.main(root,nd)
        reinject.reinject(d_vals = {'in_dir':root,'out_dir':nd})
        """continue
        l_files = []
        for fi,ext,file,path in iterIn(root,files):
            l_files.append(path)
        l_files.sort(key=suff)
        if not nd:
            os.mkdir(nd)
        print(l_files)
        doExport.merge(l_files[0],l_files[1],nd)"""
        
if __name__ == "__main__":
    main(idir,odir)