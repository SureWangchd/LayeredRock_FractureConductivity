"""Conforming local tetrahedral fit of an analytic oblate ellipsoid.
Geometry preprocessing only, NOT an electrical solver.
Every output tetrahedron is a child of one original CT tetrahedron.
Intersections cut tetrahedra; no centroid-based fracture tagging is used.
"""
from __future__ import annotations
import itertools, json, math, time
from pathlib import Path
import numpy as np
from numba import njit

EDGES=np.asarray([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],dtype=np.int64)
FACES=np.asarray([(1,2,3),(0,3,2),(0,1,3),(0,2,1)],dtype=np.int64)
U=np.asarray([.5,0.,math.sqrt(3)/2]); VV=np.asarray([0.,1.,0.]); NN=np.asarray([-math.sqrt(3)/2,0.,.5])

@njit(cache=True)
def sd_one(x,y,z,a=.25,c=.0125):
    """Exact signed Euclidean distance to a=b oblate spheroid, normalized L=1."""
    x-=.5;y-=.5;z-=.5
    u=.5*x+math.sqrt(3)/2*z
    w=abs(-math.sqrt(3)/2*x+.5*z)
    r=math.sqrt(u*u+y*y); a2=a*a;c2=c*c
    q=r*r/a2+w*w/c2
    if abs(q-1.)<=2e-14:return 0.
    if r<1e-16:return w-c
    if w<1e-16:
        transition=(a2-c2)/a
        if r>=transition:return r-a
        xx=a2*r/(a2-c2)
        zz=c*math.sqrt(max(0.,1.-xx*xx/a2))
        return -math.sqrt((xx-r)**2+zz*zz)
    lo=0. if q>1. else -c2
    hi=math.sqrt(a2*r*r+c2*w*w) if q>1. else 0.
    for _ in range(62):
        t=.5*(lo+hi)
        f=a2*r*r/((t+a2)**2)+c2*w*w/((t+c2)**2)-1.
        if f>0.:lo=t
        else:hi=t
    t=.5*(lo+hi)
    xx=a2*r/(t+a2);zz=c2*w/(t+c2)
    d=math.sqrt((xx-r)**2+(zz-w)**2)
    return d if q>1. else -d

@njit(cache=True)
def distances(v, snap=0.):
    d=np.empty(len(v))
    for i in range(len(v)):
        d[i]=sd_one(v[i,0],v[i,1],v[i,2])
        if abs(d[i])<max(2e-14,snap):d[i]=0.
    return d

@njit(cache=True)
def tetvol(v,t):
    out=np.empty(len(t))
    for i in range(len(t)):
        a=v[t[i,1]]-v[t[i,0]];b=v[t[i,2]]-v[t[i,0]];c=v[t[i,3]]-v[t[i,0]]
        out[i]=(a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))/6
    return out

@njit(cache=True)
def tri_origin_dist2(a,b,c):
    # Ericson's closest point on a triangle, including edge and vertex Voronoi regions.
    ab=b-a;ac=c-a;ap=-a;d1=np.dot(ab,ap);d2=np.dot(ac,ap)
    if d1<=0. and d2<=0.:return np.dot(a,a)
    bp=-b;d3=np.dot(ab,bp);d4=np.dot(ac,bp)
    if d3>=0. and d4<=d3:return np.dot(b,b)
    vc=d1*d4-d3*d2
    if vc<=0. and d1>=0. and d3<=0.:
        p=a+d1/(d1-d3)*ab;return np.dot(p,p)
    cp=-c;d5=np.dot(ab,cp);d6=np.dot(ac,cp)
    if d6>=0. and d5<=d6:return np.dot(c,c)
    vb=d5*d2-d1*d6
    if vb<=0. and d2>=0. and d6<=0.:
        p=a+d2/(d2-d6)*ac;return np.dot(p,p)
    va=d3*d6-d5*d4
    if va<=0. and (d4-d3)>=0. and (d5-d6)>=0.:
        p=b+(d4-d3)/((d4-d3)+(d5-d6))*(c-b);return np.dot(p,p)
    den=va+vb+vc
    if den==0.:return min(np.dot(a,a),np.dot(b,b),np.dot(c,c))
    p=a+vb/den*ab+vc/den*ac;return np.dot(p,p)

@njit(cache=True)
def contains_hidden_surface(p, tolerance):
    # AABB is only an inexpensive exclusion test; then test actual tetrahedron.
    r=np.empty((4,3))
    for j in range(4):
        x=p[j,0]-.5;y=p[j,1]-.5;z=p[j,2]-.5
        r[j,0]=(.5*x+math.sqrt(3)/2*z)/.25;r[j,1]=y/.25;r[j,2]=(-math.sqrt(3)/2*x+.5*z)/.0125
    lower=0.
    for k in range(3):
        lo=np.min(r[:,k]);hi=np.max(r[:,k])
        if lo>0.:lower+=lo*lo
        elif hi<0.:lower+=hi*hi
    if lower>=1.-1e-12:return False
    # Is origin inside tetrahedron? signed volumes of four sub-tets.
    m=np.empty((3,3))
    for j in range(3):m[:,j]=r[j+1]-r[0]
    b=np.linalg.solve(m,-r[0])
    if np.min(b)>=0. and np.sum(b)<=1.:return True
    mn=1e100
    for f in range(4):
        q=tri_origin_dist2(r[FACES[f,0]],r[FACES[f,1]],r[FACES[f,2]])
        mn=min(mn,q)
    return mn < (1.-tolerance/.25)**2

@njit(cache=True)
def mark_geometry(v,t,d,tolerance):
    mark=np.zeros(len(t),np.bool_);maxerr=0.;hidden=0;cross=0
    for i in range(len(t)):
        ids=t[i];s=d[ids];mn=np.min(s);mx=np.max(s)
        if mx<=0.:continue
        p=v[ids]
        if mn>=0.:
            # Ellipsoid is convex; tetrahedra wholly inside cannot hide an outside pocket.
            # Same-sign OUTSIDE tetrahedra may hide a thin rim: test it, do not just discard.
            if contains_hidden_surface(p,tolerance):mark[i]=True;hidden+=1
            continue
        cross+=1
        points=np.empty((10,3));n=0
        for j in range(4):
            if s[j]==0.:
                points[n]=p[j];n+=1
        for e in range(6):
            j,k=EDGES[e]
            if s[j]*s[k]<0.:
                w=s[j]/(s[j]-s[k]);points[n]=p[j]+w*(p[k]-p[j]);n+=1
        err=0.
        for j in range(n):
            err=max(err,abs(sd_one(points[j,0],points[j,1],points[j,2])))
            for k in range(j):
                z=.5*(points[j]+points[k]);err=max(err,abs(sd_one(z[0],z[1],z[2])))
        z=np.zeros(3)
        for j in range(n):z+=points[j]/n
        err=max(err,abs(sd_one(z[0],z[1],z[2])))
        maxerr=max(maxerr,err)
        if err>tolerance:mark[i]=True
    return mark,maxerr,hidden,cross

@njit(cache=True)
def edge_keys(t,nv):
    out=np.empty(len(t)*6,np.int64)
    for i in range(len(t)):
        for e in range(6):
            a=t[i,EDGES[e,0]];b=t[i,EDGES[e,1]]
            if a>b:a,b=b,a
            out[6*i+e]=a*nv+b
    return out


def make_clip_tables():
    """Topology only. Global-ID pulling triangulation is performed at runtime."""
    poly=np.full((81,2,5,4),-1,np.int64);lens=np.zeros((81,2,5),np.int64)
    p=np.asarray([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)
    edgeid={tuple(e):4+i for i,e in enumerate(EDGES)}
    def intersection(a,b,s):
        if s[a]==0:return a
        if s[b]==0:return b
        return edgeid[tuple(sorted((a,b)))]
    for state in range(81):
        s=np.asarray([(state//3**i)%3-1 for i in range(4)])
        if not (np.min(s)<0 and np.max(s)>0):continue
        coords=np.concatenate([p,np.asarray([.5*(p[a]+p[b]) for a,b in EDGES])])
        plane=[i for i in range(4) if s[i]==0]+[4+i for i,(a,b) in enumerate(EDGES) if s[a]*s[b]<0]
        pp=coords[plane];normal=np.cross(pp[1]-pp[0],pp[2]-pp[0]);normal/=np.linalg.norm(normal)
        ex=pp[0]-pp.mean(0);ex/=np.linalg.norm(ex);ey=np.cross(normal,ex)
        ang=np.arctan2((pp-pp.mean(0))@ey,(pp-pp.mean(0))@ex);plane=[plane[i] for i in np.argsort(ang)]
        for side in range(2):
            faces=[]
            for face in FACES:
                z=[]
                for i in range(3):
                    a=int(face[i]);b=int(face[(i+1)%3]);aa=s[a]<=0 if side==0 else s[a]>=0;bb=s[b]<=0 if side==0 else s[b]>=0
                    if aa:z.append(a)
                    if aa!=bb:z.append(intersection(a,b,s))
                z=list(dict.fromkeys(z))
                if len(z)>=3:faces.append(z)
            faces.append(plane)
            for f,z in enumerate(faces):poly[state,side,f,:len(z)]=z;lens[state,side,f]=len(z)
    return poly,lens
POLY,LENS=make_clip_tables()

@njit(cache=True)
def clipping_edges(t,d,nv):
    out=np.empty(len(t)*6,np.int64);n=0
    for i in range(len(t)):
        for e in range(6):
            a=t[i,EDGES[e,0]];b=t[i,EDGES[e,1]]
            if d[a]*d[b]<0.:
                if a>b:a,b=b,a
                out[n]=a*nv+b;n+=1
    return out[:n]

@njit(cache=True)
def clip_child_list(tet,d,nv,selected,poly,lens):
    loc=np.empty(10,np.int64);loc[:4]=tet;loc[4:]=-1
    state=0;mn=1e300;mx=-1e300
    for i in range(4):
        s=d[tet[i]];mn=min(mn,s);mx=max(mx,s)
        state+=(0 if s<0. else (2 if s>0. else 1))*3**i
    out=np.empty((8,4),np.int64);inside=np.zeros(8,np.uint8)
    if mn>=0. or mx<=0.:
        out[0]=tet;inside[0]=1 if mx<=0. else 0;return out,inside,1
    for e in range(6):
        a=tet[EDGES[e,0]];b=tet[EDGES[e,1]]
        if d[a]*d[b]<0.:
            if a>b:a,b=b,a
            key=a*nv+b;idx=np.searchsorted(selected,key)
            if idx>=len(selected) or selected[idx]!=key:raise ValueError('shared intersection edge missing')
            loc[4+e]=nv+idx
    nt=0
    for side in range(2):
        apex=9223372036854775807
        for f in range(5):
            for j in range(lens[state,side,f]):apex=min(apex,loc[poly[state,side,f,j]])
        for f in range(5):
            n=lens[state,side,f]
            if n<3:continue
            face=np.empty(4,np.int64);hasapex=False
            for j in range(n):
                face[j]=loc[poly[state,side,f,j]]
                if face[j]==apex:hasapex=True
            if hasapex:continue
            p=0
            for j in range(1,n):
                if face[j]<face[p]:p=j
            for j in range(1,n-1):
                if nt>=8:raise ValueError('clip template exceeded capacity')
                out[nt,0]=apex;out[nt,1]=face[p];out[nt,2]=face[(p+j)%n];out[nt,3]=face[(p+j+1)%n]
                inside[nt]=1 if side==0 else 0;nt+=1
    return out,inside,nt

@njit(cache=True)
def clip_counts(t,d,nv,keys,poly,lens):
    out=np.empty(len(t),np.int64)
    for i in range(len(t)):
        _,_,n=clip_child_list(t[i],d,nv,keys,poly,lens);out[i]=n
    return out

@njit(cache=True)
def clip_write(t,code,root,d,nv,keys,poly,lens,start):
    out=np.empty((start[-1],4),np.int32);oc=np.empty(start[-1],np.uint8);orr=np.empty(start[-1],np.int32)
    for i in range(len(t)):
        ch,inc,n=clip_child_list(t[i],d,nv,keys,poly,lens)
        for j in range(n):out[start[i]+j]=ch[j];oc[start[i]+j]=code[i]+6*inc[j];orr[start[i]+j]=root[i]
    return out,oc,orr

def clip(v,t,code,root,d,max_tets=30_000_000):
    nv=len(v);blocks=[]
    for first in range(0,len(t),100_000):blocks.append(clipping_edges(t[first:first+100_000],d,nv))
    keys=np.unique(np.concatenate(blocks));del blocks
    a=keys//nv;b=keys%nv;w=d[a]/(d[a]-d[b]);verts=np.concatenate([v,v[a]+w[:,None]*(v[b]-v[a])])
    counts=clip_counts(t,d,nv,keys,POLY,LENS);starts=np.r_[0,np.cumsum(counts)];total=int(starts[-1])
    if total>max_tets:raise MemoryError(f'cut creates {total} tets > limit {max_tets}')
    tt=np.empty((total,4),np.int32);cc=np.empty(total,np.uint8);rr=np.empty(total,np.int32)
    for first in range(0,len(t),100_000):
        end=min(len(t),first+100_000);sl=slice(first,end);start=np.r_[0,np.cumsum(counts[sl])]
        o,c,r=clip_write(t[sl],code[sl],root[sl],d,nv,keys,POLY,LENS,start)
        tt[starts[first]:starts[end]]=o;cc[starts[first]:starts[end]]=c;rr[starts[first]:starts[end]]=r
    vol=tetvol(verts,tt)
    if np.any(vol==0):raise ValueError(f'cut generated {np.count_nonzero(vol==0)} zero-volume tets; NOT deleted')
    bad=vol<0;tt[bad,0],tt[bad,1]=tt[bad,1].copy(),tt[bad,0].copy()
    return verts,tt,cc,rr


def boundary_faces(t):
    fs=t[:,FACES].reshape(-1,3);ss=np.sort(fs,axis=1)
    key=np.ascontiguousarray(ss).view(np.dtype((np.void,ss.dtype.itemsize*3))).ravel()
    _,idx,count=np.unique(key,return_index=True,return_counts=True)
    if np.any(count>2):raise ValueError('non-manifold tetrahedral faces')
    return fs[idx[count==1]]


def check_mesh(v,t,code,root,original_vol,original_code,all_faces=False):
    vols=tetvol(v,t)
    if np.any(vols<=0.):raise ValueError('non-positive tetrahedron; not allowed to delete it')
    mass=np.bincount(root,weights=vols,minlength=len(original_vol))
    parent_rel=np.max(np.abs(mass-original_vol)/original_vol)
    if parent_rel>2e-9:raise ValueError(f'child volumes do not partition parents: {parent_rel}')
    if not np.array_equal((code-1)%6+1,original_code[root]):raise ValueError('original CT material changed')
    surface=boundary_faces(t[code>6]);q=v[surface]
    # Watertight triangle surface, no zero-area artificial ellipsoid pieces.
    e=np.sort(surface[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
    ek=np.ascontiguousarray(e).view(np.dtype((np.void,e.dtype.itemsize*2))).ravel();_,ecount=np.unique(ek,return_counts=True)
    if not np.all(ecount==2):raise ValueError('ellipsoid surface is open/nonmanifold')
    samples=np.concatenate([q.reshape(-1,3),q.mean(1),(.5*(q[:,0]+q[:,1])),(.5*(q[:,1]+q[:,2])),(.5*(q[:,2]+q[:,0]))])
    maxdist=float(np.max(np.abs(distances(samples))))
    exact=4*np.pi/3*.25*.25*.0125;actual=float(vols[code>6].sum())
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    used,ii=np.unique(surface,return_inverse=True);ee=ii.reshape(-1,3)[:,[[0,1],[1,2],[2,0]]].reshape(-1,2)
    comp=connected_components(coo_matrix((np.ones(len(ee)),(ee[:,0],ee[:,1])),shape=(len(used),len(used))).tocsr(),directed=False)[0]
    if comp!=1:raise ValueError(f'artificial surface has {comp} components')
    outer_count=None
    if all_faces:
        outer=boundary_faces(t);qo=v[outer];external=np.any(np.all(np.abs(qo)<1e-12,axis=1)|np.all(np.abs(qo-1)<1e-12,axis=1),axis=1)
        if not np.all(external):raise ValueError('internal unmatched faces (hanging-node/crack)')
        outer_count=int(len(outer))
    return dict(vertices=len(v),tetrahedra=len(t),ellipsoid_triangles=len(surface),zero_volume_tets=0,
                minimum_tet_volume_over_L3=float(vols.min()),maximum_parent_volume_relative_error=float(parent_rel),
                native_material_inheritance_exact=True,ellipsoid_surface_components=int(comp),surface_closed=True,
                ellipsoid_volume_over_L3=actual,exact_volume_over_L3=exact,ellipsoid_volume_relative_error=actual/exact-1,
                sampled_max_distance_over_full_thickness=maxdist/.025,all_volume_faces_checked=all_faces,
                outer_boundary_triangles=outer_count),surface


def fit(v,t,code,tolerance=5e-5,max_rounds=11,max_tets=30_000_000,verbose=True):
    v=np.ascontiguousarray(v,dtype=np.float64);t=np.ascontiguousarray(t,dtype=np.int32);code=np.asarray(code,np.uint8)
    if len(code)!=len(t) or np.min(code)<1 or np.max(code)>6:raise ValueError('input must be ORIGINAL phase+region codes 1..6, not old fracture labels')
    origcode=code.copy();origvol=np.abs(tetvol(v,t));root=np.arange(len(t),dtype=np.int32)
    if np.any(origvol==0):raise ValueError('original CT source contains zero-volume tets')
    history=[]
    for iteration in range(max_rounds+1):
        timer=time.perf_counter();d=distances(v,snap=tolerance*.1);mark,err,hidden,cross=mark_geometry(v,t,d,tolerance)
        rec=dict(iteration=iteration,vertices=len(v),tetrahedra=len(t),marked=int(mark.sum()),mixed=cross,hidden=hidden,max_sampled_surface_distance_L=err,seconds=time.perf_counter()-timer)
        history.append(rec)
        if verbose:print(json.dumps(rec),flush=True)
        if not np.any(mark):break
        if iteration==max_rounds:raise RuntimeError('geometry tolerance not met; original size/thickness were not changed')
        v,t,code,root=refine(v,t,code,root,mark,max_tets)
    v,t,code,root=clip(v,t,code,root,d,max_tets)
    report,surface=check_mesh(v,t,code,root,origvol,origcode,all_faces=len(t)<2_000_000)
    report['geometry_refinement_history']=history
    if abs(report['ellipsoid_volume_relative_error'])>.005:raise ValueError('volume error exceeds unchanged 0.5% geometry gate')
    if report['sampled_max_distance_over_full_thickness']>.005:raise ValueError('sampled surface distance exceeds unchanged 0.5% full-thickness gate')
    return v,t,code,root,surface,report


def grid(n):
    p=np.linspace(0,1,n+1);x,y,z=np.meshgrid(p,p,p,indexing='ij');v=np.c_[x.ravel(),y.ravel(),z.ravel()]
    ids=np.arange((n+1)**3).reshape(n+1,n+1,n+1);a=ids[:-1,:-1,:-1].ravel();sy=n+1;sx=sy*sy
    corners=np.asarray([0,sx,sx+sy,sy,1,sx+1,sx+sy+1,sy+1]);cube=a[:,None]+corners
    pat=np.asarray([[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6],[0,5,1,6]])
    t=cube[:,pat].reshape(-1,4).astype(np.int32)
    return v,t

# Balanced red subdivision for all-six-edge marked cells. Transition cells
# cone deterministic face subdivisions to the original cell centroid.
# This avoids the long, oblique children produced by arbitrary sequential
# edge bisection and still guarantees identical triangulations on shared faces.
@njit(cache=True)
def triangle_subdivision(ids,mid):
    # ids=[a,b,c]; mid=[ab,bc,ca], -1 if unsplit.
    out=np.empty((4,3),np.int64);nm=0
    for k in range(3):
        if mid[k]>=0:nm+=1
    a,b,c=ids;ab,bc,ca=mid
    if nm==0:out[0]=ids;return out,1
    if nm==3:
        out[0]=np.array([a,ab,ca]);out[1]=np.array([ab,b,bc]);out[2]=np.array([ca,bc,c]);out[3]=np.array([ab,bc,ca]);return out,4
    if nm==1:
        for k in range(3):
            if mid[k]>=0:
                aa=ids[k];bb=ids[(k+1)%3];cc=ids[(k+2)%3];mm=mid[k]
                out[0]=np.array([aa,mm,cc]);out[1]=np.array([mm,bb,cc]);return out,2
    # exactly two marked edges, rotate so unmarked is ab.
    for k in range(3):
        if mid[k]<0:
            aa=ids[k];bb=ids[(k+1)%3];cc=ids[(k+2)%3];mbc=mid[(k+1)%3];mca=mid[(k+2)%3]
            out[0]=np.array([cc,mca,mbc]);quad=np.array([aa,bb,mbc,mca]);p=np.argmin(quad)
            out[1]=np.array([quad[p],quad[(p+1)%4],quad[(p+2)%4]])
            out[2]=np.array([quad[p],quad[(p+2)%4],quad[(p+3)%4]]);return out,3
    raise ValueError('invalid triangle subdivision')

@njit(cache=True)
def edge_mid_matrix(t,nv,chosen):
    mids=np.full((len(t),6),-1,np.int64);counts=np.empty(len(t),np.int64);partial=np.zeros(len(t),np.int64)
    for i in range(len(t)):
        n=0
        for e in range(6):
            a=t[i,EDGES[e,0]];b=t[i,EDGES[e,1]]
            if a>b:a,b=b,a
            key=a*nv+b;ix=np.searchsorted(chosen,key)
            if ix<len(chosen) and chosen[ix]==key:mids[i,e]=nv+ix;n+=1
        if n==0:counts[i]=1
        elif n==6:counts[i]=8
        else:
            partial[i]=1;counts[i]=0
            for f in range(4):
                nface=0
                for e in range(6):
                    a=EDGES[e,0];b=EDGES[e,1]
                    if a!=f and b!=f and mids[i,e]>=0:nface+=1
                counts[i]+=nface+1
    return mids,counts,partial

@njit(cache=True)
def write_red_transition(v,t,code,root,mid,start,centers):
    out=np.empty((start[-1],4),np.int32);oc=np.empty(start[-1],np.uint8);orr=np.empty(start[-1],np.int32)
    lookup=np.full((4,4),-1,np.int64)
    for e in range(6):lookup[EDGES[e,0],EDGES[e,1]]=e;lookup[EDGES[e,1],EDGES[e,0]]=e
    for i in range(len(t)):
        nn=start[i+1]-start[i];j=start[i]
        if nn==1:out[j]=t[i]
        elif centers[i]<0:
            for k in range(4):
                out[j+k,0]=t[i,k];p=1
                for e in range(6):
                    if EDGES[e,0]==k or EDGES[e,1]==k:out[j+k,p]=mid[i,e];p+=1
            diagonals=np.asarray([[0,5],[1,4],[2,3]])
            rings=np.asarray([[1,2,4,3],[0,2,5,3],[0,1,5,4]])
            dd=np.empty(3)
            for k in range(3):
                p=v[mid[i,diagonals[k,0]]]-v[mid[i,diagonals[k,1]]];dd[k]=np.dot(p,p)
            k=np.argmin(dd);d1=mid[i,diagonals[k,0]];d2=mid[i,diagonals[k,1]]
            for r in range(4):out[j+4+r]=np.array([d1,d2,mid[i,rings[k,r]],mid[i,rings[k,(r+1)%4]]])
        else:
            for f in range(4):
                fids=FACES[f];mm=np.empty(3,np.int64)
                for k in range(3):mm[k]=mid[i,lookup[fids[k],fids[(k+1)%3]]]
                fs,n=triangle_subdivision(t[i,fids],mm)
                for k in range(n):out[j]=np.array([centers[i],fs[k,0],fs[k,1],fs[k,2]]);j+=1
        for k in range(start[i],start[i+1]):oc[k]=code[i];orr[k]=root[i]
    return out,oc,orr

def refine(v,t,code,root,mark,max_tets=30_000_000):
    nv=len(v);chosen=np.unique(edge_keys(t[mark],nv))
    if len(chosen)==0:return v,t,code,root
    counts=np.empty(len(t),np.int64);partial=np.empty(len(t),np.int64)
    for first in range(0,len(t),100_000):
        sl=slice(first,first+100_000);_,counts[sl],partial[sl]=edge_mid_matrix(t[sl],nv,chosen)
    start=np.r_[0,np.cumsum(counts)];total=int(start[-1])
    if total>max_tets:raise MemoryError(f'local refinement would create {total} tets, limit {max_tets}')
    centers=np.full(len(t),-1,np.int64);pi=np.flatnonzero(partial);centers[pi]=nv+len(chosen)+np.arange(len(pi))
    verts=np.concatenate([v,.5*(v[chosen//nv]+v[chosen%nv]),v[t[pi]].mean(1)])
    tt=np.empty((total,4),np.int32);cc=np.empty(total,np.uint8);rr=np.empty(total,np.int32)
    for first in range(0,len(t),100_000):
        end=min(len(t),first+100_000);sl=slice(first,end);mid,_,_=edge_mid_matrix(t[sl],nv,chosen)
        o,c,r=write_red_transition(verts,t[sl],code[sl],root[sl],mid,np.r_[0,np.cumsum(counts[sl])],centers[sl])
        tt[start[first]:start[end]]=o;cc[start[first]:start[end]]=c;rr[start[first]:start[end]]=r
    return verts,tt,cc,rr
