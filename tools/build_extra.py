"""Tambah data ke WebGIS: raster DC & bahaya, riwayat terbakar MODIS, titik api, ENSO.
Pakai: python build_extra.py <grid_lama.geojson> <dir_DC> <dir_bahaya> <MODIS_burn.tif> <sebaran-hotspot.txt> <ENSO.txt> <out_dir>
Butuh: numpy pandas tifffile pillow matplotlib
Catatan: sel grid 5 km bergeser setengah piksel terhadap raster; sel diambil dari piksel yang titik
pojok kiri-atasnya = pusat piksel (korelasi bahaya~kerentanan tertinggi dibanding 3 konvensi lain)."""
import json, sys, math, warnings, collections
import numpy as np, pandas as pd, tifffile
from PIL import Image
from matplotlib.path import Path
warnings.filterwarnings('ignore')
A=6378137.0; E2=(1/298.257223563)*(2-1/298.257223563); EP2=E2/(1-E2); K0=0.9996
def fwd(lon,lat):
    lon=np.asarray(lon,float); lat=np.asarray(lat,float)
    phi=np.radians(lat); lam=np.radians(lon-99.0)
    N=A/np.sqrt(1-E2*np.sin(phi)**2); T=np.tan(phi)**2; C=EP2*np.cos(phi)**2; Aa=np.cos(phi)*lam
    e4=E2*E2; e6=E2**3
    M=A*((1-E2/4-3*e4/64-5*e6/256)*phi-(3*E2/8+3*e4/32+45*e6/1024)*np.sin(2*phi)+(15*e4/256+45*e6/1024)*np.sin(4*phi)-(35*e6/3072)*np.sin(6*phi))
    x=K0*N*(Aa+(1-T+C)*Aa**3/6+(5-18*T+T*T+72*C-58*EP2)*Aa**5/120)+500000
    y=K0*(M+N*np.tan(phi)*(Aa**2/2+(5-T+9*C+4*C*C)*Aa**4/24+(61-58*T+T*T+600*C-330*EP2)*Aa**6/720))
    return x,y
X0,Y0,PX=614281.0,325504.1846,5000.0
MON=[('Jan','Januari'),('Feb','Februari'),('Mar','Maret'),('Apr','April'),('Mei','Mei'),('Jun','Juni'),('Jul','Juli'),('Ags','Agustus'),('Sep','September')]
def lattice(x,y):
    return np.floor((x-(X0+PX/2))/PX).astype(int), np.floor(((Y0-PX/2)-y)/PX).astype(int)
def main(grid,dcdir,bhdir,burn,hot,enso,out):
    gj=json.load(open(grid)); F=gj['features']
    # indeks sel -> (kolom,baris) lattice
    cr=[]
    for f in F:
        pts=np.array([[x,y] for poly in f['geometry']['coordinates'] for ring in poly[:1] for x,y in ring])
        X,Y=fwd(pts[:,0],pts[:,1]); cr.append((int((X.min()+X.max())/2),int((Y.min()+Y.max())/2)))
    cr=np.array(cr); k,j=lattice(cr[:,0].astype(float),cr[:,1].astype(float))
    cell_of={(int(a),int(b)):i for i,(a,b) in enumerate(zip(k,j))}
    assert len(cell_of)==len(F),'sel ganda di lattice'
    # raster bulanan
    for mi,(s,l) in enumerate(MON,1):
        d=tifffile.imread(f'{dcdir}/Reclass_{s}.tif').astype(float)
        h=tifffile.imread(f'{bhdir}/Peta_Bahaya_{l}.tif').astype(float)
        for i,f in enumerate(F):
            c,r=int(k[i]),int(j[i]); p=f['properties']
            if 0<=c<85 and 0<=r<91:
                dv,hv=d[r,c],h[r,c]
                p['dc_%02d'%mi]=int(dv) if dv>0 else None
                p['bhy_%02d'%mi]=int(hv) if hv>0 else None
            else:
                p['dc_%02d'%mi]=None; p['bhy_%02d'%mi]=None
    # riwayat terbakar MODIS (500 m): persen piksel terbakar + tahun terakhir
    b=np.array(Image.open(burn)); H,W=b.shape
    gx0,gy0,gs=100.02291531028813,2.5467238304788435,0.004491576420597608
    lon=gx0+(np.arange(W)+.5)*gs; lat=gy0-(np.arange(H)+.5)*gs
    LON,LAT=np.meshgrid(lon,lat)
    for f in F:
        p=f['properties']; tot=0; nb=0; last=0
        for poly in f['geometry']['coordinates']:
            ring=np.array(poly[0]); x0,y0=ring.min(0); x1,y1=ring.max(0)
            c0,c1=max(0,int((x0-gx0)/gs)-1),min(W,int((x1-gx0)/gs)+2); r0,r1=max(0,int((gy0-y1)/gs)-1),min(H,int((gy0-y0)/gs)+2)
            if c1<=c0 or r1<=r0: continue
            sx=LON[r0:r1,c0:c1].ravel(); sy=LAT[r0:r1,c0:c1].ravel(); pts=np.c_[sx,sy]
            m=Path(ring).contains_points(pts)
            for hole in poly[1:]: m&=~Path(np.array(hole)).contains_points(pts)
            vals=b[r0:r1,c0:c1].ravel()[m]; tot+=len(vals); nb+=int((vals>0).sum())
            if (vals>0).any(): last=max(last,int(vals.max()))
        p['bk']=round(100*nb/tot,1) if tot else None; p['by']=last or None
    for f in F:
        for key in ('sip_bin','pop'): f['properties'].pop(key,None)
    json.dump({'type':'FeatureCollection','features':F},open(out+'/grid.geojson','w'),separators=(',',':'))
    # titik api
    d=pd.read_csv(hot,skipinitialspace=True); d.columns=[c.strip() for c in d.columns]
    d['t']=pd.to_datetime(d['Tanggal'],format='%d-%m-%Y',errors='coerce'); d=d.dropna(subset=['t']).drop_duplicates()
    names={f['properties']['kab'].lower():f['properties']['kab'] for f in F}; names['dumai']=names['kota dumai']
    d['kab']=d['Kab Kota'].str.strip().str.lower().map(names)
    kabs=sorted(set(names.values())); kidx={n:i for i,n in enumerate(kabs)}
    d=d.dropna(subset=['kab'])
    d['Kecamatan']=d['Kecamatan'].fillna('-').astype(str); d['Desa']=d['Desa'].fillna('-').astype(str)
    X,Y=fwd(d['Longitude'].values,d['Latitude'].values); kk,jj=lattice(X,Y)
    g=np.array([F[cell_of[(a,b)]]['properties']['g'] if (a,b) in cell_of else -1 for a,b in zip(kk,jj)])
    kec=sorted(d['Kecamatan'].str.strip().unique()); kcid={n:i for i,n in enumerate(kec)}
    des=sorted(d['Desa'].astype(str).str.strip().unique()); dsid={n:i for i,n in enumerate(des)}
    d=d.assign(g=g)
    d=d.sort_values('t')
    sat={'NASA-MODIS':0,'NASA-NOAA20':1,'NASA-NOAA21':2}; conf={'Low':0,'Medium':1,'High':2}
    hr=d['Waktu'].str.extract(r'(\d+):(\d+)').astype(float)
    cols={'x':[round(v,4) for v in d['Longitude']],'y':[round(v,4) for v in d['Latitude']],
      'th':[int(v) for v in d.t.dt.year],'bl':[int(v) for v in d.t.dt.month],'hr':[int(v) for v in d.t.dt.day],
      'k':[kidx[v] for v in d.kab],'kc':[kcid[v.strip()] for v in d['Kecamatan']],'ds':[dsid[str(v).strip()] for v in d['Desa']],
      'c':[conf[v] for v in d['Confidence'].str.strip()],'s':[sat[v] for v in d['Satelit'].str.strip()],'g':[int(v) for v in d.g],
      'jm':[int(v) if v==v else -1 for v in hr[0]]}
    json.dump({'kab':kabs,'kec':kec,'desa':des,'sat':['MODIS','NOAA-20','NOAA-21'],'conf':['Rendah','Sedang','Tinggi'],
               'tgl_akhir':str(d.t.max().date()),'n':len(d),**cols},open(out+'/titik_api.json','w'),separators=(',',':'))
    # ENSO: rata-rata anomali ONI musim ASO per tahun
    e=pd.read_csv(enso,sep=r'\s+'); res={}
    for y,gp in e.groupby('YR'):
        if 2014<=y<=2026:
            v=gp[gp.SEAS.isin(['JAS','ASO'])].ANOM
            if len(v): res[int(y)]=round(float(v.mean()),2)
    json.dump(res,open(out+'/enso.json','w'),separators=(',',':'))
    print('sel',len(F),'titik',len(d),'tanpa grid',int((g<0).sum()),'enso',res)
if __name__=='__main__': main(*sys.argv[1:])
