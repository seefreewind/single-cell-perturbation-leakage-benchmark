"""Train-only top-2000 sensitivity using the original full-gene Sci-Plex cells."""
from pathlib import Path
import sys
import json
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
OUT=Path(__file__).resolve().parent
ROOT=OUT.parent
sys.path.insert(0,str(OUT/'runtime'))
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, linalg
from threadpoolctl import threadpool_limits
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

def read_col(group,name):
    v=group[name]
    if isinstance(v,h5py.Group):
        cats=v['categories'].asstr()[:] if h5py.check_string_dtype(v['categories'].dtype) else v['categories'][:]
        codes=v['codes'][:]
        return np.array([cats[i] if i>=0 else None for i in codes])
    return v.asstr()[:] if h5py.check_string_dtype(v.dtype) else v[:]

def pearson(y,p):
    a=y-y.mean(1,keepdims=True); b=p-p.mean(1,keepdims=True)
    den=np.sqrt((a*a).sum(1)*(b*b).sum(1))
    return np.divide((a*b).sum(1),den,out=np.full(len(y),np.nan),where=den>0)

def main():
    meta=pd.read_csv(ROOT/'metadata/sciplex3_24h_top2000_response_metadata.csv')
    f=h5py.File(ROOT/'data/raw/scperturb/SrivatsanTrapnell2020_sciplex3.h5ad','r')
    obs=f['obs']; names=['perturbation','cell_line','dose_value','time','ncounts']
    cells=pd.DataFrame({k:read_col(obs,k) for k in names})
    keys=['perturbation','cell_line','dose_value','time']
    record_lookup={tuple(row):i for i,row in enumerate(meta[keys].itertuples(index=False,name=None))}
    lines=['A549','K562','MCF7']; n=len(meta)
    ids=np.array([record_lookup.get(tuple(row),-1) for row in cells[keys].itertuples(index=False,name=None)],dtype=np.int32)
    for j,line in enumerate(lines):
        mask=cells.perturbation.astype(str).str.lower().eq('control') & cells.cell_line.eq(line) & cells.time.eq(24)
        ids[mask.to_numpy()]=n+j
    ids[cells.ncounts.isna().to_numpy()]=-1
    counts=np.bincount(ids[ids>=0],minlength=n+3)
    assert np.array_equal(counts[:n],meta.n_cells.to_numpy()),'Raw-cell membership disagrees with fixed pseudobulk metadata'
    ng=int(f['X'].attrs['shape'][1])
    raw=np.memmap(OUT/'raw_gene_sums.dat',dtype='float32',mode='w+',shape=(n+3,ng)); raw[:]=0
    norm=np.memmap(OUT/'normalized_gene_sums.dat',dtype='float32',mode='w+',shape=(n+3,ng)); norm[:]=0
    ptr=f['X/indptr'][:]; indices=f['X/indices']; data=f['X/data']
    for start in range(0,len(cells),5000):
        end=min(start+5000,len(cells)); keep=ids[start:end]>=0
        if not keep.any(): continue
        a,b=int(ptr[start]),int(ptr[end]); x=sparse.csr_matrix((data[a:b].astype(np.float32),indices[a:b],ptr[start:end+1]-a),shape=(end-start,ng))
        x=x[keep]; gi=ids[start:end][keep]
        assign=sparse.csr_matrix((np.ones(len(gi),dtype=np.float32),(gi,np.arange(len(gi)))),shape=(n+3,len(gi)))
        block=(assign@x).tocoo(); np.add.at(raw,(block.row,block.col),block.data)
        denom=cells.ncounts.to_numpy()[start:end][keep].astype(float); denom[denom<=0]=1
        x.data=np.log1p(x.data*np.repeat(10000/denom,np.diff(x.indptr))).astype(np.float32)
        block=(assign@x).tocoo(); np.add.at(norm,(block.row,block.col),block.data)
        print(f'Aggregated cells through {end}/{len(cells)}',flush=True)
    raw.flush(); norm.flush()
    ensembl=read_col(f['var'],'ensembl_id'); symbols=read_col(f['var'],'gene_symbol')
    f.close()
    control_ids=np.array([n+lines.index(line) for line in meta.cell_line])
    # Validate aggregation against the existing fixed-panel response targets.
    old=h5py.File(ROOT/'data/processed/sciplex3/sciplex3_24h_top2000_response.h5ad','r')
    old_gene=read_col(old['var'],'ensembl_id'); lookup={str(v):i for i,v in enumerate(ensembl)}
    fixed=np.array([lookup[str(v)] for v in old_gene]); y_old=old['X'][:]
    y_check=norm[:n,fixed]/counts[:n,None]-norm[control_ids[:,None],fixed[None,:]]/counts[control_ids,None]
    max_error=float(np.max(np.abs(y_check-y_old)))
    assert max_error<1e-4,f'Pseudobulk reconstruction error {max_error}'
    old.close()
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024,includeChirality=False)
    fps=np.vstack([np.asarray(gen.GetFingerprint(Chem.MolFromSmiles(str(s))),dtype=float) for s in meta.canonical_smiles])
    splitroot=ROOT/'results/sciplex3_24h_top2000_multiseed_30'
    results=[]; gene_rows=[]
    with threadpool_limits(limits=1):
        for seed in range(1,31):
            manifest=pd.read_csv(splitroot/f'seed_{seed:03d}/splits/candidate_splits.csv').set_index('sample_id').loc[meta.sample_id]
            for split in ['split_random','split_scaffold_heldout','split_cell_scaffold_heldout']:
                tr=np.flatnonzero(manifest[split].eq('train').to_numpy()); te=np.flatnonzero(manifest[split].eq('test').to_numpy())
                sums=np.asarray(raw[tr]).sum(0,dtype=np.float64)
                selected=np.argsort(-sums,kind='stable')[:2000]
                gene_rows.extend(dict(seed=seed,split=split,rank=j+1,ensembl_id=str(ensembl[i]),gene_symbol=str(symbols[i]),train_raw_count=float(sums[i])) for j,i in enumerate(selected))
                y=np.asarray(norm[:n,selected])/counts[:n,None]-np.asarray(norm[control_ids[:,None],selected[None,:]])/counts[control_ids,None]
                known=sorted(meta.iloc[tr].cell_line.unique()); cell=np.column_stack([meta.cell_line.eq(v).to_numpy(float) for v in known])
                dose=np.log10(meta.dose_value.to_numpy(float)); dose=(dose-dose[tr].mean())/dose[tr].std()
                x=np.column_stack([np.ones(n),cell,dose,fps]); xt=x[tr]; penalty=np.eye(x.shape[1])*10; penalty[0,0]=0
                coef=linalg.solve(xt.T@xt+penalty,xt.T@y[tr],assume_a='pos')
                preds={'global_train_mean':np.repeat(y[tr].mean(0,keepdims=True),len(te),0),'ridge_cell_dose_drug_fp':x[te]@coef}
                for model,pred in preds.items():
                    score=pearson(y[te],pred)
                    results.append(dict(seed=seed,split=split,model=model,n_train=len(tr),n_test=len(te),mean_rowwise_pearson=float(np.nanmean(score)),rmse=float(np.sqrt(((y[te]-pred)**2).mean())),n_genes=2000))
            print(f'Completed train-only gene selection sensitivity seed {seed}',flush=True)
    pd.DataFrame(results).to_csv(OUT/'train_only_gene_selection_seed_level.csv',index=False)
    pd.DataFrame(gene_rows).to_csv(OUT/'train_only_selected_genes.csv',index=False)
    (OUT/'train_only_gene_selection_validation.json').write_text(json.dumps(dict(max_fixed_panel_reconstruction_error=max_error,
        raw_gene_count=ng,n_treated_records=n,n_seeds=30,selection='top 2000 by summed raw counts in manifest training treated records only',
        control_definition='same public matched cell-line control references as the original response target; sensitivity addresses gene-selection frequencies only'),indent=2))
    del raw,norm
    (OUT/'raw_gene_sums.dat').unlink(); (OUT/'normalized_gene_sums.dat').unlink()
    print('Completed all 30 seeds and three split settings',flush=True)

if __name__=='__main__': main()
