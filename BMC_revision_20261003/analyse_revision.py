"""Recompute reviewer-requested summaries from immutable completed outputs."""
from pathlib import Path
import hashlib
import json
import sys
import platform
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
MASTER = ROOT / 'results/bmc_master'
SPLITS = ['split_random', 'split_scaffold_heldout', 'split_cell_scaffold_heldout']

def interval(values):
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return dict(n_seeds=0, mean=np.nan, sd=np.nan, low=np.nan, high=np.nan)
    boot = np.random.default_rng(20261003).choice(x, (10000, len(x)), replace=True).mean(axis=1)
    return dict(n_seeds=len(x), mean=float(x.mean()), sd=float(x.std(ddof=1)) if len(x)>1 else np.nan,
                low=float(np.quantile(boot,.025)), high=float(np.quantile(boot,.975)))

def main():
    p = pd.read_csv(MASTER/'master_performance_seed_level.csv')
    p = p[p.split.isin(SPLITS) & ~p.model.str.contains('prnet', case=False)].copy()
    baseline = p[p.model.eq('global_train_mean')][['dataset','seed','split','metric','value']].rename(columns={'value':'baseline_value'})
    adj = p.merge(baseline, on=['dataset','seed','split','metric'], validate='many_to_one')
    adj['gain_over_global_mean'] = adj.value-adj.baseline_value
    adj.to_csv(OUT/'baseline_adjusted_seed_level.csv',index=False)
    summary=[]
    for key,g in adj.groupby(['dataset','model','split','metric']):
        summary.append(dict(zip(['dataset','model','split','metric'],key)) | interval(g.gain_over_global_mean))
    pd.DataFrame(summary).to_csv(OUT/'baseline_adjusted_summary.csv',index=False)

    # Score both fitted models on the exact same record intersection. All
    # selected records were held out by both original fits; no refitting occurs.
    per=[]
    for seed in range(1,101):
        base=ROOT/f'results/openproblems_multiseed_100_de/seed_{seed:03d}'
        for name in ['mean_baselines/simple_baseline_per_row.csv','ridge_baselines/ridge_baseline_per_row.csv']:
            f=pd.read_csv(base/name)
            f['seed']=seed
            per.append(f)
    op=pd.concat(per,ignore_index=True)
    sci=pd.read_csv(ROOT/'results/sciplex3_24h_top2000_multiseed_30/baseline_per_row_all.csv')
    intersections=[]
    for dataset,f in [('OpenProblems',op),('Sci-Plex 3',sci)]:
        for (seed,model),g in f.groupby(['seed','baseline']):
            a=g[g.split.eq('split_random')]
            b=g[g.split.eq('split_cell_scaffold_heldout')]
            merged=a.merge(b,on='sample_id',suffixes=('_random','_joint'),validate='one_to_one')
            metrics=[m for m in ['rowwise_pearson','de_rowwise_pearson','rmse','de_rmse'] if m in g.columns]
            for metric in metrics:
                pair=merged[[metric+'_random',metric+'_joint']].dropna()
                intersections.append(dict(dataset=dataset,seed=seed,model=model,metric=metric,n_common_records=len(pair),
                    random_mean=pair.iloc[:,0].mean(),joint_mean=pair.iloc[:,1].mean(),
                    difference=(pair.iloc[:,0]-pair.iloc[:,1]).mean()))
    inter=pd.DataFrame(intersections)
    inter.to_csv(OUT/'common_test_records_seed_level.csv',index=False)
    summaries=[]
    for key,g in inter.groupby(['dataset','model','metric']):
        usable=g[g.n_common_records.gt(0)]
        summaries.append(dict(zip(['dataset','model','metric'],key)) | interval(usable.difference) |
            dict(mean_common_records=usable.n_common_records.mean(),min_common_records=usable.n_common_records.min(),
                 max_common_records=usable.n_common_records.max(),random_mean=usable.random_mean.mean(),joint_mean=usable.joint_mean.mean()))
    pd.DataFrame(summaries).to_csv(OUT/'common_test_records_summary.csv',index=False)

    # Chemical and annotation counts are computed on the full final cohort,
    # separately from seed-dependent annotated-test denominators.
    sys.path.insert(0,str(OUT/'runtime'))
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdFingerprintGenerator
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024,includeChirality=False)
    moa=pd.read_csv(ROOT/'data/raw/openproblems_neurips2023/moa_annotations.csv')
    moa_map=moa.dropna(subset=['sm_name','moa']).assign(key=lambda x:x.sm_name.astype(str).str.strip().str.lower()).groupby('key').moa.agg(lambda x:set(x.astype(str))).to_dict()
    op_meta=pd.read_csv(ROOT/'results/openproblems_multiseed_100_de/seed_001/splits/candidate_splits.csv')
    op_meta=op_meta[op_meta.condition.eq('treated')]
    sci_meta=pd.read_csv(ROOT/'metadata/sciplex3_24h_top2000_response_metadata.csv').rename(columns={'perturbation':'drug_name','canonical_smiles':'smiles'})
    rows=[]
    for dataset,meta in [('OpenProblems',op_meta),('Sci-Plex 3',sci_meta)]:
        for drug,g in meta.groupby('drug_name',dropna=False):
            smi=g.smiles.iloc[0]
            missing=pd.isna(smi) or not str(smi).strip()
            mol=None if missing else Chem.MolFromSmiles(str(smi))
            labels=moa_map.get(str(drug).strip().lower(),set()) if dataset=='OpenProblems' else set()
            rows.append(dict(dataset=dataset,drug_name=drug,n_records=len(g),missing_smiles=missing,invalid_smiles=not missing and mol is None,
                zero_fingerprint=mol is None or gen.GetFingerprint(mol).GetNumOnBits()==0,
                missing_scaffold=pd.isna(g.scaffold.iloc[0]) or not str(g.scaffold.iloc[0]).strip(),
                n_distinct_smiles=g.smiles.nunique(dropna=True),moa_annotated=bool(labels),n_moa_labels=len(labels)))
    chem=pd.DataFrame(rows)
    chem.to_csv(OUT/'chemical_annotation_coverage_by_drug.csv',index=False)
    counts=[]
    for dataset,g in chem.groupby('dataset'):
        row=dict(dataset=dataset,n_drugs=len(g),n_records=int(g.n_records.sum()))
        source_meta=op_meta if dataset=='OpenProblems' else sci_meta
        row['unique_scaffolds']=int(source_meta.scaffold.replace('',np.nan).nunique())
        for field in ['missing_smiles','invalid_smiles','zero_fingerprint','missing_scaffold','moa_annotated']:
            row[field+'_drugs']=int(g[field].sum())
            row[field+'_records']=int(g.loc[g[field],'n_records'].sum())
        counts.append(row)
    pd.DataFrame(counts).to_csv(OUT/'chemical_annotation_coverage_summary.csv',index=False)
    inputs=[MASTER/'master_performance_seed_level.csv',ROOT/'data/raw/openproblems_neurips2023/moa_annotations.csv',
        ROOT/'results/sciplex3_24h_top2000_multiseed_30/baseline_per_row_all.csv']
    manifest={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs}
    (OUT/'analysis_provenance.json').write_text(json.dumps(dict(inputs=manifest,bootstrap_resamples=10000,bootstrap_seed=20261003,
        uncertainty='conditional split-construction sensitivity on the fixed observed cohorts',python=sys.version,
        platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__,rdkit=rdBase.rdkitVersion),indent=2))
    print('Completed baseline-adjusted scores, common-test-record contrasts, and representation/annotation coverage.')
    print(pd.DataFrame(counts).to_string(index=False))

if __name__=='__main__':
    main()
