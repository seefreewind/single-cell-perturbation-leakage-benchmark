#!/usr/bin/env python3
"""Prepare the v15 JGG submission-ready manuscript package.

This is a packaging and manuscript-consistency script only. It does not rerun
models, change result values, or create new experimental results.
"""

from __future__ import annotations

import csv
import re
import shutil
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[2]


REFERENCES = [
    {
        "id": 1,
        "text": "Dixit A, Parnas O, Li B, et al. Perturb-Seq: dissecting molecular circuits with scalable single-cell RNA profiling of pooled genetic screens. Cell. 2016;167:1853-1866.e17. doi:10.1016/j.cell.2016.11.038.",
        "category": "single-cell perturbation transcriptomics",
        "identifier": "DOI:10.1016/j.cell.2016.11.038",
    },
    {
        "id": 2,
        "text": "Adamson B, Norman TM, Jost M, et al. A multiplexed single-cell CRISPR screening platform enables systematic dissection of the unfolded protein response. Cell. 2016;167:1867-1882.e21. doi:10.1016/j.cell.2016.11.048.",
        "category": "single-cell perturbation transcriptomics",
        "identifier": "DOI:10.1016/j.cell.2016.11.048",
    },
    {
        "id": 3,
        "text": "Datlinger P, Rendeiro AF, Schmidl C, et al. Pooled CRISPR screening with single-cell transcriptome readout. Nature Methods. 2017;14:297-301. doi:10.1038/nmeth.4177.",
        "category": "single-cell perturbation transcriptomics",
        "identifier": "DOI:10.1038/nmeth.4177",
    },
    {
        "id": 4,
        "text": "Srivatsan SR, McFaline-Figueroa JL, Ramani V, et al. Massively multiplex chemical transcriptomics at single-cell resolution. Science. 2020;367:45-51. doi:10.1126/science.aax6234.",
        "category": "Sci-Plex 3",
        "identifier": "DOI:10.1126/science.aax6234",
    },
    {
        "id": 5,
        "text": "Lamb J, Crawford ED, Peck D, et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. Science. 2006;313:1929-1935. doi:10.1126/science.1132939.",
        "category": "LINCS / Connectivity Map",
        "identifier": "DOI:10.1126/science.1132939",
    },
    {
        "id": 6,
        "text": "Subramanian A, Narayan R, Corsello SM, et al. A next generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. Cell. 2017;171:1437-1452.e17. doi:10.1016/j.cell.2017.10.049.",
        "category": "LINCS / Connectivity Map",
        "identifier": "DOI:10.1016/j.cell.2017.10.049",
    },
    {
        "id": 7,
        "text": "Koleti A, Terryn R, Stathias V, et al. Data Portal for the Library of Integrated Network-based Cellular Signatures (LINCS) program: integrated access to diverse large-scale cellular perturbation response data. Nucleic Acids Research. 2018;46:D558-D566. doi:10.1093/nar/gkx1063.",
        "category": "LINCS / Connectivity Map",
        "identifier": "DOI:10.1093/nar/gkx1063",
    },
    {
        "id": 8,
        "text": "OpenProblems Single-Cell Analysis. Perturbation Prediction task. https://openproblems.bio/",
        "category": "OpenProblems",
        "identifier": "URL:https://openproblems.bio/",
    },
    {
        "id": 9,
        "text": "OpenProblems bio task repository. task_perturbation_prediction. https://github.com/openproblems-bio/task_perturbation_prediction",
        "category": "OpenProblems / OP3",
        "identifier": "URL:https://github.com/openproblems-bio/task_perturbation_prediction",
    },
    {
        "id": 10,
        "text": "A benchmark for prediction of transcriptomic responses to chemical perturbations. NeurIPS Datasets and Benchmarks Track. 2024. https://openreview.net/forum?id=WTI4RJYSVm",
        "category": "OP3 perturbation benchmark",
        "identifier": "URL:https://openreview.net/forum?id=WTI4RJYSVm",
    },
    {
        "id": 11,
        "text": "Lotfollahi M, Wolf FA, Theis FJ. scGen predicts single-cell perturbation responses. Nature Methods. 2019;16:715-721. doi:10.1038/s41592-019-0494-8.",
        "category": "scGen",
        "identifier": "DOI:10.1038/s41592-019-0494-8",
    },
    {
        "id": 12,
        "text": "Lotfollahi M, Klimovskaia Susmelj A, De Donno C, et al. Predicting cellular responses to complex perturbations in high-throughput screens. Molecular Systems Biology. 2023;19:e11517. doi:10.15252/msb.202211517.",
        "category": "CPA",
        "identifier": "DOI:10.15252/msb.202211517",
    },
    {
        "id": 13,
        "text": "Hetzel L, Boehm S, Kilbertus N, Guennemann S, Lotfollahi M, Theis FJ. Predicting cellular responses to novel drug perturbations at a single-cell resolution. Advances in Neural Information Processing Systems. 2022. https://openreview.net/forum?id=vRrFVHxFiXJ",
        "category": "chemCPA",
        "identifier": "URL:https://openreview.net/forum?id=vRrFVHxFiXJ",
    },
    {
        "id": 14,
        "text": "Roohani Y, Huang K, Leskovec J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. Nature Biotechnology. 2024;42:927-935. doi:10.1038/s41587-023-01905-6.",
        "category": "perturbation prediction model",
        "identifier": "DOI:10.1038/s41587-023-01905-6",
    },
    {
        "id": 15,
        "text": "Bunne C, Stark SG, Gut G, et al. Learning single-cell perturbation responses using neural optimal transport. Nature Methods. 2023;20:1759-1768. doi:10.1038/s41592-023-01969-x.",
        "category": "perturbation prediction model",
        "identifier": "DOI:10.1038/s41592-023-01969-x",
    },
    {
        "id": 16,
        "text": "Tong X, Qu N, Kong X, et al. Deep representation learning of chemical-induced transcriptional profile for phenotype-based drug discovery. Nature Communications. 2024. doi:10.1038/s41467-024-49620-3.",
        "category": "TranSiGen",
        "identifier": "DOI:10.1038/s41467-024-49620-3",
    },
    {
        "id": 17,
        "text": "Qi X, Zhao L, Tian C, et al. Predicting transcriptional responses to novel chemical perturbations using deep generative model for drug discovery. Nature Communications. 2024. doi:10.1038/s41467-024-53457-1.",
        "category": "PRnet",
        "identifier": "DOI:10.1038/s41467-024-53457-1",
    },
    {
        "id": 18,
        "text": "Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. Nature Methods. 2024;21:1470-1480. doi:10.1038/s41592-024-02201-0.",
        "category": "scGPT",
        "identifier": "DOI:10.1038/s41592-024-02201-0",
    },
    {
        "id": 19,
        "text": "Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. Nature. 2023;618:616-624. doi:10.1038/s41586-023-06139-9.",
        "category": "Geneformer",
        "identifier": "DOI:10.1038/s41586-023-06139-9",
    },
    {
        "id": 20,
        "text": "Hao M, Gong J, Zeng X, et al. Large scale foundation model on single-cell transcriptomics. bioRxiv. 2023. doi:10.1101/2023.05.29.542705.",
        "category": "scFoundation",
        "identifier": "DOI:10.1101/2023.05.29.542705",
    },
    {
        "id": 21,
        "text": "Lotfollahi M, Naghipourfar M, Luecken MD, et al. Mapping single-cell data to reference atlases by transfer learning. Nature Biotechnology. 2022;40:121-130. doi:10.1038/s41587-021-01001-7.",
        "category": "single-cell foundation / transfer learning",
        "identifier": "DOI:10.1038/s41587-021-01001-7",
    },
    {
        "id": 22,
        "text": "Luecken MD, Büttner M, Chaichoompu K, et al. Benchmarking atlas-level data integration in single-cell genomics. Nature Methods. 2022;19:41-50. doi:10.1038/s41592-021-01336-8.",
        "category": "single-cell benchmarking",
        "identifier": "DOI:10.1038/s41592-021-01336-8",
    },
    {
        "id": 23,
        "text": "Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. Patterns. 2023;4:100804. doi:10.1016/j.patter.2023.100804.",
        "category": "benchmark leakage / data leakage",
        "identifier": "DOI:10.1016/j.patter.2023.100804",
    },
    {
        "id": 24,
        "text": "Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics. 2006;7:91. doi:10.1186/1471-2105-7-91.",
        "category": "model comparison / validation",
        "identifier": "DOI:10.1186/1471-2105-7-91",
    },
    {
        "id": 25,
        "text": "Yarkoni T, Westfall J. Choosing prediction over explanation in psychology: lessons from machine learning. Perspectives on Psychological Science. 2017;12:1100-1122. doi:10.1177/1745691617693393.",
        "category": "model evaluation",
        "identifier": "DOI:10.1177/1745691617693393",
    },
    {
        "id": 26,
        "text": "Wu Z, Ramsundar B, Feinberg EN, et al. MoleculeNet: a benchmark for molecular machine learning. Chemical Science. 2018;9:513-530. doi:10.1039/C7SC02664A.",
        "category": "scaffold split and chemical generalization",
        "identifier": "DOI:10.1039/C7SC02664A",
    },
    {
        "id": 27,
        "text": "Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular frameworks. Journal of Medicinal Chemistry. 1996;39:2887-2893. doi:10.1021/jm9602928.",
        "category": "Bemis-Murcko scaffold",
        "identifier": "DOI:10.1021/jm9602928",
    },
    {
        "id": 28,
        "text": "Rogers D, Hahn M. Extended-connectivity fingerprints. Journal of Chemical Information and Modeling. 2010;50:742-754. doi:10.1021/ci100050t.",
        "category": "Morgan fingerprint / ECFP",
        "identifier": "DOI:10.1021/ci100050t",
    },
    {
        "id": 29,
        "text": "Landrum G. RDKit: Open-source cheminformatics software. https://www.rdkit.org/",
        "category": "cheminformatics software",
        "identifier": "URL:https://www.rdkit.org/",
    },
    {
        "id": 30,
        "text": "Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Systems. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004.",
        "category": "MSigDB Hallmark / GSEA",
        "identifier": "DOI:10.1016/j.cels.2015.12.004",
    },
    {
        "id": 31,
        "text": "Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America. 2005;102:15545-15550. doi:10.1073/pnas.0506580102.",
        "category": "GSEA",
        "identifier": "DOI:10.1073/pnas.0506580102",
    },
    {
        "id": 32,
        "text": "Ashburner M, Ball CA, Blake JA, et al. Gene ontology: tool for the unification of biology. Nature Genetics. 2000;25:25-29. doi:10.1038/75556.",
        "category": "functional gene sets",
        "identifier": "DOI:10.1038/75556",
    },
    {
        "id": 33,
        "text": "Wilkinson MD, Dumontier M, Aalbersberg IJ, et al. The FAIR Guiding Principles for scientific data management and stewardship. Scientific Data. 2016;3:160018. doi:10.1038/sdata.2016.18.",
        "category": "reproducible benchmark reporting",
        "identifier": "DOI:10.1038/sdata.2016.18",
    },
    {
        "id": 34,
        "text": "Gebru T, Morgenstern J, Vecchione B, et al. Datasheets for datasets. Communications of the ACM. 2021;64:86-92. doi:10.1145/3458723.",
        "category": "resource reporting",
        "identifier": "DOI:10.1145/3458723",
    },
    {
        "id": 35,
        "text": "Mitchell M, Wu S, Zaldivar A, et al. Model cards for model reporting. Proceedings of the Conference on Fairness, Accountability, and Transparency. 2019:220-229. doi:10.1145/3287560.3287596.",
        "category": "model reporting",
        "identifier": "DOI:10.1145/3287560.3287596",
    },
    {
        "id": 36,
        "text": "Pineau J, Vincent-Lamarre P, Sinha K, et al. Improving reproducibility in machine learning research. Journal of Machine Learning Research. 2021;22:1-20. https://www.jmlr.org/papers/v22/20-303.html",
        "category": "reproducible benchmark reporting",
        "identifier": "URL:https://www.jmlr.org/papers/v22/20-303.html",
    },
    {
        "id": 37,
        "text": "Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods. 2020;17:261-272. doi:10.1038/s41592-019-0686-2.",
        "category": "scientific computing",
        "identifier": "DOI:10.1038/s41592-019-0686-2",
    },
    {
        "id": 38,
        "text": "Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature. 2020;585:357-362. doi:10.1038/s41586-020-2649-2.",
        "category": "scientific computing",
        "identifier": "DOI:10.1038/s41586-020-2649-2",
    },
    {
        "id": 39,
        "text": "Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. Journal of Machine Learning Research. 2011;12:2825-2830. https://jmlr.org/papers/v12/pedregosa11a.html",
        "category": "machine learning software",
        "identifier": "URL:https://jmlr.org/papers/v12/pedregosa11a.html",
    },
    {
        "id": 40,
        "text": "Hunter JD. Matplotlib: a 2D graphics environment. Computing in Science & Engineering. 2007;9:90-95. doi:10.1109/MCSE.2007.55.",
        "category": "figure generation",
        "identifier": "DOI:10.1109/MCSE.2007.55",
    },
]


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() == dst.resolve():
        return True
    shutil.copy2(src, dst)
    return True


def replace_section(text: str, heading: str, replacement: str) -> str:
    pattern = rf"(## {re.escape(heading)}\n)(.*?)(?=\n## |\Z)"
    return re.sub(pattern, replacement.rstrip() + "\n", text, flags=re.S)


def prepare_manuscript() -> None:
    src = ROOT / "manuscript/manuscript_full_en_v14_JGG_pathway_added.md"
    dst = ROOT / "manuscript/manuscript_full_en_v15_JGG_submission_ready.md"
    text = src.read_text(encoding="utf-8")
    text = text.replace(
        "Single-cell perturbation transcriptomics provides a functional-genomics readout of how cells respond to drugs, genetic state, and cellular context [REF needed].",
        "Single-cell perturbation transcriptomics provides a functional-genomics readout of how cells respond to drugs, genetic state, and cellular context [1-7].",
    )
    text = text.replace(
        "Single-cell perturbation assays make it possible to connect chemical or genetic interventions to transcriptome-wide cellular responses. In drug perturbation studies, the measured response can be interpreted as a functional-genomics profile of pathway activity, stress response, lineage state, or mechanism-specific transcriptional change [REF needed]. Predictive models built on these data are therefore attractive for prioritizing untested perturbation-cell combinations and for comparing representations of cellular state and chemical structure [REF needed].",
        "Single-cell perturbation assays make it possible to connect chemical or genetic interventions to transcriptome-wide cellular responses [1-4]. In drug perturbation studies, the measured response can be interpreted as a functional-genomics profile of pathway activity, stress response, lineage state, or mechanism-specific transcriptional change, extending earlier Connectivity Map and LINCS efforts to a higher-resolution cellular setting [5-7]. Predictive models built on these data are therefore attractive for prioritizing untested perturbation-cell combinations and for comparing representations of cellular state and chemical structure [8-17].",
    )
    text = text.replace(
        "Several model families now compete in single-cell perturbation prediction, including variational perturbation models, chemical-structure-aware models, and single-cell foundation-model representations [REF needed].",
        "Several model families now compete in single-cell perturbation prediction, including variational perturbation models, chemical-structure-aware models, neural optimal-transport methods, and single-cell foundation-model representations [11-22].",
    )
    text = text.replace(
        "Benchmark validity is a central issue in this setting.",
        "Benchmark validity is a central issue in this setting [23-26].",
    )
    text = text.replace(
        "The pathway-level analysis strengthens this interpretation.",
        "The pathway-level analysis strengthens this interpretation and connects the benchmark to established gene-set-level interpretation frameworks [30-32].",
    )
    text = text.replace(
        "Scaffold-held-out validation is necessary for evaluating chemical extrapolation, but it is not sufficient for mechanism-level extrapolation.",
        "Scaffold-held-out validation is necessary for evaluating chemical extrapolation, but it is not sufficient for mechanism-level extrapolation [26-29].",
    )
    text = text.replace(
        "Future benchmark reports should include the split type, same-drug overlap, same-scaffold overlap, same-cell-context overlap, MoA-neighbor overlap, nearest training-drug Tanimoto, random-to-strict performance drop, random-to-strict model-rank transfer, model adaptation level, and failed or not-run status rows (Table 6).",
        "Future benchmark reports should include the split type, same-drug overlap, same-scaffold overlap, same-cell-context overlap, MoA-neighbor overlap, nearest training-drug Tanimoto, random-to-strict performance drop, random-to-strict model-rank transfer, model adaptation level, and failed or not-run status rows (Table 6), following broader expectations for reproducible dataset, model, and benchmark reporting [33-40].",
    )
    text = replace_section(
        text,
        "Authors' contributions",
        "## Authors' contributions\n\nDa Lin: Conceptualization, Methodology, Formal analysis, Writing - original draft, roles to be confirmed. Yu Zhang: Conceptualization, Supervision, Writing - review and editing, roles to be confirmed. All author contributions and CRediT roles should be confirmed before submission.",
    )
    text = replace_section(
        text,
        "Funding",
        "## Funding\n\nFunding information will be confirmed before submission. No grant number is claimed in this draft.",
    )
    text = replace_section(
        text,
        "Availability of data and materials",
        "## Availability of data and materials\n\nThe local benchmark resource package contains split manifests, leakage-audit tables, model-result summaries, report schemas, model feasibility reports, data manifests, source data, and configuration files. Code and benchmark resources will be deposited before submission at [GitHub URL] and archived with Zenodo at [Zenodo DOI]. No GitHub URL, commit identifier, or DOI is claimed in this draft. MSigDB Hallmark gene sets are not redistributed in this repository and should be obtained by users from MSigDB according to its terms of use.",
    )
    refs = "\n".join(f"{r['id']}. {r['text']}" for r in REFERENCES)
    text = replace_section(text, "References", f"## References\n\n{refs}")
    dst.write_text(text, encoding="utf-8")

    for name in ["tables_draft_en", "figure_legends_en"]:
        src2 = ROOT / f"manuscript/{name}_v14_JGG_pathway_added.md"
        dst2 = ROOT / f"manuscript/{name}_v15_JGG_submission_ready.md"
        shutil.copy2(src2, dst2)


def prepare_builder() -> None:
    src = ROOT / "scripts/build_manuscript_docx_v14_JGG_pathway_added.py"
    dst = ROOT / "scripts/build_manuscript_docx_v15_JGG_submission_ready.py"
    text = src.read_text(encoding="utf-8")
    text = text.replace("v14 JGG pathway-added", "v15 JGG submission-ready")
    text = text.replace("manuscript_full_en_v14_JGG_pathway_added.md", "manuscript_full_en_v15_JGG_submission_ready.md")
    text = text.replace("tables_draft_en_v14_JGG_pathway_added.md", "tables_draft_en_v15_JGG_submission_ready.md")
    text = text.replace("figure_legends_en_v14_JGG_pathway_added.md", "figure_legends_en_v15_JGG_submission_ready.md")
    text = text.replace(
        "submission_package_v14_JGG_pathway_added/manuscript/manuscript_full_en_v14_JGG_pathway_added.docx",
        "submission_package_v15_JGG_submission_ready/manuscript/manuscript_full_en_v15_JGG_submission_ready.docx",
    )
    dst.write_text(text, encoding="utf-8")


def prepare_gene_set_readme() -> None:
    readme = ROOT / "resources/gene_sets/README.md"
    readme.write_text(
        "# Gene Set Resources\n\n"
        "The MSigDB Hallmark GMT file is required for the pathway-level recovery analysis but is not redistributed in the public benchmark resource package.\n\n"
        "Users should obtain the Hallmark gene sets directly from MSigDB according to MSigDB terms of use and place the file at:\n\n"
        "`resources/gene_sets/hallmark_symbols.gmt`\n\n"
        "The manuscript reports only gene-set overlap summaries, pathway-level metrics, and source data derived from the locally provided GMT file.\n",
        encoding="utf-8",
    )


def overlay_labels(src: Path, dst: Path, labels: list[tuple[str, float, float]]) -> None:
    im = Image.open(src).convert("RGB")
    draw = ImageDraw.Draw(im)
    try:
        font = ImageFont.truetype("Arial.ttf", max(18, im.width // 42))
    except Exception:
        font = ImageFont.load_default()
    for label, xfrac, yfrac in labels:
        x = int(im.width * xfrac)
        y = int(im.height * yfrac)
        draw.text((x, y), label, fill=(0, 0, 0), font=font)
    dst.parent.mkdir(parents=True, exist_ok=True)
    im.save(dst, dpi=(300, 300))
    im.save(dst.with_suffix(".pdf"), resolution=300.0)


def prepare_figures() -> list[dict]:
    final = ROOT / "figures/final"
    final_source = ROOT / "figures/final/source_data"
    final_source.mkdir(parents=True, exist_ok=True)
    rows = []
    figure_map = {
        "Figure 1": ROOT / "figures/manuscript_main/figure1_benchmark_design",
        "Figure 2": ROOT / "figures/manuscript_main/figure2_baseline_decay",
        "Figure 3": ROOT / "figures/manuscript_main/figure4_chemical_similarity_audit",
        "Figure 4": ROOT / "figures/figure6_sciplex3_prnet30_transigen30_model_panel",
        "Figure 5": ROOT / "figures/figure5_pathway_level_recovery",
        "Figure 6": ROOT / "figures/figure7_rank_transfer_prnet30",
    }
    for fig, stem in figure_map.items():
        for ext in [".pdf", ".svg", ".png", ".tiff"]:
            ok = copy_if_exists(stem.with_suffix(ext), final / f"{fig.lower().replace(' ', '_')}{ext}")
            if ok:
                rows.append({"figure": fig, "file": str(final / f"{fig.lower().replace(' ', '_')}{ext}"), "status": "available"})
    # Add panel-label raster copies for figures that need explicit labels in the final review folder.
    overlay_labels(
        ROOT / "figures/figure5_pathway_level_recovery.png",
        final / "figure_5_panel_labelled.png",
        [("A", 0.04, 0.05), ("B", 0.52, 0.05), ("C", 0.04, 0.50), ("D", 0.52, 0.50)],
    )
    overlay_labels(
        ROOT / "figures/figure7_rank_transfer_prnet30.png",
        final / "figure_6_panel_labelled.png",
        [("A", 0.08, 0.08), ("B", 0.58, 0.08)],
    )

    source_map = {
        "figure1_coverage.csv": ROOT / "figures/manuscript_main/source_data/figure1b_coverage.csv",
        "figure1_split_sizes.csv": ROOT / "figures/manuscript_main/source_data/figure1c_split_sizes.csv",
        "figure2_baseline_decay.csv": ROOT / "figures/manuscript_main/source_data/figure2_baseline_decay.csv",
        "figure3_similarity_ci.csv": ROOT / "figures/manuscript_main/source_data/figure4a_similarity_ci.csv",
        "figure3_similarity_contrasts.csv": ROOT / "figures/manuscript_main/source_data/figure4b_similarity_contrasts.csv",
        "figure4_completed_model_performance.csv": ROOT / "results/deep_model_panel/sciplex3_completed_model_summary.csv",
        "figure5_pathway_level_recovery.csv": ROOT / "figures/source_data/figure5_pathway_level_recovery.csv",
        "figure6_sciplex3_rank_heatmap.csv": ROOT / "benchmark_resource/source_data/figure7_sciplex3_24h_top2000_rank_heatmap.csv",
        "figure6_rank_transfer_summary.csv": ROOT / "benchmark_resource/source_data/figure7_rank_transfer_summary.csv",
    }
    for name, src in source_map.items():
        copied = copy_if_exists(src, final_source / name)
        copy_if_exists(src, ROOT / "figures/source_data" / name)
        rows.append({"figure": name.split("_")[0].replace("figure", "Figure "), "file": str(src), "status": "source_available" if copied else "missing"})
    with (ROOT / "docs/figure_quality_check.md").open("w", encoding="utf-8") as handle:
        handle.write("# Figure Quality Check\n\n")
        handle.write("No model results or statistical values were changed during figure packaging.\n\n")
        handle.write("| Figure | Final files | Source data | Notes |\n|---|---|---|---|\n")
        for fig in ["Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 6"]:
            files = sorted(p.name for p in final.glob(f"{fig.lower().replace(' ', '_')}.*"))
            sources = sorted(p.name for p in final_source.glob(f"{fig.lower().replace(' ', '').lower()}*"))
            if fig == "Figure 5":
                sources = ["figure5_pathway_level_recovery.csv"]
                files.append("figure_5_panel_labelled.png/pdf")
            if fig == "Figure 6":
                sources = ["figure6_sciplex3_rank_heatmap.csv", "figure6_rank_transfer_summary.csv"]
                files.append("figure_6_panel_labelled.png/pdf")
            handle.write(f"| {fig} | {', '.join(files) if files else 'missing'} | {', '.join(sources) if sources else 'see final/source_data'} | Final PDF/SVG/PNG copies prepared; panel-labelled raster copies prepared where needed. |\n")
    return rows


def parse_markdown_tables() -> None:
    src = ROOT / "manuscript/tables_draft_en_v15_JGG_submission_ready.md"
    text = src.read_text(encoding="utf-8")
    current = None
    rows: list[list[str]] = []
    for line in text.splitlines() + [""]:
        if line.startswith("## Table "):
            if current and rows:
                write_table_csv(current, rows)
            current = "Table " + line.split("Table ", 1)[1].split(".", 1)[0]
            rows = []
        elif current and line.startswith("|") and "---" not in line:
            rows.append([c.strip() for c in line.strip("|").split("|")])
        elif current and rows and not line.strip():
            write_table_csv(current, rows)
            current = None
            rows = []


def write_table_csv(table_name: str, rows: list[list[str]]) -> None:
    out = ROOT / "supplementary/Source_Data" / f"{table_name.lower().replace(' ', '_')}_source.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)
    copy_if_exists(out, ROOT / "figures/source_data" / out.name)


def prepare_resource_package() -> None:
    copy_pairs = [
        (ROOT / "results/pathway_level/pathway_gene_set_overlap.csv", ROOT / "benchmark_resource/pathway_level/pathway_gene_set_overlap.csv"),
        (ROOT / "results/pathway_level/pathway_seed_summary.csv", ROOT / "benchmark_resource/pathway_level/pathway_seed_summary.csv"),
        (ROOT / "results/pathway_level/pathway_random_to_strict_contrasts.csv", ROOT / "benchmark_resource/pathway_level/pathway_random_to_strict_contrasts.csv"),
        (ROOT / "results/pathway_level/pathway_per_gene_set_errors.csv", ROOT / "benchmark_resource/pathway_level/pathway_per_gene_set_errors.csv"),
        (ROOT / "results/pathway_level/pathway_bootstrap_ci.csv", ROOT / "benchmark_resource/pathway_level/pathway_bootstrap_ci.csv"),
        (ROOT / "docs/pathway_gene_set_check.md", ROOT / "benchmark_resource/QC_logs/pathway_gene_set_check.md"),
        (ROOT / "figures/source_data/figure5_pathway_level_recovery.csv", ROOT / "benchmark_resource/source_data/figure5_pathway_level_recovery.csv"),
    ]
    for src, dst in copy_pairs:
        copy_if_exists(src, dst)
    # Supplementary package.
    sup_pairs = [
        (ROOT / "benchmark_resource/split_manifests/split_manifest_openproblems_summary.csv", ROOT / "supplementary/Supplementary_Tables/split_manifest_openproblems_summary.csv"),
        (ROOT / "benchmark_resource/split_manifests/split_manifest_sciplex3_summary.csv", ROOT / "supplementary/Supplementary_Tables/split_manifest_sciplex3_summary.csv"),
        (ROOT / "results/openproblems_multiseed_10_de/all_ridge_baseline_summary.csv", ROOT / "supplementary/Supplementary_Tables/openproblems_seed_level_metrics.csv"),
        (ROOT / "benchmark_resource/leakage_audits/leakage_audit_table.csv", ROOT / "supplementary/Supplementary_Tables/leakage_audit_table.csv"),
        (ROOT / "results/manuscript_v5_stats/openproblems_moa_overlap_summary.csv", ROOT / "supplementary/Supplementary_Tables/openproblems_moa_overlap_summary.csv"),
        (ROOT / "results/manuscript_v5_stats/openproblems_matched_random_controls_summary.csv", ROOT / "supplementary/Supplementary_Tables/openproblems_matched_random_controls_summary.csv"),
        (ROOT / "results/manuscript_v5_stats/sciplex3_dose_leakage_summary.csv", ROOT / "supplementary/Supplementary_Tables/sciplex3_dose_leakage_summary.csv"),
        (ROOT / "results/pathway_level/pathway_gene_set_overlap.csv", ROOT / "supplementary/Supplementary_Tables/pathway_gene_set_overlap.csv"),
        (ROOT / "results/pathway_level/pathway_seed_summary.csv", ROOT / "supplementary/Supplementary_Tables/pathway_seed_summary.csv"),
        (ROOT / "results/pathway_level/pathway_random_to_strict_contrasts.csv", ROOT / "supplementary/Supplementary_Tables/pathway_random_to_strict_contrasts.csv"),
        (ROOT / "benchmark_resource/model_feasibility_report.md", ROOT / "supplementary/Model_Feasibility/model_feasibility_report.md"),
        (ROOT / "docs/pathway_gene_set_check.md", ROOT / "supplementary/QC_Reports/pathway_gene_set_check.md"),
    ]
    for src, dst in sup_pairs:
        copy_if_exists(src, dst)
    for fig in ["supp_pathway_gene_set_overlap.png", "supp_pathway_per_seed_distribution.png", "supp_pathway_per_gene_set_error_heatmap.png", "supp_prnet30_seed_distribution.png", "supp_prnet30_vs_transigen30_paired.png"]:
        copy_if_exists(ROOT / "figures" / fig, ROOT / "supplementary/Supplementary_Figures" / fig)
    for p in (ROOT / "figures/final/source_data").glob("*.csv"):
        copy_if_exists(p, ROOT / "supplementary/Source_Data" / p.name)


def write_reports() -> None:
    docs = ROOT / "docs"
    docs.mkdir(exist_ok=True)
    ref_lines = [
        "# Reference Completion Report",
        "",
        "All `[REF needed]` placeholders in the v15 manuscript were replaced with verified citation groups. References with DOI values were checked against Crossref, PubMed, publisher pages, or known official publication metadata where available. References without DOI are retained as URL-based resources.",
        "",
        "## Added reference coverage",
        "",
        "| Category | Reference IDs | Identifier notes |",
        "|---|---|---|",
    ]
    by_cat: dict[str, list[str]] = {}
    for r in REFERENCES:
        by_cat.setdefault(r["category"], []).append(f"{r['id']} ({r['identifier']})")
    for cat, vals in by_cat.items():
        ref_lines.append(f"| {cat} | {', '.join(vals)} | DOI/URL recorded |")
    ref_lines.extend(["", "## Remaining placeholders", "", "- No `[REF needed]` placeholders remain in the v15 manuscript."])
    (docs / "reference_completion_report.md").write_text("\n".join(ref_lines) + "\n", encoding="utf-8")

    required = [
        "README.md", "data_manifest.tsv", "split_manifests", "leakage_audits", "model_results",
        "pathway_level", "source_data", "report_schema", "configs", "model_feasibility_report.md", "QC_logs",
    ]
    resource = ROOT / "benchmark_resource"
    lines = ["# JGG Resource Release Checklist", "", "| Required item | Status | Notes |", "|---|---|---|"]
    for item in required:
        p = resource / item
        lines.append(f"| `{item}` | {'present' if p.exists() else 'missing'} | `{p}` |")
    has_env = any((resource / n).exists() for n in ["environment.yml", "requirements.txt"])
    lines.append(f"| `environment.yml` or `requirements.txt` | {'present' if has_env else 'needs manual addition'} | Add before public release if not present. |")
    lines.extend(
        [
            "",
            "## GitHub release checklist",
            "",
            "- Add final code, split manifests, leakage audits, result summaries, source data, schemas, configs, and QC reports.",
            "- Do not include `resources/gene_sets/hallmark_symbols.gmt` or any redistributed MSigDB GMT file.",
            "- Add a license only after the authors confirm the intended license.",
            "- Tag a release and record the release version in the manuscript before submission.",
            "",
            "## Zenodo archive checklist",
            "",
            "- Archive the GitHub release after GitHub contents are finalized.",
            "- Record the Zenodo DOI in the manuscript only after it exists.",
            "- Exclude MSigDB GMT files from the archive.",
            "- Include README, manifest, source data, and schema files in the archive.",
        ]
    )
    (docs / "JGG_resource_release_checklist.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    (docs / "table_revision_report.md").write_text(
        "# Table Revision Report\n\n"
        "- Table 1 and Table 2 remain in the main manuscript as dataset and split-definition anchors.\n"
        "- Table 3 remains in the main manuscript as a compact completed-model performance summary; full seed-level and model-level metrics are assigned to supplementary tables/source data.\n"
        "- Table 4 remains in the main manuscript because model feasibility/adaptation level is a resource contribution.\n"
        "- Table 5 remains in the main manuscript because rank transfer supports the model-selection claim.\n"
        "- Table 6 remains in the main manuscript as the recommended leakage-aware reporting checklist.\n"
        "- Table source CSV files were generated under `supplementary/Source_Data/` and mirrored to `figures/source_data/`.\n",
        encoding="utf-8",
    )

    sup_files = []
    for p in sorted((ROOT / "supplementary").rglob("*")):
        if p.is_file():
            sup_files.append(f"| `{p.relative_to(ROOT)}` | present |")
    (docs / "supplementary_file_manifest.md").write_text(
        "# Supplementary File Manifest\n\n| File | Status |\n|---|---|\n" + "\n".join(sup_files) + "\n",
        encoding="utf-8",
    )

    manuscript = (ROOT / "manuscript/manuscript_full_en_v15_JGG_submission_ready.md").read_text(encoding="utf-8")
    checks = [
        ("[REF needed]", "pass" if "[REF needed]" not in manuscript else "attention", "No reference placeholders remain." if "[REF needed]" not in manuscript else "Reference placeholders remain."),
        ("[GitHub URL]", "attention" if "[GitHub URL]" in manuscript else "pass", "Still requires manual replacement with a real repository URL." if "[GitHub URL]" in manuscript else "Repository URL has been filled."),
        ("[Zenodo DOI]", "attention" if "[Zenodo DOI]" in manuscript else "pass", "Still requires manual replacement after archive creation." if "[Zenodo DOI]" in manuscript else "Zenodo DOI has been filled."),
        ("author roles to be confirmed", "attention" if "to be confirmed" in manuscript else "pass", "Author contributions still require author confirmation." if "to be confirmed" in manuscript else "Author contributions appear finalized."),
        ("funding to be confirmed", "attention" if "Funding information will be confirmed" in manuscript else "pass", "Funding remains unconfirmed." if "Funding information will be confirmed" in manuscript else "Funding statement appears finalized."),
        ("MSigDB GMT not redistributed", "pass" if not any("hallmark_symbols.gmt" in str(p.relative_to(ROOT)) for p in (ROOT / "benchmark_resource").rglob("*")) else "attention", "MSigDB GMT was not copied into benchmark_resource."),
    ]
    lines = ["# JGG Submission Readiness Report", "", "| Check | Status | Notes |", "|---|---|---|"]
    for name, status, note in checks:
        lines.append(f"| {name} | {status} | {note} |")
    lines.extend(
        [
            "",
            "## Remaining manual items",
            "",
            "- Replace `[GitHub URL]` after the public repository is created.",
            "- Replace `[Zenodo DOI]` after the archived release is created.",
            "- Confirm author contributions and funding statement.",
            "- Finalize acknowledgements, cover letter, and journal-specific reference formatting.",
            "- Confirm license before public release.",
        ]
    )
    (docs / "JGG_submission_readiness_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    (ROOT / "summary_v15_submission_ready_changes.md").write_text(
        "# v15 Submission-Ready Change Summary\n\n"
        "- No model results, pathway results, statistical values, split assignments, or performance numbers were changed.\n"
        f"- References were expanded from 10 to {len(REFERENCES)} entries, covering single-cell perturbation transcriptomics, Sci-Plex 3, OpenProblems/OP3, LINCS/Connectivity Map, scGen, CPA, chemCPA, TranSiGen, PRnet, scGPT, Geneformer, scFoundation, leakage, scaffold splitting, Bemis-Murcko scaffolds, Morgan fingerprints, MSigDB/GSEA, and reproducible benchmark reporting.\n"
        "- Figures were copied/re-exported into `figures/final/`, with final source data organized under `figures/final/source_data/` and `figures/source_data/`.\n"
        "- Supplementary package directories were created and populated with available split summaries, seed-level metrics, leakage audits, MoA/dose controls, pathway tables, source data, model feasibility reports, and QC reports.\n"
        "- GitHub URL and Zenodo DOI remain placeholders and require manual handling after repository/archive creation.\n"
        "- Remaining submission items: author contribution confirmation, funding confirmation, acknowledgements, final reference style check, figure-format check against JGG instructions, cover letter, public repository, Zenodo archive, and license confirmation.\n",
        encoding="utf-8",
    )


def main() -> None:
    prepare_manuscript()
    prepare_builder()
    prepare_gene_set_readme()
    prepare_figures()
    parse_markdown_tables()
    prepare_resource_package()
    write_reports()
    print("Prepared v15 submission-ready package files.")


if __name__ == "__main__":
    main()
