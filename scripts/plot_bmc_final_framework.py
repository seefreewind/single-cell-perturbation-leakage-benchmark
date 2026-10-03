"""Reproduce the text-only final Figure 1 revision without model fitting."""
from pathlib import Path
import importlib.util
import inspect

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / 'scripts/plot_bmc_revision_figures.py'
spec = importlib.util.spec_from_file_location('existing_figure_layout', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
assets = ROOT / 'BMC_revision_20261003/final_polish_assets'
assets.mkdir(parents=True, exist_ok=True)
module.FIG = module.SRC = assets
source = inspect.getsource(module.fig1)
changes = {
    'Leakage-aware benchmark evaluation defines what each score can claim':
        'Evaluation workflow for structured perturbation-response data',
    'Generalization is interpreted as a model--data--split--claim system':
        'Different splits evaluate different prediction settings',
    'D. Portable evidence': 'D. Reproducibility outputs',
    'familiar drug/cell neighborhoods': 'familiar chemical neighborhoods',
    'new chemical framework': 'unseen scaffold',
    'new cell + new scaffold intersection': 'unseen cell + unseen scaffold',
    'Output: interpretable scores, not a single undifferentiated claim of model generalization':
        'Outputs: audit fields, predictions, metrics, and source data',
    'fontsize=6.8, va="center"': 'fontsize=6.4, va="center"',
}
for old, new in changes.items():
    assert old in source, old
    source = source.replace(old, new)
exec(compile(source, str(path), 'exec'), module.__dict__)
module.fig1()
