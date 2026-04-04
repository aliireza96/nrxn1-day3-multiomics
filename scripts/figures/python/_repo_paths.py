from pathlib import Path


def find_project_root() -> Path:
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        if (parent / ".nrxn1_project_root").exists():
            return parent
    raise RuntimeError("Could not find project root sentinel .nrxn1_project_root")


PROJECT_ROOT = find_project_root()
OUTPUTS_DIR = PROJECT_ROOT / "outputs"
