# src/project_package/config.py

from dataclasses import dataclass, field, asdict
from pathlib import Path

@dataclass
class Config:
    # ---- paths ----
    project_root: Path = Path(__file__).parent.parent.parent
    data_dir:     Path = field(init=False)
    data_dir_dhs001z: Path = field(init=False)
    npz_files:    list = field(default_factory=list)
    output_dir:   Path = field(init=False)
    data_test_dir: Path = field(init=False)

    # ---- image encoding ----
    img_height:   int = 5
    img_width:    int = 66
    img_channels: int = 2

    # ---- model architecture ----
    n_filters:        int   = 32
    dropout_rate:     float = 0.25
    fc_units:         int   = 128
    allow_negative_output: bool = False

    # ---- training ----
    batch_size:    int   = 128
    epochs:        int   = 15
    learning_rate: float = 1e-3
    momentum:      float = 0.9
    patience:      int   = 10
    lr_patience:   int   = 5
    val_fraction:  float = 0.30
    #test_fraction: float = 0.10
    # Paper-faithful: list of panel IDs to hold out for final test eval.
    seed:          int   = 42

    # ---- DataLoader / hardware ----
    n_workers:    int  = 4
    pin_memory:   bool = True
    device:       str  = "auto"      # 'auto' | 'cpu' | 'cuda' | 'mps'

    # ---- reproducibility (PyTorch reproducibility docs) ----
    deterministic:  bool = True      # cudnn.deterministic + use_deterministic_algorithms
    cudnn_benchmark: bool = False    # must be False for full determinism

    # ---- cross-validation ----
    n_folds:       int = 3

    # ---- evaluation ----
    n_bootstrap_iterations: int = 1000
    bootstrap_confidence:   float = 0.95

    # ---- experiment tracking ----
    experiment_name:  str = "pordle_cnn"
    run_name:         str = "baseline"
    mlflow_tracking_uri: str = "sqlite:///mlflow.db"   # empty = use local ./mlruns folder

    # What to log
    log_model:           bool = True      # log the full model artifact
    log_dataset_info:    bool = True      # log shapes, primer counts, hash
    log_environment:     bool = True      # log pip freeze + python version
    log_git_state:       bool = True      # log current commit SHA
    log_per_epoch_plots: bool = False     # heavy; only for short debug runs
 
    # Free-form tags applied to every run (used for filtering in the UI)
    tags: dict = field(default_factory=lambda: {
        "model_type": "CNN",
        "encoding":   "2D_5bit",
        "framework":  "pytorch",
        "project":    "primer_offtarget",
    })


    def __post_init__(self):
        self.data_test_dir = self.project_root / "data" / "hold_out_test_data"  
        self.data_dir = self.project_root / "data" / "processed" / "images_labels_CollapsedAlign01f_npz"
        self.data_dir_dhs001z = self.project_root / "data" / "need_to_be_combined" / "processed" /"images_labels_CollapsedAlign_npz" 
        self.output_dir = self.project_root / "outputs"
        self.reports_dir = self.project_root / "reports"
        self.report_test_dir = self.project_root / "report_test"
        self.project_root = self.project_root 
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.reports_dir.mkdir(exist_ok=True, parents=True)

    def to_dict(self):
        """Flat dict of all fields. Paths converted to strings for MLflow."""
        return {k: str(v) if isinstance(v, Path) else v for k, v in asdict(self).items()}