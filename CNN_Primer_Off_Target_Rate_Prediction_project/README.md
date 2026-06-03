# CNN — Primer Off-Target Rate Prediction

## Introduction

This project builds and trains a Convolutional Neural Network (CNN) for the prediction of primer off-target rate, given a DNA target and a primer. The model is tracked with MLflow and trained on alignments encoded as images, using the number of UMI counts at each site, from which the off-target rate is derived.

Throughout this document, `$PROJECT_ROOT` refers to the root directory of the repository (the folder containing `pyproject.toml` and `README.md`). Replace it with the actual path on your machine, or export it once:

```bash
export PROJECT_ROOT=~/CNN_Primer_Off_Target_Rate_Prediction_project
```

---

## Getting started

> **Important:** the `pysam` Python library is not supported on Windows. On Windows, the project must be run under **WSL (Windows Subsystem for Linux)**.

### Windows + WSL notes

If you are working on Windows using WSL, a few considerations ensure the project runs correctly.

**Repository location.** Place the repository inside the WSL filesystem:

- Copy the HTTPS or SSH URL from GitHub and `git clone` it into `~/` (i.e. `/home/<your-user>/`).
- Open the repository in the WSL environment: in Visual Studio Code, type `>` in the command bar and select **WSL: Connect to WSL** (the WSL extension must be installed).

**Line endings (CRLF vs LF).** Operating systems use different line endings — Windows uses CRLF, Linux (WSL) uses LF. **Scripts must use LF endings to execute correctly in WSL.** After scripts are executed and outputs are generated, files may revert to CRLF when pushed from a Windows environment; Git handles this automatically and prints a warning on the terminal.

---

## 1. Installation

This repository uses the **`uv`** environment manager.

1. Install `uv` (if not already installed):
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh   # Linux / macOS
   ```

2. **First time only** — create the environment. Skip this step if the repository already contains a `.venv` and a `uv.lock`, and go straight to step 3.
   ```bash
   uv init        # create pyproject.toml (skip if it already exists)
   uv venv        # create the virtual environment
   uv add -r requirements.txt   # install from requirements.txt (skip if pyproject.toml and uv.lock exist)
   ```

3. Install dependencies from the locked specification:
   ```bash
   uv sync        # uses pyproject.toml (declared deps) + uv.lock (pinned versions)
   ```

4. Activate the environment:
   ```bash
   source .venv/bin/activate
   ```

5. Add further libraries as needed:
   ```bash
   uv add <library>     # e.g. uv add scikit-learn
   ```

---

## 2. Data-processing pipeline (alignments → images)

Ensure the `reports/` and `report_test/` folders exist (create them if not).

### `01f_align_collap_encode_and_metada_1.py`

```bash
cd $PROJECT_ROOT/scripts
chmod +x 01f_align_collap_encode_and_metada_1.py   # skip if already executable
./01f_align_collap_encode_and_metada_1.py
```

For each raw `{name}_alignments.csv` file, this generates a metadata file `{name}_metadata.json` and an encoded array `{name}_images_and_labels_uint8_.npz`.

### `01f_verify_npz_imglab.py`

```bash
chmod +x 01f_verify_npz_imglab.py   # skip if already executable
./01f_verify_npz_imglab.py
```

For each `{name}_images_and_labels_uint8_.npz`, this writes a per-file folder inside `reports/`, and a data-quality report (in both Markdown and JSON) summarising the encoded data.

The terminal output of this step is recorded in `logs/terminal_print_01f.txt`.

---

## 3. Training a model

> **Note:** the data used to develop this project (and to train the CNN) is confidential. The repository is therefore prepared to run a mock model trained on a few images in a single `.npz` file; the models developed during the project are evaluated on `test_set_2_images_and_labels_uint8_.npz`.

Move the test set into the held-out test folder and remove it from the processed directory:

```bash
cd $PROJECT_ROOT/dada/processed/images_labels_CollapsedAlign01f_npz
cp test_set_2_images_and_labels_uint8_.npz $PROJECT_ROOT/dada/hold_out_test_data/
rm test_set_2_images_and_labels_uint8_.npz
cd $PROJECT_ROOT
```

Then run one of the three experiment scripts from the project root:

```bash
chmod +x one_run_bypanel.py one_run_byprimer.py one_run_DHS001z.py   # skip if already executable

./one_run_bypanel.py     # panel-disjoint split (generalisation to unseen panels)
./one_run_byprimer.py    # primer-disjoint split (generalisation to unseen primers)
./one_run_DHS001z.py     # single-panel pilot (DHS-001Z)
```

Each run trains the model with MLflow tracking and evaluates it once on the held-out test data. Outputs are written to `reports/` (training), `report_test/` (held-out test) and `outputs/` (figures), and logged to the MLflow store (`mlflow.db`).

---

## Repository layout

```
$PROJECT_ROOT/
├── dada/                              data
│   ├── raw/                           {name}_alignments.csv source files
│   ├── processed/
│   │   └── images_labels_CollapsedAlign01f_npz/   encoded .npz inputs
│   └── hold_out_test_data/            held-out test .npz
├── scripts/                           preprocessing (encode + verify)
├── src/project_package/               core library
│   ├── config.py  data.py  model.py  train.py  evaluate.py  test.py
│   ├── cross_validation.py  utils.py
│   └── tests/                         unit / smoke tests
├── one_run_bypanel.py                 panel-disjoint experiment
├── one_run_byprimer.py                primer-disjoint experiment
├── one_run_DHS001z.py                 single-panel pilot
├── models/                            saved trained models (.pt)
├── reports/  report_test/  outputs/   run artifacts and figures
├── logs/                              terminal logs
├── mlflow.db                          MLflow tracking store
├── pyproject.toml  uv.lock            environment specification
├── requirements.txt                   dependency list (used at first setup)
└── README.md
```
