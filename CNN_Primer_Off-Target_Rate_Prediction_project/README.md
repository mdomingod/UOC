# Introduction 
This project focuses on building and training a Convolutional Neural Network (CNN) for the prediction of primer off-target rate given a DNA target and a primer.
The model is track by MLflow and trained by the alignments encoded into an images and the number of umi counts which are used to get the off-target rate.

# Getting Started
Note: pysam python library is not supported by windows so for Windows operation systems is required the use of wsl (Windows Subsystem for Linux).

## 1. Installation Process

### ⚠️ Windows + WSL Notes
If you are working on Windows using WSL (Windows Subsystem for Linux), there are a few important considerations to ensure the project runs correctly.

#### 📁 Repository Location
When cloning this repository on Windows, it is strongly recommended to place it inside the WSL filesystem:  
- Copy the https direction from GitHub and git clone into "home/<your-user>/".

- Open the repository directory in WSL enviroment, [WSL:Ubuntu-24.04] terminal.
The use of Visual Code Studio is good option.


- The first time to `git pull`, `git add`, `git commit -m` and `git push` on git bash, a message will appear asking the ownership exception to execute git commands to explicitly trust that repository so Git can allow operations on it.

- Line Ending Differences (CRLF vs LF). Operating systems use different line endings:
Windows → CRLF
Linux (WSL) → LF
**Scripts must use LF endings to execute correctly in WS**  

-  After script execution and outputs generation, scripts my revert to CRLF when pushed from Windows environments. 
This will be done automatically by git and a warning message will be displayed on the terminal.

### This repo uses UV enviroment:
- Open the repository folder by an interface as Visual Studio code.
1. install uv () (in case you don't have it already installed)
    - linux/mac: `curl -LsSf https://astral.sh/uv/install.sh | sh`  


2. Create a uv enviroment (ONLY THE FIRST TIME). 
Skip this step if the repo has already an enviroment created and go ahead to the next step.
    - Create pyproject.toml: `uv init`
    - Create the environment: `uv venv` 
    - Install the requirments.txt: `uv add -r requirements.txt` 

3. Install dependencies: `uv sync`
Only if the project has already pyproject.toml
    - pyproject.toml (declares dependencies)
    - uv.lock (freezes exact versions)      

4. Activate venv: `source .venv/bin/activate`   

5. Install more libraries: `uv add <library>` (e.g. `uv add scikit-learn`)  


## Software dependencies and Latest releases

Check "pyproject.toml" file:

dependencies = [
    "biopython==1.86",
    "contourpy==1.3.0",
    "cycler==0.12.1",
    "fonttools==4.60.2",
    "importlib-resources==6.5.2",
    "iprogress>=0.4",
    "kiwisolver==1.4.7",
    "matplotlib>=3.10.8",
    "mlflow>=1.27.0",
    "nbconvert>=7.17.0",
    "notebook>=7.5.5",
    "numpy>=2.4.3",
    "packaging==25.0",
    "pandas>=3.0.1",
    "pillow>=12.1.1",
    "pyparsing==3.3.1",
    "pysam>=0.23.3",
    "python-dateutil==2.9.0.post0",
    "pytz==2025.2",
    "regex>=2026.4.4",
    "seaborn>=0.13.2",
    "six==1.17.0",
    "torch>=2.12.0",
    "torchvision>=0.27.0",
    "tqdm>=4.67.3",
    "tzdata==2025.3",
    "zipp==3.23.0",
]

Make sure you already added this to pyproject.toml:
[tool.setuptools.packages.find]
where = ["src"]
Otherwise uv won’t know where your code is.

Install the project in editable mode with `uv pip install -e .` so Python can import your package directly from the src directory and reflect changes instantly.




## API references

 MLflow API 
`mlflow.set_experiment("name")`           # Group runs together (e.g., all baseline runs)
`mlflow.start_run(run_name="...")`        # Start tracking a run
`mlflow.log_param(key, value)`            # Log one parameter
`mlflow.log_params(dict)`                 # Log many parameters at once
`mlflow.log_metric(key, value, step=N)`   # Log one metric, optionally at a specific step (epoch)
`mlflow.log_metrics(dict, step=N)`        # Log many metrics at once
`mlflow.log_artifact(filepath)`           # Log a file (model, report, plot)
`mlflow.end_run()`                        # Stop tracking (auto-called if using context manager)


Command to check MLFlow ui:  `mlflow ui --backend-store-uri sqlite:///mlflow.db --port 5000`

# Build and Test
Since 

# Contribute
TODO: Explain how other users and developers can contribute to make your code better. 

If you want to learn more about creating good readme files then refer the following [guidelines](https://docs.microsoft.com/en-us/azure/devops/repos/git/create-a-readme?view=azure-devops). You can also seek inspiration from the below readme files:
- [ASP.NET Core](https://github.com/aspnet/Home)
- [Visual Studio Code](https://github.com/Microsoft/vscode)
- [Chakra Core](https://github.com/Microsoft/ChakraCore)