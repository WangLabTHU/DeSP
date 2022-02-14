# DeSP: a systematic error Simulation Platform for DNA data storage channel

## Getting started with the web app

### Live Demo
For a live demo, http://170.106.110.86/ (Silicon Valley, US)

![app](https://github.com/WangLabTHU/DNA-D2S/blob/main/files/webapp.jpg)

In case that the web app doesn't response, you can also:

### Run the web app locally 

**Clone or download this repository, and install the dependencies:**

```bash
git clone https://github.com/WangLabTHU/DeSP

# install packages
# under python 3.x
pip install matplotlib numpy plotly seaborn streamlit scipy reedsolo prettytable

```

**Navigate to the project folder,  and run the app:**

```bash

# move to the project folder
cd download/DeSP main

streamlit run main.py
```

Navigate to [localhost:8501](https://localhost:8501/). You should see the app running in your broswer :)


## Embed the channel model in your research pipeline
You can construct a model, use the model to generate simulated sequencing readouts, and analysis the sequencing reults following the scripts
below:
```python
from Model.Model import *
import Model.config as config

# editting parameters of the channel
arg = config.DEFAULT_PASSER
arg.seq_depth = 10

# construct a channel by linking modules
Modules = [
  ('synthesizing',Synthesizer(arg)),
  ('decaying',Decayer(arg)),
  ('pcring',PCRer(arg = arg)),
  ('sampling',Sampler(arg = arg)),
  ('sequencing',Sequencer(arg))
]
Model = DNA_Channel_Model(Modules)

# load the data, and use the model to generate simulated sequencing results
with open('files/lena.dna') as f:
dnas = f.readlines()
in_dnas = [dna.split('\n')[0] for dna in dnas]
out_dnas = Model(in_dnas)

# examine the output dnas
from Analysis.Analysis import inspect_distribution, examine_strand
inspect_distribution(out_dnas, show = True) # oligo number and error number distribution of the entire sequencing results
examine_strand(out_dnas, index = index) # sequencing readouts and voting result of a single sequence.
```

## Jupyter notebooks and experiment data
### Notebooks
Notebooks to reproduce figures in the paper are also provided.

* **Analysis of Individual Stages of DNA data channel.ipynb** (Corresponding to part 2.1 of the paper)
* **Validation.ipynb** (Correspond to part 3.1 of the paper)
* **Step-by-step simulation.ipynb** (Part 3.2 of the paper)
* **Encoding Design.ipynb** (Part 3.3 of the paper)

To use the notebooks, please download the project and install the dependencie following the instructions above first.

### Data
The experimental sequencing data used in the paper are provided in seq_data.csv and seq_data.pkl. Data in the two files are the same,
the pickle file is just for fast loading in python scripts. Please download the data from https://cloud.tsinghua.edu.cn/d/157100ca4c8e4fd387e1/, and move the data to files/ before running the validation notebook.


## Acknowledgements
Code about the DNA Fountain code is adopted from https://github.com/TeamErlich/dna-fountain.



