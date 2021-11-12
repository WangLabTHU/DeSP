# SOS-DNA: a systematic error Simulation and encoding Optimization System for DNA data storage channel 

## Getting started

### Live Demo
For a live demo, visit http://192.144.226.158/ (server located in Beijing, China) or http://170.106.110.86/ (Silicon Valley, US)

![app](https://github.com/WangLabTHU/SOSDNA/blob/main/webapp.jpg)

In case that the web app doesn't response, you can also:

### Run the web app locally 

**Clone or download this repository, and install the dependencies:**

```bash
git clone https://github.com/WangLabTHU/SOSDNA

# install packages
# under python 3.x
pip install matplotlib numpy plotly seaborn streamlit scipy reedsolo prettytable

```

**Navigate to the project folder,  and run the app:**

```bash

# move to the project folder
cd download/SOSDNA main

streamlit run main.py
```

Navigate to [localhost:8501](https://localhost:8501/). You should see the app running in your broswer :)


### Use jupyter notebooks

Download the project and install the dependencie following the instructions above, and run the notebooks. Notebooks include:

* **Analysis of Individual Stages of DNA data channel.ipynb** (Corresponding to 2.2 of the paper)
* **Step-by-step simulation.ipynb** (Part 3.1 of the paper)
* **Encoding Design.ipynb** (Part 3.2 of the paper)

## Acknowledgements
Code about the DNA Fountain code is adopted from https://github.com/TeamErlich/dna-fountain.



