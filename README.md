# Numerical-model

This repository contains a Python-based model for analyzing the behaviour of radio emission from protostellar jets. The numerical model incorporates both thermal free–free and non-thermal synchrotron emission mechanisms within the jet geometry. This is designed to explore the effect of various parameters on the turnover frequencies and the radio spectral indices (between 10 MHz and 300 GHz) of protostellar jets, which are the primary observational features. For showcasing this model here, I have included a simple case with a sample parameter-set.  
If you're interested in understanding more about the effects of the parameters on the model and its applications, please view the following publication:  
[**Link to Publication**](https://academic.oup.com/mnras/article/514/3/3709/6577131)

Additionally, for easier understanding of the model, a simple presentation is included in the repository:  
- **`presentation.pdf`** - A presentation that explains the model and its functionality in detail.

The main script (`main.py`) runs the model, and various supporting scripts handle calculations and plotting.  

## Folder Structure  

📁 **Project Root**  
├── `main.py` - The main script to run the program  
├── `constants.py` - Contains a sample set of free parameters for the model  
├── `plot.py` - Generates a plot based on the parameters in `constants.py`  
├── `output.eps` - Output plot saved in PNG format  
├── `TotalFlux_x_const.eps` - Output plot saved in EPS format  
├── `total_flux_qxprime0.txt` - Sample output from `main.py`  
└── Other `.py` files - Helper functions called within `main.py`  

## Usage & Notes  

Since some functions are private, **`main.py` cannot be executed directly**.  
However, you can use `total_flux_qxprime0.txt` as input for `plot.py` to visualize a sample model:  

```
python plot.py
```

## Access to Private Functions  

Some files are excluded from this repository for privacy reasons.  
If you need access to the complete code, please feel free to reach out.  
