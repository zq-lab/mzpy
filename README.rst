# mzpy

Mzpy is a comprehensive toolkit for LC-MS data processing and metabolomics-related tasks. The toolkit includes a variety of Python scripts, each focusing on specific aspects of data handling, enrichment analysis, visualization, and integration with external databases. This Toolkit requires Python version 3.11. Please ensure that you have this version installed before setting up the environment.

## Environment Setup

### 1. **Using Conda (environment.yml)**
   ```bash
   conda env create -f environment.yml
   ```

### 2. **Using Pip (requirements.txt)**
   ```bash
   pip install -r requirements.txt
   ```

## Scripts Overview

### 1. **hmdb.py: HMDB Data Processor**
   - **Functionality:** Processes Human Metabolome Database (HMDB) XML data sheets.
   - **Key Features:**
     - Extraction of physiological effects and disease information.
     - Collection of compound details (e.g., name, SMILES, InChIKey).
     - Conversion of XML data into Pandas DataFrames.
     - Enrichment analysis for disease data.

### 2. **metab.py: Metabolomics Data Module**
   - **Functionality:** Handles metabolomics data, including processing, imputation, and enrichment analysis.
   - **Key Classes:**
     - `Metab`: Subclass of Pandas DataFrame for storing and manipulating metabolomics data.
     - `Enrichment`: Subclass for conducting enrichment analysis.
     - `RaMP`: Class for interacting with the Rat Metabolomic Portal (RaMP) database.

### 3. **mzFrame.py: Mass Spectrometry Data Processor**
   - **Functionality:** Focuses on processing mass spectrometry data, particularly MS/MS data.
   - **Key Features:**
     - Conversion of profile-style data to centroid-style.
     - Handling and analyzing precursor ions.

### 4. **NP_classify.py: Natural Products Classification**
   - **Functionality:** Obtains classification information for natural products.
   - **Key Features:**
     - Classify natural products with [NP-Classifier](https://npclassifier.ucsd.edu/)

### 5. **plotfine.py: Plotting Toolkit**
   - **Functionality:** Offers a versatile toolkit for creating bioinformatics plots with a unified style.
   - **Key Features:**
     - Bubble plots, lollipop plots, PCA plots, volcano plots, Venn diagrams.
     - Customizable parameters for each plot type.

### 6. **pubchem.py: PubChem Compound Information Finder**
   - **Functionality:** Finds compound information based on PubChem Compound ID, compound name, InChIKey, or other identifiers.
   - **Key Features:**
     - Utilizes the PubChem REST API for data retrieval.

## Contributors and Acknowledgments

- This toolkit was developed and is maintained by Zhang Qiang.
- Contributions and bug reports are welcome. Feel free to create issues or submit pull requests.

We hope this guide helps you set up the necessary environment for utilizing the Metabolomics Data Processing Toolkit. If you encounter any issues or have suggestions, please don't hesitate to reach out.