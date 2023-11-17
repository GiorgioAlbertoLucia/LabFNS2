# README for 'utils' Folder

This README provides an overview and usage instructions for the 'utils' folder. This folder contains classes intended to be used for the analysis of data measured with a silicon pixel sensor. A brief description of the experiment is also provided. In principle, this code can be used for any sufficiently similar experiment.

## Table of Contents

1. [Introduction](#introduction)
    - [The experimemt](#the-experiment)
2. [Prerequisites](#prerequisites)
3. [Included Scripts](#included-scripts)
   - [dataPreprocessor.py](#datapreprocessorpy)
   - [patternIdentifier.py](#patternidentifierpy)
   - [spectrumAnalyzer.py](#spectrumanalyzerpy)
   - [energyCal.py](#energycalpy)
   - [peakMonitoring.py](#peakmonitoringpy)
   - [email_manager.py](#email_managerpy)


---

## Introduction

The 'utils' folder is written for data analysis with a pixel silicon sensor. It provides a set of Python scripts and classes designed to streamline various data processing tasks, making it easier to work with the sensor data. This README will guide you through the folder's structure, the prerequisites, and how to utilize the included scripts effectively.

### The experiment

This code is used to study data from a measurement with a Analog Pixel test Structure. A ^{55}Fe source was used.
The APTS is a Monolithic Active Pixel Sensor, which is a silicon sensor that has the readout electronics integrated on the sensitive silicon wafer. The measurements were taken by exposing the sensor to the source and recording an event each time the signal from a pixel would go over threshold. More on the structure of the sensor, on the physical processes the source undergoes and on the spectrum expected and observed can be found [here](www.add-presentation-link.com). Results of the analysis are also presented there.

Some things you might want to keep in mind when using the scripts:
- the experiment used an internal trigger using a threshold value. In order to obtain valuable informations from the code, some threshold value should be provided;
- the sensor was a 4x4 pixel matrix. Although the number of pixels in a matrix is not an issue (the code will automatically adapt to different matrices), the information should have some specific format and information for a matrix is expected. If not working with pixel sensor, you might want to customize your dataset in order to be properly uploaded;
- the code was written to be integrated in the [APTS software](www.ass-apts-software.com). Although everything works independently, data produced with that software are expected. You might want to take a look at how the .npy files are produced (some basic information on their structure is provided in this guide as well).

---

## Prerequisites

Before using the scripts in 'utils', ensure you have the following prerequisites:

- Python packages: `yaml`, `uproot`, `cv2`, `numpy`, `pandas`, `scipy`, `alive_progress`.
- ROOT installation: A ROOT installation is required (default or custom). Verify that all used ROOT classes are included.

---

## Included Scripts

### dataPreprocessor.py

**Description**: This script contains the `NpyFileManager` class, which loads sensor data, converts values, and creates a TTree with various information. The loaded file should be a .npy file. It is expected to be a n x r x c x t array, where n is the number of events recorded, r is the number of rows  and c is the number of columns of the pixel matrix, t is the number of timeframes recorded. 
NOTE: Since the PatternIdentifier class is not ready yet, a more naive approach is used to define cluster shapes in this class. 

**Usage**: To use this script, import the `NpyFileManager` class into your Python script. Detailed usage instructions and examples can be found in the script itself.

### patternIdentifier.py

**Description**: This script contains the `PatternIdentifier` class, designed for classifying different cluster shapes in the dataset. WORK IN PROGRESS - DO NOT USE AS IT IS!

**Usage**: Import the `PatternIdentifier` class into your Python script for cluster shape classification. Refer to the script for usage details and examples.

### spectrumAnalyzer.py

**Description**: This script contains the `SpectrumAnalyzer` class, which loads and analyzes data distributions from a TTree branch. Different structures can be analysed and fitted at the same time. Configurations have to be set from a configuration file. The basic structure of the configuration file is explained in the README inside the 'configs' folder.

**Usage**: Import the `SpectrumAnalyzer` class into your Python script for data distribution analysis. Refer to the script for usage details and examples.

### energyCal.py

**Description**: This script contains the `EnergyCalibration` class, which draws energy calibration curves from CSV data. The calibration is drawn as a n order polynomial, where n is decided by the user. All configurations are to be set from a yaml configuration file. The basic structure of the configuration file is explained in the README inside the 'configs' folder.

**Usage**: Import the `EnergyCalibration` class into your Python script for energy calibration. Refer to the script for usage details and examples.

### peakMonitoring.py

**Description**: This script contains the `PeakMonitoring` class, used to monitor data evolution over time, fit data subsets, and track position changes. All configurations are to be set from a yaml configuration file. The basic structure of the configuration file is explained in the README inside the 'configs' folder.

**Usage**: Import the `PeakMonitoring` class into your Python script for data monitoring. Refer to the script for usage details and examples.

### email_manager.py

**Description**: This script provides a simple function for sending notification emails based on information from `email_info.py`.

**Usage**: Import the `send_email` function into your Python script for email notifications. Refer to the script for usage details and examples.

---

## Integration with Analysis

To effectively integrate the 'utils' scripts into your analysis workflow, consider the following steps:

1. Import the required classes from 'utils' into your analysis scripts.
2. Use configuration files located in the 'config' folder to customize script behavior.
3. Follow the examples and usage instructions provided in each script's documentation.
4. Monitor and adapt the scripts to your specific analysis needs, especially when working with ROOT classes and data.

---

