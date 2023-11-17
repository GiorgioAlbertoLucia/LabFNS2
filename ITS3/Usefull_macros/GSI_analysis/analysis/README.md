# README for 'analysis' Folder

This README provides an overview and usage instructions for the 'analysis' folder. This folder contains scripts intended to be used for the analysis of data measured with a silicon pixel sensor. A brief description of the experiment is also provided. In principle, this code can be used for any sufficiently similar experiment.

## Table of Contents

1. [Introduction](#introduction)
    - [The experimemt](#the-experiment)
2. [Prerequisites](#prerequisites)
3. [Included Scripts](#included-scripts)
   - [customAnalysis.py](#customanalysispy)
   - [mergeAnalysis.py](#mergeanalysispy)

---

## Introduction

The 'analysis' folder is written for data analysis with a pixel silicon sensor. It provides a set of Python scripts and classes designed to streamline various data processing tasks, making it easier to work with the sensor data. This README will guide you through the folder's structure, the prerequisites, and how to utilize the included scripts effectively.

### The experiment

This code is used to study data from a measurement with a Analog Pixel test Structure. A ^{55}Fe source was used.
The APTS is a Monolithic Active Pixel Sensor, which is a silicon sensor that has the readout electronics integrated on the sensitive silicon wafer. The measurements were taken by exposing the sensor to the source and recording an event each time the signal from a pixel would go over threshold. More on the structure of the sensor, on the physical processes the source undergoes and on the spectrum expected and observed can be found [here](www.add-presentation-link.com). Results of the analysis are also presented there.

Some things you might want to keep in mind when using the scripts:
- the experiment used an internal trigger using a threshold value. In order to obtain valuable informations from the code, some threshold value should be provided;
- the sensor was a 4x4 pixel matrix. Although the number of pixels in a matrix is not an issue (the code will automatically adapt to different matrices), the information should have some specific format and information for a matrix is expected. If not working with pixel sensor, you might want to customize your dataset in order to be properly uploaded;
- the code was written to be integrated in the [APTS software](www.ass-apts-software.com). Although everything works independently, data produced with that software are expected. You might want to take a look at how the .npy files are produced (some basic information on their structure is provided in this guide as well).

---

## Prerequisites

Before using the scripts in 'analysis', ensure you have the following prerequisites:

- Python packages: `yaml`, `uproot`, `cv2`, `numpy`, `pandas`, `scipy`, `alive_progress`.
- ROOT installation: A ROOT installation is required (default or custom). Verify that all used ROOT classes are included.

---

## Included scripts

The scripts included perform analysis using classes and function from the 'utils' folder. They are needed for a proper operation of these scripts. These are also an example to write your own analysis scripts and incorporate the classes and functions from 'utils' in it.

### customAnalysis.py

**Description**: This script can perform different data analysis steps in sequence and lets you use all the previously presented classes. Requires a yaml configuration file.

**Usage**: Import the `CustomAnalysis` class into your Python script for custom data analysis. Refer to the script for usage details and examples.

### mergeAnalysis.py

**Description**: This script performs a complete analysis of a spectrum obtained by merging different datasets. It is possible to perform a voltage to energy conversion inside by explicitly providing the conversion values (the linear relation Voltage = A + B * Energy is assumed, A, B should be provided for all the datasets). Requires a yaml configuration file. The basic structure of the configuration file is explained in the README inside the 'configs' folder.

**Usage**: Import the `MergeAnalysis` class into your Python script for data merging and analysis. Refer to the script for usage details and examples. The basic structure of the configuration file is explained in the README inside the 'configs' folder.

---