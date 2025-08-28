# Voice Recognition Calculator

## Overview
This project is a complete, **voice recognition system built in C++**. It is designed to recognize spoken digits (0–9) and function as a simple, voice-operated calculator.

The application operates in two main phases:
1. **Training Phase**  
   - Processes a dataset of pre-recorded audio files for each digit.  
   - Extracts key features from the audio.  
   - Generates a universal codebook of sound characteristics.  
   - Trains a distinct Hidden Markov Model (HMM) for each digit (0–9).  

2. **Live Calculator Phase**  
   - After training and testing are complete, the program launches an interactive command-line calculator.  
   - Uses the computer's microphone to record the user's voice.  
   - Recognizes the spoken command using the trained HMMs.  
   - Performs the requested calculation.  

---

## ✨ Features
- **End-to-End Speech Recognition**: From raw audio signal to digit recognition.  
- **Signal Processing**: DC offset correction, normalization, and Hamming windowing.  
- **Advanced Feature Extraction**: Linear Predictive Coding (LPC) and Cepstral Coefficient analysis.  
- **Vector Quantization**: LBG algorithm for efficient speech vector representation.  
- **Hidden Markov Models (HMM)**: Training, testing, and recognition.  
- **Interactive Voice Calculator**: Command-line interface for multi-digit input and chained operations (`+`, `-`, `*`, `/`).  

---

## 🛠 Technology Stack
- **Language**: C++  
- **Compiler**: Visual Studio 2010 (tested)  
- **Core Algorithms**:  
  - Linear Predictive Coding (LPC)  
  - Levinson–Durbin Recursion  
    

--
