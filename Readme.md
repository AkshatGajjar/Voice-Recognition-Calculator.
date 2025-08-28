Voice Recognition Calculator
Overview
This project is a complete, voice recognition system built in C++. It is designed to recognize spoken digits (0-9) and function as a simple, voice-operated calculator.

The application operates in two main phases:

Training Phase: The program first processes a dataset of pre-recorded audio files for each digit. It extracts key features from the audio, generates a universal codebook of sound characteristics, and then trains a distinct Hidden Markov Model (HMM) for each digit (0 through 9).

Live Calculator Phase: After training and testing are complete, the program launches an interactive command-line calculator. It uses the computer's microphone to record the user's voice, recognizes the spoken command using the trained HMMs, and performs the requested calculation.

Features
End-to-End Speech Recognition: Implements the full pipeline from raw audio signal to digit recognition.

Signal Processing: Includes DC offset correction, normalization, and Hamming windowing.

Advanced Feature Extraction: Uses Linear Predictive Coding (LPC) and Cepstral Coefficient analysis to extract robust features from speech.

Vector Quantization: Employs the LBG algorithm to create a codebook for efficient representation of speech vectors.

Hidden Markov Models (HMM): Trains, tests, and uses HMMs for the core recognition task.

Interactive Voice Calculator: A command-line interface that allows for multi-digit number entry and chained arithmetic operations (+, -, *, /) using voice commands.

Technology Stack
Language: C++

Compiler: Built and tested with Visual Studio 2010

Core Algorithms: Linear Predictive Coding (LPC), Levinson-Durbin Recursion, K-Means/LBG for Vector Quantization, and the Baum-Welch algorithm for HMM training.

IDE: Microsoft Visual Studio 2010 (or a later version with the C++ Desktop Development workload).

Setup and Installation
To get the project running, you must set up the directory structure and data correctly.

1. Create the Project Folder Structure
The program requires a specific set of folders to read data from and write its output to. In your main project directory, find the Debug folder (or create it if it doesn't exist). Inside the Debug folder, create the following subfolders:

txt

OUTPUT_CEPSTRAL

OUTPUT_OBSERVATION

OUTPUT_a_b_final

RECORDING

OUTPUTRECORDING

RECORDINGSEQUENCE

2. Add the Training Data
Place your 400 audio training files inside the txt folder.

The files must be named in the format 244101003_E_D_U.txt, where D is the digit (0-9) and U is the utterance number (1-40).

Example: 244101003_E_0_1.txt, 244101003_E_0_2.txt, ..., 244101003_E_9_40.txt.

How to Build and Run
Open the Project: Open your project in Visual Studio 2010.

Add the Code: Copy the entire C++ code into your main source file (e.g., Voice_calculator.cpp).

Build the Solution: Press F7 or go to the Build menu and select Build Solution. This will compile the code and create Voice_calculator.exe inside the Debug folder.

Run the Program:

Navigate to the Debug folder (C:\SP\Voice_calculator\Debug\).

Ensure all the required subfolders (like txt, OUTPUT_CEPSTRAL, etc.) are present.

Double-click on Voice_calculator.exe to run it.

The program will first execute the lengthy training and testing process. Once complete, the interactive calculator will start automatically in the same command window.

How to Use the Calculator
The calculator uses specific spoken digits as commands for operations.

Speak Digits (0, 6, 7, 8, 9): Use these to form numbers. You can speak them sequentially to create multi-digit numbers (e.g., speaking "8" then "9" creates the number 89).

Speak '1': Enters the PLUS (+) operator.

Speak '2': Enters the MINUS (-) operator.

Speak '3': Enters the MULTIPLY (*) operator.

Speak '4': Enters the DIVIDE (/) operator.

Speak '5': Calculates the final EQUALS (=) result.

The calculator will display the current operation on screen and update it after each command