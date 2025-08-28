// VoiceCalculator.cpp : Defines the entry point for the console application.
// This program first trains a voice recognition model for digits 0-9
// and then launches a command-line voice-operated calculator.

#include "stdafx.h" // This MUST be the very first line for Visual Studio 2010

#define _CRT_SECURE_NO_WARNINGS // Disables security warnings for fopen, etc.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <stack>
#include "conio.h"
#include <Windows.h>
#pragma comment(lib, "winmm.lib") // Link the Windows Multimedia library

using namespace std;

// ================================================================================= //
//                            CONFIGURATION & CONSTANTS                              //
// ================================================================================= //

#define LPC_ORDER 12           // Order of LPC (Linear Predictive Coding)
#define PI 3.1415926535        // More precise value of PI
#define FRAME_SIZE 320         // Number of samples per frame
#define MAX_LINE_LENGTH 1024
#define TOTAL_FILES 400        // 10 digits * 40 utterances each
#define K 32                   // Codebook size
#define DIMENSION 12           // Vector dimensions (same as LPC_ORDER)
#define DELTA 0.00001          // Convergence threshold for K-Means
#define MAX_VECTORS 15207      // Maximum number of vectors expected in the universe
#define EPSILON 0.03           // Epsilon for centroid splitting in LBG algorithm
#define N 5                    // Number of states in HMM
#define O 32                   // Number of observation symbols (same as K)
#define T_MAX 1000             // Maximum number of observations in a sequence

// ================================================================================= //
//                               GLOBAL VARIABLES                                    //
// ================================================================================= //

// --- Signal Processing ---
int totalSamples;
int TOTAL_FRAMES;
double rawSignal[22000] = {0.0};
double rawSignalData[65][FRAME_SIZE];
double autocorrelation[65][LPC_ORDER + 1];
double lpcCoefficients[65][LPC_ORDER + 1];
double cepstralCoefficients[65][LPC_ORDER + 1];
double universe[MAX_VECTORS][DIMENSION];
double codebook[K][DIMENSION];

// --- Observation Sequence Generation ---
double framestore[100][DIMENSION];
int frameonedigit;
int ansstore[100];

// --- HMM Training & Testing ---
long double A[N + 1][N + 1];
long double B[N + 1][O + 1];
long double pi[N + 1];
long double A_comp[N + 1][N + 1];
long double B_comp[N + 1][O + 1];
long double pi_comp[N + 1];
long int obs[T_MAX] = {0};
long double alpha[T_MAX][N + 1];
long double beta[T_MAX][N + 1];
long double gamma_hmm[T_MAX][N + 1]; 
long double delta[T_MAX][N + 1];
long double zai[T_MAX][N + 1][N + 1];
int shai[T_MAX][N + 1];
int qstar[T_MAX];
int T; // Number of observations in a sequence
long int observation[100];
long double prob[10] = {0};
int correct_predictions = 0;
int wrong_predictions = 0;


// ================================================================================= //
//                        PART 1: SIGNAL PROCESSING FUNCTIONS                        //
// ================================================================================= //

// Reads raw signal data from a text file, skipping the header.
void ReadSignalDataFromFile(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char tempStr[150];
    for (int skip = 0; skip < 5; skip++) { // Skip header lines
        fgets(tempStr, sizeof(tempStr), file);
    }

    totalSamples = 0; // Reset for each file
    while (!feof(file) && totalSamples < 22000) {
        if (fscanf(file, "%*d %lf", &rawSignal[totalSamples]) == 1) {
            totalSamples++;
        } else {
            break;
        }
    }
    fclose(file);
}

// Removes DC offset from the signal by subtracting the mean.
void CorrectDCOffsetInSignal() {
    double mean = 0.0;
    if (totalSamples == 0) return;
    for (int i = 0; i < totalSamples; i++) {
        mean += rawSignal[i];
    }
    mean /= totalSamples;
    for (int i = 0; i < totalSamples; i++) {
        rawSignal[i] -= mean;
    }
}

// Normalizes the signal to a range of [-5000, 5000].
void NormalizeSignalValues() {
    if (totalSamples == 0) return;
    double minVal = rawSignal[0];
    double maxVal = rawSignal[0];

    for (int i = 1; i < totalSamples; i++) {
        if (rawSignal[i] < minVal) minVal = rawSignal[i];
        if (rawSignal[i] > maxVal) maxVal = rawSignal[i];
    }

    double range = maxVal - minVal;
    if (range == 0) return; // Avoid division by zero

    for (int i = 0; i < totalSamples; i++) {
        rawSignal[i] = -5000.0 + ((rawSignal[i] - minVal) / range) * 10000.0;
    }
}

// Segments the raw signal into frames.
void SelectSteadyStateFramesFromSignal() {
    for (int frame = 0; frame < TOTAL_FRAMES; frame++) {
        for (int sample = 0; sample < FRAME_SIZE; sample++) {
            rawSignalData[frame][sample] = rawSignal[frame * FRAME_SIZE + sample];
        }
    }
}

// Applies a Hamming window to a frame to reduce spectral leakage.
void ApplyHammingWindowToFrame(int frameIndex) {
    for (int sample = 0; sample < FRAME_SIZE; sample++) {
        rawSignalData[frameIndex][sample] *= 0.54 - 0.46 * cos(2 * PI * sample / (FRAME_SIZE - 1));
    }
}

// Computes autocorrelation coefficients for a given frame.
void ComputeAutoCorrelationForFrame(int frameIndex) {
    for (int lag = 0; lag <= LPC_ORDER; lag++) {
        autocorrelation[frameIndex][lag] = 0.0;
        for (int sample = 0; sample < FRAME_SIZE - lag; sample++) {
            autocorrelation[frameIndex][lag] += rawSignalData[frameIndex][sample] * rawSignalData[frameIndex][sample + lag];
        }
    }
}

// Computes LPC coefficients using the Levinson-Durbin algorithm.
void ComputeLPCoefficientsForFrame(int frameIndex) {
    double error[LPC_ORDER + 1] = {0.0};
    double reflectionCoeffs[LPC_ORDER + 1] = {0.0};
    double predictionErrors[LPC_ORDER + 1][LPC_ORDER + 1] = {0.0};

    error[0] = autocorrelation[frameIndex][0];
    if (error[0] == 0) return;

    for (int i = 1; i <= LPC_ORDER; i++) {
        double sum = 0.0;
        for (int j = 1; j < i; j++) {
            sum += predictionErrors[j][i - 1] * autocorrelation[frameIndex][i - j];
        }
        reflectionCoeffs[i] = (autocorrelation[frameIndex][i] - sum) / error[i - 1];
        predictionErrors[i][i] = reflectionCoeffs[i];
        for (int j = 1; j < i; j++) {
            predictionErrors[j][i] = predictionErrors[j][i - 1] - reflectionCoeffs[i] * predictionErrors[i - j][i - 1];
        }
        error[i] = (1 - reflectionCoeffs[i] * reflectionCoeffs[i]) * error[i - 1];
    }

    for (int i = 1; i <= LPC_ORDER; i++) {
        lpcCoefficients[frameIndex][i] = predictionErrors[i][LPC_ORDER];
    }
}

// Computes Cepstral Coefficients from LPC coefficients.
void ComputeCepstralCoefficientsForFrame(int frameIndex) {
    cepstralCoefficients[frameIndex][0] = 0.0;
    for (int n = 1; n <= LPC_ORDER; n++) {
        cepstralCoefficients[frameIndex][n] = lpcCoefficients[frameIndex][n];
        double sum = 0.0;
        for (int k = 1; k < n; k++) {
            sum += (double)k * cepstralCoefficients[frameIndex][k] * lpcCoefficients[frameIndex][n - k];
        }
        cepstralCoefficients[frameIndex][n] += sum / (double)n;
    }
}

// Writes the computed cepstral coefficients to a file.
void WriteCepstralCoeffsToFile(const char* filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening cepstral output file");
        exit(EXIT_FAILURE);
    }
    for (int frame = 0; frame < TOTAL_FRAMES; frame++) {
        for (int coef = 1; coef <= LPC_ORDER; coef++) {
            fprintf(file, "%lf ", cepstralCoefficients[frame][coef]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

// Merges all individual cepstral files into one large "universe" file.
void mergeFiles(const char *outputFile) {
    FILE *outFile = fopen(outputFile, "w");
    if (outFile == NULL) {
        perror("Error opening universe file for writing");
        exit(EXIT_FAILURE);
    }
    char fileName[100] = "";
    for (int i = 1; i <= TOTAL_FILES; i++) {
        sprintf(fileName, "./OUTPUT_CEPSTRAL/244101003_%d.txt", i);
        FILE *inFile = fopen(fileName, "r");
        if (inFile == NULL) {
            fprintf(stderr, "Warning: Could not open input file %s. Skipping.\n", fileName);
            continue;
        }
        char line[MAX_LINE_LENGTH];
        while (fgets(line, sizeof(line), inFile) != NULL) {
            fputs(line, outFile);
        }
        fclose(inFile);
    }
    fclose(outFile);
    printf("Cepstral files merged successfully into '%s'.\n", outputFile);
}


// ================================================================================= //
//                      PART 2: VECTOR QUANTIZATION (CODEBOOK)                       //
// ================================================================================= //

// Loads the universe of vectors from the merged file.
void store_all_value_of_universe(const char *filename, double universe[][DIMENSION]) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening universe file for reading.\n");
        exit(1);
    }
    int k = 0;
    while (k < MAX_VECTORS && !feof(file)) {
        int items_read = fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &universe[k][0], &universe[k][1], &universe[k][2], &universe[k][3],
               &universe[k][4], &universe[k][5], &universe[k][6], &universe[k][7],
               &universe[k][8], &universe[k][9], &universe[k][10], &universe[k][11]);
        if(items_read == DIMENSION) k++;
    }
    fclose(file);
}

// Initializes the codebook with the single centroid of the entire universe.
void initialize_codebook(double codebook[1][DIMENSION], double universe[][DIMENSION], int M) {
    for (int i = 0; i < DIMENSION; i++) {
        double sum = 0;
        for (int j = 0; j < M; j++) {
            sum += universe[j][i];
        }
        codebook[0][i] = sum / M;
    }
}

// Doubles the codebook size by splitting each existing centroid.
void split_codebook(double current_codebook[][DIMENSION], int *current_size) {
    int new_size = *current_size * 2;
    for (int i = 0; i < *current_size; i++) {
        for (int j = 0; j < DIMENSION; j++) {
            current_codebook[i + *current_size][j] = current_codebook[i][j] * (1 + EPSILON);
            current_codebook[i][j] = current_codebook[i][j] * (1 - EPSILON);
        }
    }
    *current_size = new_size;
}

// Calculates the weighted Tokhura distance between two vectors.
double tokhura_distance(double vec1[DIMENSION], double vec2[DIMENSION]) {
    double weights[DIMENSION] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
    double dist = 0;
    for (int i = 0; i < DIMENSION; i++) {
        dist += weights[i] * pow(vec1[i] - vec2[i], 2);
    }
    return dist;
}

// Assigns each universe vector to the nearest centroid and calculates total distortion.
void form_cluster(double universe[][DIMENSION], int M, double current_codebook[][DIMENSION], int current_size, int region[], double *dist) {
    *dist = 0;
    for (int i = 0; i < M; i++) {
        double min_tok_dist = DBL_MAX;
        int nearest_centroid = -1;
        for (int j = 0; j < current_size; j++) {
            double tok_dist = tokhura_distance(universe[i], current_codebook[j]);
            if (tok_dist < min_tok_dist) {
                min_tok_dist = tok_dist;
                nearest_centroid = j;
            }
        }
        region[i] = nearest_centroid;
        *dist += min_tok_dist;
    }
}

// Updates centroids to be the mean of all vectors assigned to them.
void update_centroids(double universe[][DIMENSION], int M, double current_codebook[][DIMENSION], int current_size, int region[]) {
    double store_sum[K][DIMENSION] = {{0}};
    int current_cluster_size[K] = {0};

    for (int i = 0; i < M; i++) {
        int reg = region[i];
        if (reg != -1) {
            for (int j = 0; j < DIMENSION; j++) {
                store_sum[reg][j] += universe[i][j];
            }
            current_cluster_size[reg]++;
        }
    }

    for (int i = 0; i < current_size; i++) {
        if (current_cluster_size[i] > 0) {
            for (int j = 0; j < DIMENSION; j++) {
                current_codebook[i][j] = store_sum[i][j] / current_cluster_size[i];
            }
        }
    }
}


// ================================================================================= //
//                     PART 3: OBSERVATION SEQUENCE GENERATION                       //
// ================================================================================= //

// Reads a cepstral file for a single utterance.
void readfileforstep3(const char* filename3){
    FILE *file1 = fopen(filename3, "r");
    if (!file1) {
        printf("Error opening file for observation sequence generation: %s\n", filename3);
        exit(1);
    }
    frameonedigit = 0;
    while (frameonedigit < 100 && !feof(file1)) {
        int items_read = fscanf(file1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &framestore[frameonedigit][0], &framestore[frameonedigit][1], &framestore[frameonedigit][2], &framestore[frameonedigit][3],
               &framestore[frameonedigit][4], &framestore[frameonedigit][5], &framestore[frameonedigit][6], &framestore[frameonedigit][7],
               &framestore[frameonedigit][8], &framestore[frameonedigit][9], &framestore[frameonedigit][10], &framestore[frameonedigit][11]);
        if(items_read == DIMENSION) frameonedigit++;
    }
    fclose(file1);
}

// Calculates Tokhura distance (duplicate, can be removed).
double tokhura(double vec1[DIMENSION], double vec2[DIMENSION]) {
    return tokhura_distance(vec1, vec2);
}

// Writes the generated observation sequence to a file.
void writeobservation(const char* filename2){
    FILE *file2 = fopen(filename2, "w");
    if (file2 == NULL) {
        perror("Error opening observation output file");
        exit(EXIT_FAILURE);
    }
    for (int c = 0; c < frameonedigit; c++) {
        fprintf(file2, "%d ", ansstore[c]);
    }
    fclose(file2);
}


// ================================================================================= //
//                         PART 4: HMM TRAINING & TESTING                            //
// ================================================================================= //

// Reads an observation sequence from a file.
void readfileforstep4(const char* filename4){
    FILE *file1 = fopen(filename4, "r");
    if (!file1) {
        printf("Error opening observation file for training: %s\n", filename4);
        exit(1);
    }
    T = 0;
    while (T < T_MAX && fscanf(file1, "%ld", &obs[T + 1]) == 1) {
        T++;
    }
    fclose(file1);
}

// HMM Forward procedure (Alpha pass).
long double forward_procedure() {
    for (int i = 1; i <= N; i++)
        alpha[1][i] = pi[i] * B[i][obs[1]];

    for (int t = 1; t <= T - 1; t++) {
        for (int j = 1; j <= N; j++) {
            long double temp = 0;
            for (int i = 1; i <= N; i++)
                temp += alpha[t][i] * A[i][j];
            alpha[t + 1][j] = temp * B[j][obs[t + 1]];
        }
    }

    long double probability = 0;
    for (int i = 1; i <= N; i++)
        probability += alpha[T][i];
    return probability;
}

// HMM Backward procedure (Beta pass).
void backward_procedure() {
    for (int i = 1; i <= N; i++)
        beta[T][i] = 1;

    for (int t = T - 1; t >= 1; t--) {
        for (int i = 1; i <= N; i++) {
            long double temp = 0;
            for (int j = 1; j <= N; j++) {
                temp += A[i][j] * B[j][obs[t + 1]] * beta[t + 1][j];
            }
            beta[t][i] = temp;
        }
    }
}

// Calculates gamma values.
void calculate_gamma() {
    for (int t = 1; t <= T; t++) {
        long double denominator = 0;
        for (int j = 1; j <= N; j++) {
            denominator += alpha[t][j] * beta[t][j];
        }
        if (denominator == 0) continue;
        for (int i = 1; i <= N; i++) {
            gamma_hmm[t][i] = (alpha[t][i] * beta[t][i]) / denominator;
        }
    }
}

// Viterbi algorithm to find the most likely state sequence.
long double viterbi_algo() {
    for (int i = 1; i <= N; i++) {
        delta[1][i] = pi[i] * B[i][obs[1]];
        shai[1][i] = 0;
    }

    for (int t = 2; t <= T; t++) {
        for (int j = 1; j <= N; j++) {
            long double max_val = 0;
            int max_state = 0;
            for (int i = 1; i <= N; i++) {
                long double temp = delta[t - 1][i] * A[i][j];
                if (temp > max_val) {
                    max_val = temp;
                    max_state = i;
                }
            }
            delta[t][j] = max_val * B[j][obs[t]];
            shai[t][j] = max_state;
        }
    }

    long double Pstar = 0;
    int q_T = 0;
    for (int i = 1; i <= N; i++) {
        if (delta[T][i] > Pstar) {
            Pstar = delta[T][i];
            q_T = i;
        }
    }
    qstar[T] = q_T;
    return Pstar;
}

// Re-estimates model parameters (A, B, pi).
void re_estimation() {
    // Re-estimate pi
    for (int i = 1; i <= N; i++) {
        pi[i] = gamma_hmm[1][i];
    }

    // Re-estimate A
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            long double numerator = 0, denominator = 0;
            for (int t = 1; t <= T - 1; t++) {
                numerator += zai[t][i][j];
                denominator += gamma_hmm[t][i];
            }
            A[i][j] = (denominator == 0) ? 0 : numerator / denominator;
        }
    }

    // Re-estimate B
    for (int j = 1; j <= N; j++) {
        long double denominator = 0;
        for (int t = 1; t <= T; t++) {
            denominator += gamma_hmm[t][j];
        }
        if (denominator == 0) continue;
        for (int k = 1; k <= O; k++) {
            long double numerator = 0;
            for (int t = 1; t <= T; t++) {
                if (obs[t] == k)
                    numerator += gamma_hmm[t][j];
            }
            B[j][k] = numerator / denominator;
        }
    }
}

// Ensures model parameters remain stochastic (rows sum to 1).
void maintain_stochastic() {
    // For B matrix
    for (int i = 1; i <= N; i++) {
        long double row_sum = 0;
        for (int j = 1; j <= O; j++) {
            if (B[i][j] < 1e-30) B[i][j] = 1e-30;
            row_sum += B[i][j];
        }
        if (row_sum == 0) continue; // Avoid division by zero if all elements are zero
        for (int j = 1; j <= O; j++) {
            B[i][j] /= row_sum;
        }
    }
}

// Calculates zai values and performs re-estimation.
void baum_welch() {
    long double denominator_zai[T_MAX] = {0};
    for(int t=1; t<=T-1; ++t){
        for(int i=1; i<=N; ++i){
            for(int j=1; j<=N; ++j){
                denominator_zai[t] += alpha[t][i] * A[i][j] * B[j][obs[t+1]] * beta[t+1][j];
            }
        }
    }

    for (int t = 1; t <= T - 1; t++) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                if(denominator_zai[t] == 0) zai[t][i][j] = 0;
                else zai[t][i][j] = (alpha[t][i] * A[i][j] * B[j][obs[t + 1]] * beta[t + 1][j]) / denominator_zai[t];
            }
        }
    }
    re_estimation();
    maintain_stochastic();
}

// Main training loop for a single model, iterates until convergence.
void converge_model() {
    long double old_Pstar = -1, new_Pstar = 0;
    int countt = 0;
    do {
        countt++;
        old_Pstar = new_Pstar;
        forward_procedure();
        backward_procedure();
        calculate_gamma();
        new_Pstar = viterbi_algo();
        baum_welch();
    } while (countt < 100 && old_Pstar < new_Pstar);
}

// Saves the trained A and B matrices to a file.
void store_in_file(string file) {
    ofstream store(file);
    for (int var = 1; var <= N; var++) {
        for (int val = 1; val <= N; val++) {
            store << A[var][val] << " ";
        }
        store << endl;
    }
    store << endl;
    for (int var = 1; var <= N; var++) {
        for (int val = 1; val <= O; val++) {
            store << B[var][val] << " ";
        }
        store << endl;
    }
    store.close();
}

// Accumulates model parameters over multiple training utterances.
void adding_values() {
    for (int i = 1; i <= N; i++) pi_comp[i] += gamma_hmm[1][i];
    for (int i = 1; i <= N; i++) for (int j = 1; j <= N; j++) A_comp[i][j] += A[i][j];
    for (int i = 1; i <= N; i++) for (int j = 1; j <= O; j++) B_comp[i][j] += B[i][j];
}

// Averages the accumulated model parameters to create the final model.
void avg_values() {
    for (int i = 1; i <= N; i++) {
        pi_comp[i] /= 30.0;
        pi[i] = pi_comp[i];
    }
    for (int var = 1; var <= N; var++) {
        for (int val = 1; val <= N; val++) {
            A_comp[var][val] /= 30.0;
            A[var][val] = A_comp[var][val];
        }
    }
    for (int var = 1; var <= N; var++) {
        for (int val = 1; val <= O; val++) {
            B_comp[var][val] /= 30.0;
            B[var][val] = B_comp[var][val];
        }
    }
    maintain_stochastic();
}

// Recognizes a digit by finding the model with the highest probability.
int recognizer(int val) {
    int idx = -1;
    long double maxi = -DBL_MAX;
    for (int i = 0; i <= 9; i++) {
        if (prob[i] > maxi) {
            idx = i;
            maxi = prob[i];
        }
    }
    if(val != -1) { // -1 is used for live recognition where we don't know the answer
        printf("Test Utterance (Digit %d) -> Recognized as %d. ", val, idx);
        if (val == idx) {
            correct_predictions++;
            printf("(Correct)\n");
        } else {
            wrong_predictions++;
            printf("(Wrong)\n");
        }
    }
    return idx;
}

// Calculates the probability of an observation sequence given a model.
int find_probability(int val) {
    for (int k = 0; k <= 9; k++) {
        char model_path[100] = "";
        sprintf(model_path, "./OUTPUT_a_b_final/244101003_%d.txt", k);
        ifstream fp(model_path);
        if (!fp.is_open()) {
            fprintf(stderr, "Error: Cannot open final model file %s\n", model_path);
            continue;
        }

        for (int i = 1; i <= N; ++i) for (int j = 1; j <= N; ++j) fp >> A[i][j];
        for (int i = 1; i <= N; ++i) for (int j = 1; j <= O; ++j) fp >> B[i][j];
        fp.close();

        prob[k] = forward_procedure();
    }
    return recognizer(val);
}


// ================================================================================= //
//                      PART 5: LIVE RECORDING & CALCULATOR                          //
// ================================================================================= //

#pragma comment(lib, "winmm.lib")
short int waveIn[32050];

// Records 2 seconds of audio from the microphone.
void record_audio() {
    const int NUMPTS = 16025 * 2;
    int sampleRate = 16025;
    HWAVEIN hWaveIn;
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;
    pFormat.nChannels = 1;
    pFormat.nSamplesPerSec = sampleRate;
    pFormat.nAvgBytesPerSec = sampleRate * 2;
    pFormat.nBlockAlign = 2;
    pFormat.wBitsPerSample = 16;
    pFormat.cbSize = 0;

    waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
    WAVEHDR WaveInHdr;
    WaveInHdr.lpData = (LPSTR)waveIn;
    WaveInHdr.dwBufferLength = NUMPTS * 2;
    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;
    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    waveInStart(hWaveIn);
    
    cout << "Recording..." << flush;
    Sleep(2000); // Record for 2 seconds
    cout << " Finished." << endl;

    waveInClose(hWaveIn);

    // Save the recorded data to a temporary file
    FILE *rec_file = fopen("./RECORDING/live_recording.txt", "w");
    if(!rec_file) {
        perror("Could not create recording file");
        return;
    }
    fprintf(rec_file, "Header\nHeader\nHeader\nHeader\nHeader\n"); // Dummy header
    for (int i = 0; i < NUMPTS; i++) {
        fprintf(rec_file, "%d %d\n", i, (int)waveIn[i]);
    }
    fclose(rec_file);
}

// Full pipeline to record and recognize a single spoken digit.
int recognize_live_digit() {
    record_audio();

    // 1. Process the recorded signal
    ReadSignalDataFromFile("./RECORDING/live_recording.txt");
    CorrectDCOffsetInSignal();
    NormalizeSignalValues();
    TOTAL_FRAMES = totalSamples / FRAME_SIZE;
    if (TOTAL_FRAMES <= 0) {
        cout << "Could not detect sufficient audio. Please speak louder." << endl;
        return -1;
    }
    SelectSteadyStateFramesFromSignal();
    for (int j = 0; j < TOTAL_FRAMES; j++) {
        ApplyHammingWindowToFrame(j);
        ComputeAutoCorrelationForFrame(j);
        ComputeLPCoefficientsForFrame(j);
        ComputeCepstralCoefficientsForFrame(j);
    }
    WriteCepstralCoeffsToFile("./OUTPUTRECORDING/live_cepstral.txt");

    // 2. Generate observation sequence
    readfileforstep3("./OUTPUTRECORDING/live_cepstral.txt");
    for (int t = 0; t < frameonedigit; t++) {
        double min_dist = DBL_MAX;
        int nearest = 0;
        for (int u = 0; u < K; u++) {
            double tok_dist = tokhura(framestore[t], codebook[u]);
            if (tok_dist < min_dist) {
                min_dist = tok_dist;
                nearest = u;
            }
        }
        ansstore[t] = nearest + 1;
    }
    writeobservation("./RECORDINGSEQUENCE/live_observation.txt");

    // 3. Recognize using HMMs
    readfileforstep4("./RECORDINGSEQUENCE/live_observation.txt");
    return find_probability(-1); // -1 indicates live recognition
}

// **** NEW AND IMPROVED CALCULATOR LOGIC ****
void run_calculator() {
    stack<double> operands;
    stack<char> operators;
    string current_number_str = "";
    bool expecting_number = true; // The calculator always starts by expecting a number

    cout << "\n\n--- Voice Calculator Initialized ---" << endl;
    cout << "====================================" << endl;

    while (true) {
        cout << "\n------------------------------------" << endl;
        // Display the current state of the calculation
        cout << "Current Entry: " << (current_number_str.empty() ? "..." : current_number_str) << endl;

        cout << "CHOOSE ACTION: Speak '1' for Digit | '0' for Operator | '2' for Equals" << endl;
        int mode_choice = recognize_live_digit();
        cout << "Recognized Action: " << mode_choice << endl;

        if (mode_choice == 1) { // User wants to enter a digit
            if (!expecting_number) {
                cout << "Illegal Syntax! Expected an operator but you chose to enter a digit. Please try again." << endl;
                continue;
            }
            cout << "Speak a digit (0-9):" << endl;
            int digit = recognize_live_digit();
            cout << "Recognized Digit: " << digit << endl;
            if (digit >= 0 && digit <= 9) {
                current_number_str += to_string((long long)digit);
                expecting_number = false; // After entering a digit, we expect an operator or another digit
            } else {
                cout << "Invalid digit recognized. Please try again." << endl;
            }

        } else if (mode_choice == 0) { // User wants to enter an operator
            if (expecting_number || current_number_str.empty()) {
                cout << "Illegal Syntax! You must enter a number before an operator." << endl;
                continue;
            }
            operands.push(atof(current_number_str.c_str()));
            current_number_str = "";

            cout << "Speak an operator: '1' for + | '2' for - | '3' for * | '4' for /" << endl;
            int op_choice = recognize_live_digit();
            cout << "Recognized Operator Command: " << op_choice << endl;

            if (op_choice >= 1 && op_choice <= 4) {
                char op = ' ';
                if (op_choice == 1) op = '+';
                else if (op_choice == 2) op = '-';
                else if (op_choice == 3) op = '*';
                else if (op_choice == 4) op = '/';
                operators.push(op);
                expecting_number = true; // After an operator, we expect a number
            } else {
                cout << "Invalid operator recognized. Please try again." << endl;
                // Give the user their number back if the operator was bad
                current_number_str = to_string((long double)operands.top());
                operands.pop();
            }

        } else if (mode_choice == 2) { // User wants to calculate equals
            if (operators.empty() || current_number_str.empty()) {
                cout << "Illegal Syntax! Incomplete expression." << endl;
                continue;
            }
            operands.push(atof(current_number_str.c_str()));
            
            // This logic only handles a single operation, e.g., A + B
            if (operands.size() == 2 && operators.size() == 1) {
                double num2 = operands.top(); operands.pop();
                double num1 = operands.top(); operands.pop();
                char op = operators.top(); operators.pop();
                double result = 0;

                if (op == '+') result = num1 + num2;
                else if (op == '-') result = num1 - num2;
                else if (op == '*') result = num1 * num2;
                else if (op == '/') {
                    if (num2 == 0) { cout << "Error: Division by zero." << endl; result = 0; }
                    else result = num1 / num2;
                }
                
                cout << "====================================" << endl;
                cout << "FINAL RESULT: " << result << endl;
                cout << "====================================" << endl;
                return; // Exit calculator
            } else {
                cout << "Calculation error: This calculator only supports single operations (e.g., 12 + 34)." << endl;
                // Clear stacks to restart
                while(!operands.empty()) operands.pop();
                while(!operators.empty()) operators.pop();
                current_number_str = "";
                expecting_number = true;
            }
        } else {
            cout << "Invalid action. Please try again." << endl;
        }
    }
}


// ================================================================================= //
//                                  MAIN FUNCTION                                    //
// ================================================================================= //

int main(int argc, char* argv[]) {
    
    cout << "=======================================================" << endl;
    cout << "      VOICE RECOGNITION SYSTEM & CALCULATOR" << endl;
    cout << "=======================================================" << endl;
    
    // -----------------------------------------------------------------
    // STEP 1: Process raw audio files into Cepstral Coefficients
    // -----------------------------------------------------------------
    cout << "\nSTEP 1: Processing audio files into Cepstral Coefficients..." << endl;
    char inputFilePath[100] = "";
    char outputFilePath[100] = "";
    for (int v = 0; v < 10; v++) { // For each digit 0-9
        for (int i = 1; i <= 40; i++) { // For each utterance 1-40
            sprintf(inputFilePath, "./txt/244101003_E_%d_%d.txt", v, i);
            ReadSignalDataFromFile(inputFilePath);
            CorrectDCOffsetInSignal();
            NormalizeSignalValues();
            TOTAL_FRAMES = totalSamples / FRAME_SIZE;
            if (TOTAL_FRAMES <= 0) continue;
            SelectSteadyStateFramesFromSignal();
            for (int j = 0; j < TOTAL_FRAMES; j++) {
                ApplyHammingWindowToFrame(j);
                ComputeAutoCorrelationForFrame(j);
                ComputeLPCoefficientsForFrame(j);
                ComputeCepstralCoefficientsForFrame(j);
            }
            int q = 40 * v + i;
            sprintf(outputFilePath, "./OUTPUT_CEPSTRAL/244101003_%d.txt", q);
            WriteCepstralCoeffsToFile(outputFilePath);
        }
    }
    mergeFiles("UNIVERSE.csv");
    cout << "STEP 1 finished." << endl;

    // -----------------------------------------------------------------
    // STEP 2: Generate Codebook using K-Means/LBG
    // -----------------------------------------------------------------
    cout << "\nSTEP 2: Generating Codebook from Universe..." << endl;
    int M = MAX_VECTORS;
    int region[MAX_VECTORS];
    store_all_value_of_universe("./UNIVERSE.csv", universe);
    initialize_codebook(codebook, universe, M);
    int current_size = 1;
    while (current_size < K) {
        split_codebook(codebook, &current_size);
        double distortion = 0, avg_dist;
        int m = 0;
        do {
            distortion = (m > 0) ? avg_dist : 0;
            form_cluster(universe, M, codebook, current_size, region, &avg_dist);
            avg_dist /= M;
            update_centroids(universe, M, codebook, current_size, region);
            m++;
        } while (fabs(avg_dist - distortion) > DELTA);
    }
    FILE *file = fopen("codebook.csv", "w");
    if (file) {
        for (int q = 0; q < K; q++) {
            for (int w = 0; w < DIMENSION; w++) {
                fprintf(file, "%lf ", codebook[q][w]);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }
    cout << "STEP 2 finished. Codebook saved to codebook.csv" << endl;

    // -----------------------------------------------------------------
    // STEP 3: Generate Observation Sequences for each utterance
    // -----------------------------------------------------------------
    cout << "\nSTEP 3: Generating Observation Sequences..." << endl;
    char inputFilePath3[100] = "";
    char outputFilePath3[100] = "";
    for (int q = 1; q <= 400; q++) {
        sprintf(inputFilePath3, "./OUTPUT_CEPSTRAL/244101003_%d.txt", q);
        readfileforstep3(inputFilePath3);
        for (int t = 0; t < frameonedigit; t++) {
            double min_dist = DBL_MAX;
            int nearest = 0;
            for (int u = 0; u < K; u++) {
                double tok_dist = tokhura(framestore[t], codebook[u]);
                if (tok_dist < min_dist) {
                    min_dist = tok_dist;
                    nearest = u;
                }
            }
            ansstore[t] = nearest + 1;
        }
        sprintf(outputFilePath3, "./OUTPUT_OBSERVATION/244101003_%d.txt", q);
        writeobservation(outputFilePath3);
    }
    cout << "STEP 3 finished." << endl;

    // -----------------------------------------------------------------
    // STEP 4: Train HMM for each digit
    // -----------------------------------------------------------------
    cout << "\nSTEP 4: Training Hidden Markov Models for each digit..." << endl;
    char inputFilePath4[100] = "";
    char outputFilePath4[100] = "";
    char finalModelPath[100] = "";
    for (int z = 0; z < 10; z++) { // For each digit
        cout << "  Training model for digit " << z << "..." << endl;
        // Reset accumulators
        memset(A_comp, 0, sizeof(A_comp));
        memset(B_comp, 0, sizeof(B_comp));
        memset(pi_comp, 0, sizeof(pi_comp));

        for (int y = 1; y <= 30; y++) { // Using 30 utterances for training
            int num = z * 40 + y;
            sprintf(inputFilePath4, "./OUTPUT_OBSERVATION/244101003_%d.txt", num);
            readfileforstep4(inputFilePath4);

            // Initialize A, B, pi for each training run
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) A[i][j] = 0;
                if (i == 1) pi[i] = 1.0; else pi[i] = 0.0;
                A[i][i] = (i == N) ? 1.0 : 0.8;
                if (i < N) A[i][i+1] = 0.2;
                for (int j = 1; j <= O; j++) B[i][j] = 1.0 / O;
            }
            
            converge_model();
            adding_values();
        }
        avg_values();
        sprintf(finalModelPath, "./OUTPUT_a_b_final/244101003_%d.txt", z);
        store_in_file(finalModelPath);
    }
    cout << "STEP 4 finished. All models trained." << endl;

    // -----------------------------------------------------------------
    // STEP 5: Test the models on unseen data
    // -----------------------------------------------------------------
    cout << "\nSTEP 5: Testing models on unseen data..." << endl;
    char testFilePath[100] = "";
    for (int i = 0; i <= 9; i++) {
        for (int k = 31; k <= 40; k++) { // Using last 10 utterances for testing
            int num = i * 40 + k;
            sprintf(testFilePath, "./OUTPUT_OBSERVATION/244101003_%d.txt", num);
            readfileforstep4(testFilePath);
            find_probability(i);
        }
    }
    cout << "STEP 5 finished." << endl;
    cout << "\n----------------------------------------" << endl;
    cout << "          TESTING ACCURACY" << endl;
    cout << "----------------------------------------" << endl;
    cout << "Correctly Recognized: " << correct_predictions << endl;
    cout << "Incorrectly Recognized: " << wrong_predictions << endl;
    double accuracy = (correct_predictions + wrong_predictions == 0) ? 0 : (double)correct_predictions / (correct_predictions + wrong_predictions) * 100.0;
    cout << "Accuracy: " << fixed << setprecision(2) << accuracy << "%" << endl;
    cout << "----------------------------------------" << endl;


    // -----------------------------------------------------------------
    // STEP 6: Launch the Live Voice Calculator
    // -----------------------------------------------------------------
    run_calculator();

    return 0;
}
