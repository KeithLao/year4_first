#include <iostream>
#include <vector>
#include <cmath>
#include "support_cpp/base_matrix.cpp"
#include "support_cpp/ini_seq.cpp"
#include "support_cpp/decoder.cpp"
#include "support_cpp/performance.cpp"
using namespace std;

int main() {
    // Parameter settings
    int max_iterations = 100;
    //vector<double> snr_range = {2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0};
    vector<double> snr_range = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0};
    int max_frames = 1000; // 1000 - 1min6s 10000 - 15min
    int seed = 42;
    
    // Load data
    string algorithm;
    string type = "PEGReg504x1008.txt";
    BaseMatrix H(type);

    int ini_length = H.n - H.m;
    IniSeq iniSeq(H.n, H.m, ini_length); 
    double bitrate = (double)H.m / (double)H.n;

    Decoder DD(H);
    Performance perf;
    
    // Pre-generate data
    vector<vector<int>> ini_save(max_frames, vector<int> (ini_length) );
    vector<vector<int>> messages_save(max_frames, vector<int> (H.n) );

    for (int i = 0 ; i < max_frames; i++){
        int bias = i;
        ini_save[i] = iniSeq.generateRandom(seed + bias);   // (0, 1) sequences k length
        messages_save[i] = iniSeq.Encoder(H, ini_save[i]);  // (1,-1) sequences n length
    }

    // Test for each SNR
    for (double snr : snr_range) {
        cout << "\nTesting SNR = " << snr << " dB" << endl;
        
        int frame_count = 0;
        int error_frames = 0;
        int total_errors = 0;
        int hard_total_count = 0;
        int hard_max = 0, hard_min = 1000;

        while (frame_count < max_frames) {
            // Get random information sequence
            int bias = frame_count;
            vector<int> message = messages_save[frame_count];   // (1,-1) sequences n length

            // Add AWGN noise
            vector<double> received = iniSeq.addAWGN(message, snr, seed + bias, bitrate);
            vector<double> noise = iniSeq.generate_noise(snr);

            // Hard decision on received signal
            vector<int> hard_received = iniSeq.hard_decision(received);

            // Check hard decision pefromance
            int hard_count = 0;
            for (int i = 0; i < message.size(); i++) {
                if (message[i] != hard_received[i]) {
                    hard_count++;
                    hard_total_count++;
                }
            }
            if(hard_max < hard_count){ hard_max = hard_count; }
            if(hard_min > hard_count){ hard_min = hard_count; }


            // Decoding
            algorithm = "GDBF";
            vector<int> message_after = DD.GDBF(algorithm, received, hard_received, max_iterations, perf, noise);  // Using GDBF algorithm
            //vector<int> message_after = hard_received;

            // Check if decoding is successful
            int bit_errors = 0;
            for (int i = 0; i < message.size(); i++) {
                if (message[i] != message_after[i]) {
                    bit_errors++;
                }
            }
            if (bit_errors > 0) {
                error_frames++;
                total_errors += bit_errors;
            }

            frame_count++;
        }
        
        double hard_average_error = (double)hard_total_count / frame_count;

        cout << "  Hard decesion average bit error: " << fixed         
             << setprecision(5) << hard_average_error << endl;

        cout << "  max bit error: " << hard_max 
             << ", min bit error: " << hard_min << endl;

        // Calculate final performance
        double fer = (double)error_frames / frame_count;
        double ber = (double)total_errors / (frame_count * H.n);
        
        cout << "  Final result - FER: " << scientific << fer 
             << ", BER: " << ber << endl;

        cout << "  error frames: " << error_frames << endl;
        
        perf.addResult(algorithm, snr, fer, ber);
    }
    
    // Display and save results
    perf.saveToFile(type);
    
    return 0;
}