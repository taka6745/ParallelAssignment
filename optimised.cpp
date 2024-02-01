#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include <memory>
#include <cstdio>
int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = {0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
#define encode(ch) code[ch - 'A']
#define LEN 6
#define AA_NUMBER 20
#define EPSILON 1e-010

void Init() {
    M2 = 1;
    for (int i = 0; i < LEN - 2; i++)    // M2 = AA_NUMBER ^ (LEN-2);
        M2 *= AA_NUMBER;
    M1 = M2 * AA_NUMBER;        // M1 = AA_NUMBER ^ (LEN-1);
    M = M1 * AA_NUMBER;            // M  = AA_NUMBER ^ (LEN);
}

class Bacteria {
private:
    std::vector<long> second;
    long one_l[AA_NUMBER];
    long indexs;
    long total;
    long total_l;
    long complement;

public:
    std::vector<long> vector;

    Bacteria(const std::vector<char> &data) : vector(M, 0), second(M1, 0) {
        memset(one_l, 0, AA_NUMBER * sizeof(long));
        total = 0;
        total_l = 0;
        complement = 0;
        LoadData(data);
    }

    ~Bacteria() {}

    void LoadData(const std::vector<char> &data) {
        InitVectors();

        size_t data_size = data.size();
        for (size_t i = 0; i < data_size; ++i) {
            char ch = data[i];
            if (ch == '>') {
                while (data[++i] != '\n' && i < data_size); // skip rest of line
                char buffer[LEN - 1];
                for (int j = 0; j < LEN - 1 && i < data_size; ++j)
                    buffer[j] = data[++i];
                init_buffer(buffer);
            } else if (ch != '\n') {
                cont_buffer(ch);
            }
        }
    }

    double stochastic_compute(long i) {
        double p1 = (double) second[i / AA_NUMBER] / (total + complement);
        double p2 = (double) one_l[i % AA_NUMBER] / total_l;
        double p3 = (double) second[i % M1] / (total + complement);
        double p4 = (double) one_l[i / M1] / total_l;
        return total * (p1 * p2 + p3 * p4) / 2;
    }

private:
    void InitVectors() {}

    void init_buffer(char *buffer) {
        complement++;
        indexs = 0;
        for (int i = 0; i < LEN - 1; i++) {
            short enc = encode(buffer[i]);
            one_l[enc]++;
            total_l++;
            indexs = indexs * AA_NUMBER + enc;
        }
        second[indexs]++;
    }

    void cont_buffer(char ch) {
        short enc = encode(ch);
        one_l[enc]++;
        total_l++;
        long index = indexs * AA_NUMBER + enc;
        vector[index]++;
        total++;
        indexs = (indexs % M2) * AA_NUMBER + enc;
        second[indexs]++;
    }
};

void ReadInputFile(const char *input_name, std::vector<std::vector<char>> &bacteria_data) {
    FILE *input_file = fopen(input_name, "r");
    if (input_file == nullptr) {
        fprintf(stderr, "Error: failed to open file %s\n", input_name);
        exit(EXIT_FAILURE);
    }

    if (fscanf(input_file, "%d", &number_bacteria) != 1) {
        fprintf(stderr, "Error: failed to read number of bacteria\n");
        fclose(input_file);
        exit(EXIT_FAILURE);
    }
    bacteria_name = new char *[number_bacteria];
    bacteria_data.resize(number_bacteria);

    for (long i = 0; i < number_bacteria; i++) {
        char name[10];
        if (fscanf(input_file, "%9s", name) != 1) { // Ensure buffer is not overrun
            fprintf(stderr, "Error: failed to read bacteria name\n");
            fclose(input_file);
            exit(EXIT_FAILURE);
        }
        bacteria_name[i] = new char[20];
        snprintf(bacteria_name[i], 20, "data/%s.faa", name);

        FILE *bacteria_file = fopen(bacteria_name[i], "r");
        if (bacteria_file == nullptr) {
            fprintf(stderr, "Error: failed to open file %s\n", bacteria_name[i]);
            exit(EXIT_FAILURE);
        }

        fseek(bacteria_file, 0, SEEK_END);
        size_t file_size = ftell(bacteria_file);
        fseek(bacteria_file, 0, SEEK_SET);
        bacteria_data[i].resize(file_size);
        if (fread(bacteria_data[i].data(), 1, file_size, bacteria_file) != file_size) {
            fprintf(stderr, "Error: failed to read bacteria data from %s\n", bacteria_name[i]);
            fclose(bacteria_file);
            exit(EXIT_FAILURE);
        }
        fclose(bacteria_file);
    }

    fclose(input_file);
}

double CompareBacteria(Bacteria *b1, Bacteria *b2) {
    double correlation = 0;
    double vector_len1 = 0;
    double vector_len2 = 0;

    #pragma omp parallel for reduction(+:correlation, vector_len1, vector_len2)
    for (long i = 0; i < M; i++) {
        double stochastic1 = b1->stochastic_compute(i);
        double t1 = (stochastic1 > EPSILON) ? (b1->vector[i] - stochastic1) / stochastic1 : 0;
        vector_len1 += (t1 * t1);

        double stochastic2 = b2->stochastic_compute(i);
        double t2 = (stochastic2 > EPSILON) ? (b2->vector[i] - stochastic2) / stochastic2 : 0;
        vector_len2 += (t2 * t2);

        correlation += t1 * t2;
    }

    vector_len1 = sqrt(vector_len1);
    vector_len2 = sqrt(vector_len2);

    correlation = correlation / (vector_len1 * vector_len2);
    return correlation;
}

void CompareAllBacteria(std::vector<std::unique_ptr<Bacteria>> &bacteria) {
    for (int i = 0; i < number_bacteria - 1; i++) {
        for (int j = i + 1; j < number_bacteria; j++) {
            double correlation = CompareBacteria(bacteria[i].get(), bacteria[j].get());
            printf("%03d %03d -> %.10lf\n", i, j, correlation);
        }
    }
}

int main(int argc, char *argv[]) {

    time_t t1 = time(NULL);
    Init();

    std::vector<std::vector<char>> bacteria_data;
    ReadInputFile("list.txt", bacteria_data);

    std::vector<std::unique_ptr<Bacteria>> bacteria(number_bacteria);
    #pragma omp parallel for
    for (int i = 0; i < number_bacteria; i++)
        bacteria[i] = std::make_unique<Bacteria>(bacteria_data[i]);

    CompareAllBacteria(bacteria);

    time_t t2 = time(NULL);
    printf("time elapsed: %ld seconds\n", t2 - t1);
    return 0;
}
