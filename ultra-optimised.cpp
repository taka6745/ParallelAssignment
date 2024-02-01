#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <memory>

int number_bacteria;
char** bacteria_name;
alignas(64) long M, M1, M2;
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
    std::unordered_map<long, long> second;
    alignas(64) long one_l[AA_NUMBER];
    alignas(64) long indexs;
    alignas(64) long total;
    alignas(64) long total_l;
    alignas(64) long complement;

public:
    std::unordered_map<long, long> vector;

    Bacteria(const std::vector<char> &data) : total(0), total_l(0), complement(0) {
        memset(one_l, 0, AA_NUMBER * sizeof(long));
        LoadData(data);
    }

    ~Bacteria() {}

    void LoadData(const std::vector<char> &data) {
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

    inline double stochastic_compute(long i) const {
        double total_plus_complement_reciprocal = 1.0 / (total + complement + EPSILON);
        double total_l_reciprocal = 1.0 / (total_l + EPSILON);
        
        auto find_second_i_div_AA_NUMBER = second.find(i / AA_NUMBER);
        double p1 = (find_second_i_div_AA_NUMBER != second.end() ? find_second_i_div_AA_NUMBER->second : 0) * total_plus_complement_reciprocal;
        
        double p2 = one_l[i % AA_NUMBER] * total_l_reciprocal;
        
        auto find_second_i_mod_M1 = second.find(i % M1);
        double p3 = (find_second_i_mod_M1 != second.end() ? find_second_i_mod_M1->second : 0) * total_plus_complement_reciprocal;
        
        double p4 = one_l[i / M1] * total_l_reciprocal;
        
        return total * (p1 * p2 + p3 * p4) / 2;
    }

private:
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

    #pragma omp parallel for
    for (long i = 0; i < number_bacteria; i++) {
        char name[10];
        #pragma omp critical
        {
            if (fscanf(input_file, "%9s", name) != 1) { // Ensure buffer is not overrun
                fprintf(stderr, "Error: failed to read bacteria name\n");
                fclose(input_file);
                exit(EXIT_FAILURE);
            }
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

double CompareBacteria(const Bacteria *b1, const Bacteria *b2) {
    double correlation = 0;
    double vector_len1 = 0;
    double vector_len2 = 0;

    #pragma omp parallel for reduction(+:correlation, vector_len1, vector_len2) schedule(static)
    for (long i = 0; i < M; i++) {
        double stochastic1 = b1->stochastic_compute(i);
        auto b1_vector_find = b1->vector.find(i);
        double b1_vector_i = (b1_vector_find != b1->vector.end()) ? b1_vector_find->second : 0;
        double t1 = (stochastic1 > EPSILON) ? (b1_vector_i - stochastic1) * (1.0 / stochastic1) : 0;
        vector_len1 += (t1 * t1);

        double stochastic2 = b2->stochastic_compute(i);
        auto b2_vector_find = b2->vector.find(i);
        double b2_vector_i = (b2_vector_find != b2->vector.end()) ? b2_vector_find->second : 0;
        double t2 = (stochastic2 > EPSILON) ? (b2_vector_i - stochastic2) * (1.0 / stochastic2) : 0;
        vector_len2 += (t2 * t2);

        correlation += t1 * t2;
    }

    vector_len1 = sqrt(vector_len1);
    vector_len2 = sqrt(vector_len2);

    correlation = correlation / (vector_len1 * vector_len2);
    return correlation;
}

void CompareAllBacteria(std::vector<std::unique_ptr<Bacteria>> &bacteria) {
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < number_bacteria - 1; i++) {
        for (int j = i + 1; j < number_bacteria; j++) {
            double correlation = CompareBacteria(bacteria[i].get(), bacteria[j].get());
            #pragma omp critical
            {
                printf("%03d %03d -> %.10lf\n", i, j, correlation);
            }
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
