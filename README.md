# Bioinformatics - Genome Similarity Application

## Overview

### What is it?
The Genome Similarity application is a specialized tool in C++ for analyzing amino acid sequences in genomes, using frequency vectors to compare genomic data. It's essential for research in biology, disease, and viral studies.

### What does it do?
The application calculates a correlation score between -1 and 1 for amino acids in protein sequences, using k-mers. It provides comprehensive comparisons crucial for various bioinformatics analytics.

### How does it work?
The application initializes by reading a text file of genome names in FASTA format. Each sequence is processed into a frequency vector, which is then used to calculate correlation scores using a stochastic approach. The program is designed to optimize memory usage and efficiency.

## Architecture

1. **Initialization:** Sets up data structures, encoding for amino acids, and initializes the Bacteria class for frequency vector calculations.
2. **Input Reading:** Reads from 'list.txt', listing genome files in .faa format.
3. **Genome Processing:** Processes each genome file, creating frequency vectors from amino acid sequences.
4. **Correlation Calculation:** Compares frequency vectors of genomes to calculate correlation scores.
5. **Output:** Outputs correlation scores and performance metrics.

## Files and Modules

- `tidy.cpp`: Main application file.
- `list.txt`: Input file listing genome files to process.
- `data` folder: Contains .faa genome files.

## Functions

1. `Init()`: Initializes global variables for array sizes.
2. `Bacteria Class`: Represents a bacteria species, processes genomic data.
   - `Constructor: Bacteria(char* filename)`
   - `double stochastic_compute(long i)`
3. `ReadInputFile(const char* input_name)`: Reads genome filenames.
4. `double CompareBacteria(Bacteria* b1, Bacteria* b2)`: Calculates correlation between two bacteria species.
5. `void CompareAllBacteria()`: Compares all bacteria pairs listed in the input file.

## Performance Optimization and Parallelization

### Key Optimizations:
- **CompareAllBacteria Function:** Identified as a major bottleneck, optimized using parallelization.
- **CompareBacteria:** Significant time spent in this function, optimized by reducing computational overhead and parallelizing loops.
- **Stochastic_compute:** Optimized by simplifying mathematical calculations.

### Parallelization Approach:
- Utilized OpenMP on an 8-core M2 Mac for efficient task distribution and synchronization.
- Specific focus on parallelizing genome comparisons and computational functions.

### Memory Management:
- Code optimized for on-the-fly file reading to save memory, with adjustments for better performance on hardware with low memory.

## Performance Evaluation
- Significant reduction in runtime, from 858 seconds to 98 seconds.

## Challenges and Solutions

- Addressed various challenges such as I/O strategy, memory management issues, and optimization balance.
- Solutions included adjustments in data reading strategies and memory management using std::vector and smart pointers.

## Reflection and Conclusions

- The optimization process highlighted the importance of balancing memory and hardware optimization.
- The final optimized code achieved significant speed improvements, demonstrating the effectiveness of the optimization strategies.

## Instructions

This repository includes three versions of a bacteria sequence comparison program:

- `original.cpp`: The original, unoptimized version.
- `optimised.cpp`: An optimized version for better performance.
- `ultra-optimised.cpp`: Highly optimized with parallel processing and memory alignment. May slow down your computer.

## Compilation Instructions

### Windows/macOS

- Compile using g++ with OpenMP support.
- Minimum 8GB of RAM recommended.
- Command: `g++-11 -fopenmp -Ofast <filename>.cpp -o <outputname>`
- Set OpenMP threads: 
  - Windows: `set OMP_NUM_THREADS=8`
  - macOS: `export OMP_NUM_THREADS=8`

### macOS Specific

- Requires LLVM with OpenMP (install via Homebrew).
- Command: `/opt/homebrew/opt/llvm/bin/clang++ -Xpreprocessor -fopenmp -O3 -march=native -funroll-loops <filename>.cpp -o <outputname> -I/opt/homebrew/opt/llvm/include -L/opt/homebrew/opt/llvm/lib -lomp`
- Set OpenMP threads: `export OMP_NUM_THREADS=8`

## Required Header Files and Libraries

Includes standard C/C++ libraries such as `stdio.h`, `stdlib.h`, `string.h`, `time.h`, `math.h`, `vector`, `omp.h`, and `memory`.

## Running the Programs

- Execute the compiled file (e.g., `./optimised`).
- Adjust the number of threads according to your CPU capabilities for best performance.

## Important Note

- Not tested on Windows; developed on macOS.
- The C++ code is basic and should be portable. For any issues, please contact.


