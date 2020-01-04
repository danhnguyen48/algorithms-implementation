#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define BINARY_TEST_FILE "results/binari_test_file.bin"
#define RESULT_FILE "results/result.txt"

// Write and read binary file
void ExitError(const char *miss, int errcode);
void read_binary_file_for_testing(unsigned int ***ways, unsigned int *amount_nodes, unsigned short **nsuccs, char *binary_file);

int haveWay(unsigned int ***ways, unsigned short **nsuccs, unsigned long firstIndex, unsigned long secondIndex) {

    for (int i=0; i<*(*nsuccs + firstIndex); i++) {

        if (*(*(*ways + firstIndex) + i) == secondIndex)
            return 1;

    }

    return 0;

}

int test(unsigned int ***ways, unsigned int *amount_nodes, unsigned short **nsuccs, char *file_name) {

    FILE *fp;
    char *endptr;
    char *firstPair, *lastPair;
    char chunk[12];
    int success = 1;
    firstPair = (char *) malloc(12*sizeof(char));
    lastPair = (char *) malloc(12*sizeof(char));
    if (file_name == NULL) {
        printf("Uploaded file not found!\n");
        exit(0);
    } else
    if ((fp = fopen(file_name, "r")) == NULL) {
        printf("Cannot read the uploaded file!\n");
        exit(0);
    }
    fgets(chunk, sizeof(chunk), fp);
    memcpy(firstPair, chunk, strlen(chunk));
    while (fgets(chunk, sizeof(chunk), fp) != NULL) {

        memcpy(lastPair, chunk, strlen(chunk));
        unsigned long firstIndex, secondIndex;
        firstIndex = strtoul(firstPair, &endptr, 10);
        secondIndex = strtoul(lastPair, &endptr, 10);
        if (strlen(lastPair) == 1) continue;

        if (haveWay(ways, nsuccs, firstIndex, secondIndex) == 0) {
            printf("Do not find the way from index %lu to %lu\n", firstIndex, secondIndex);
            success = 0;
            break;
        }

        char *temp = firstPair; firstPair = lastPair; lastPair = temp;

    }

    fclose(fp);
    return success;

}

int main(int argc, char **argv) {

    unsigned int **ways;
    unsigned int amount_nodes = 0;
    unsigned short *nsuccs;
    char *file_name;
    char *binary_file;

    if (argc < 2) {
        printf("Parameter missing!\n");
        exit(0);
    } else {
        file_name = argv[1];
        binary_file = argv[2];
    }

    read_binary_file_for_testing(&ways, &amount_nodes, &nsuccs, binary_file);
    if (test(&ways, &amount_nodes, &nsuccs, file_name) == 1) {
        printf("All ways are connected!\n");
    }

    return 0;

}


void ExitError(const char *miss, int errcode) {
    fprintf (stderr, "\nERROR: %s.\nStopping...\n\n", miss); exit(errcode);
}

void read_binary_file_for_testing(unsigned int ***ways, unsigned int *amount_nodes, unsigned short **nsuccs, char *binary_file) {

    FILE *fout;
    unsigned short *nsucc = malloc(sizeof(unsigned short));
    if (binary_file == NULL || (fout = fopen(binary_file, "rb")) == NULL) {
        printf("Cannot read source file from system!\n");
        exit(0);
    }
    if (fread(amount_nodes, sizeof(unsigned int), 1, fout) != 1)
        ExitError("Error when reading the amount of nodes", 31);
    
    *ways = (unsigned int **) malloc(*amount_nodes * sizeof(unsigned int *));
    *nsuccs = (unsigned short *) malloc(*amount_nodes * sizeof(unsigned short));
    for (int i=0; i<*amount_nodes; i++) {
        if (fread(nsucc, sizeof(unsigned short), 1, fout) != 1)
            ExitError("Error when reading the number of successor", 32);
        *(*nsuccs + i) = *nsucc;

        *(*ways + i) = (unsigned int *) malloc(*nsucc * sizeof(unsigned int));
        if (fread(*(*ways + i), sizeof(unsigned int), *nsucc, fout) != *nsucc)
            ExitError("Error when reading the ways", 32);
    }
    fclose(fout);

}

