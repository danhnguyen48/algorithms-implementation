#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define BINARY_TEST_FILE "results/binari_test_file.bin"
#define RESULT_FILE "results/result.txt"

// Write and read binary file
void ExitError(const char *miss, int errcode);
void read_binary_file_for_testing(unsigned long ***ways, unsigned int *amount_nodes, unsigned int **nsuccs);

int haveWay(unsigned long ***ways, unsigned int **nsuccs, unsigned long firstIndex, unsigned long secondIndex) {

    for (int i=0; i<*(*nsuccs + firstIndex); i++) {

        if (*(*(*ways + firstIndex) + i) == secondIndex)
            return 1;

    }

    return 0;

}

int test(unsigned long ***ways, unsigned int *amount_nodes, unsigned int **nsuccs) {

    FILE *fp;
    char *endptr;
    char *firstPair, *lastPair;
    char chunk[32];
    int success = 1;
    firstPair = (char *) malloc(32*sizeof(char));
    lastPair = (char *) malloc(32*sizeof(char));
    fp = fopen(RESULT_FILE, "r");
    fgets(chunk, sizeof(chunk), fp);
    memcpy(firstPair, chunk, strlen(chunk));
    while (fgets(chunk, sizeof(chunk), fp) != NULL) {

        memcpy(lastPair, chunk, strlen(chunk));
        unsigned long firstIndex, secondIndex;
        firstIndex = strtoul(firstPair, &endptr, 10);
        secondIndex = strtoul(lastPair, &endptr, 10);
        if (strlen(lastPair) == 1) continue;

        if (haveWay(ways, nsuccs, firstIndex, secondIndex) == 0) {
            printf("Cannot find the way from index %lu to %lu\n", firstIndex, secondIndex);
            success = 0;
            break;
        }

        char *temp = firstPair; firstPair = lastPair; lastPair = temp;

    }

    fclose(fp);
    return success;

}

int main(int argc, char **argv) {

    unsigned long **ways;
    unsigned int amount_nodes = 0;
    unsigned int *nsuccs;

    read_binary_file_for_testing(&ways, &amount_nodes, &nsuccs);
    if (test(&ways, &amount_nodes, &nsuccs) == 1) {
        printf("All ways are connected!\n");
    }

    return 0;

}


void ExitError(const char *miss, int errcode) {
    fprintf (stderr, "\nERROR: %s.\nStopping...\n\n", miss); exit(errcode);
}

void read_binary_file_for_testing(unsigned long ***ways, unsigned int *amount_nodes, unsigned int **nsuccs) {

    FILE *fout;
    unsigned short *nsucc = malloc(sizeof(unsigned short));
    fout = fopen(BINARY_TEST_FILE, "rb");
    if (fread(amount_nodes, sizeof(unsigned int), 1, fout) != 1)
        ExitError("Error when reading the amount of nodes", 31);
    
    *ways = (unsigned long **) malloc(*amount_nodes * sizeof(unsigned long *));
    *nsuccs = (unsigned int *) malloc(*amount_nodes * sizeof(unsigned int));
    for (int i=0; i<*amount_nodes; i++) {
        if (fread(nsucc, sizeof(unsigned short), 1, fout) != 1)
            ExitError("Error when reading the number of successor", 32);
        *(*nsuccs + i) = *nsucc;

        *(*ways + i) = (unsigned long *) malloc(*nsucc * sizeof(unsigned long));
        if (fread(*(*ways + i), sizeof(unsigned long), *nsucc, fout) != *nsucc)
            ExitError("Error when reading the ways", 32);
    }
    fclose(fout);

}

