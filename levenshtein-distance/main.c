#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

void readString(char *firstString);
int min(int x, int y);
int max(int x, int y);
int levenshteinDistance(int i, int j, char * firstString, char * secondString, int * levenshteinMatrix, int * trace);
void printLevenshteinMatrix(int *levenshteinMatrix, char *firstString, char *secondString);

int main(int agrc, char** argv) {

    char *firstString, *secondString;
    int *levenshteinMatrix, *trace;
    int rows, columns;

    firstString = (char *) malloc(sizeof(char));
    secondString = (char *) malloc(sizeof(char));

    printf("The first string is: ");
    readString(firstString);
    printf("The second string is: ");
    readString(secondString);

    rows = strlen(secondString);
    columns = strlen(firstString);
    levenshteinMatrix = (int *) malloc(strlen(firstString)*strlen(secondString)*sizeof(int));

    // assumed that -1 - delete, 0 - institude, 1 - insert
    trace = (int *) malloc(strlen(firstString)*strlen(secondString)*sizeof(int));

    memset(levenshteinMatrix, -1, strlen(firstString)*strlen(secondString)*sizeof(int));
    *levenshteinMatrix = 0;

    *(levenshteinMatrix + (rows-1)*columns + columns - 1) =
                levenshteinDistance(rows - 1, columns - 1, firstString, secondString, levenshteinMatrix, trace);

    printf("The Levenshtein distance: %d\n", *(levenshteinMatrix + (rows-1)*columns + columns - 1));
    printLevenshteinMatrix(levenshteinMatrix, firstString, secondString);

}

void readString(char *firstString) {

    char c;
    int step = 1, sensor = 0;

    while (c != '\n') {
        c = getc(stdin);
        firstString = (char *) realloc(firstString, step*sizeof(char));
        *(firstString + sensor) = c;
        step++; sensor++;
    }
    *(firstString + sensor) = '\0';

}

int levenshteinDistance(int i, int j, char * firstString, char * secondString, int * levenshteinMatrix, int * trace) {

    int rows = strlen(secondString);
    int columns = strlen(firstString);

    if (min(i,j) == 0) {
        *(levenshteinMatrix + i*columns + j) = max(i,j);
        return max(i,j);
    }
    
    // left
    if (*(levenshteinMatrix + i*columns + j - 1) == -1) {
        *(levenshteinMatrix + i*columns + j - 1) = levenshteinDistance(i, j-1, firstString, secondString, levenshteinMatrix, trace);
    }
    // top
    if (*(levenshteinMatrix + (i-1)*columns + j) == -1) {
        *(levenshteinMatrix + (i-1)*columns + j) = levenshteinDistance(i-1, j, firstString, secondString, levenshteinMatrix,trace);
    }
    // topleft
    if (*(levenshteinMatrix + (i-1)*columns + j - 1) == -1) {
        *(levenshteinMatrix + (i-1)*columns + j - 1) = levenshteinDistance(i-1, j-1, firstString, secondString, levenshteinMatrix, trace);
    }

    int left = *(levenshteinMatrix + i*columns + j - 1);
    int top = *(levenshteinMatrix + (i-1)*columns + j);
    int topLeft = *(levenshteinMatrix + (i-1)*columns + j - 1);

    int minOfNeighbors = min(left, top);
    minOfNeighbors = min(minOfNeighbors, topLeft);

    if ((char)*(firstString + j - 1) == (char)*(secondString + i - 1)) 
        return minOfNeighbors;
    else
        return minOfNeighbors + 1;

}

void printLevenshteinMatrix(int *levenshteinMatrix, char *firstString, char *secondString) {

    int rows = strlen(secondString);
    int columns = strlen(firstString);
    printf("     '' ");
    for (int i=0; i<strlen(firstString); i++) 
        printf("%c ", firstString[i]);
    for (int i=0; i<rows; i++) {

        if (i == 0) printf(" ' ' ");
        else printf("   %c  ", *(secondString + i - 1));

        for (int j=0; j<columns; j++) {
            printf("%d ", *(levenshteinMatrix + i*columns + j));
        }
        printf("\n");

    }

}

int min(int x, int y) {
    return x > y ? y : x;
}

int max(int x, int y) {
    return x > y ? x : y;
}