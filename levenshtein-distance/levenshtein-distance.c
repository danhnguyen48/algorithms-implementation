#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

void readString(char *firstString);
int min(int x, int y);
int max(int i, int j, int * levenshteinMatrix, int * trace, int columns);
int levenshteinDistance(int i, int j, char * firstString, char * secondString, int * levenshteinMatrix, int * trace);
int minOfNeighbors(int left, int top, int topleft, int * traceElement);
void printLevenshteinMatrix(int *levenshteinMatrix, int *trace, char *firstString, char *secondString);
void traceTransformation(char * firstString, char * secondString, int * trace);
void insertCharacter(char *string, char character, int pos);
void substituteCharacter(char *string, char character, int pos);
void delete(char *string, int pos);

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
    *levenshteinMatrix = 0; // First element should be 0

    *(levenshteinMatrix + (rows-1)*columns + columns - 1) =
                levenshteinDistance(rows - 1, columns - 1, firstString, secondString, levenshteinMatrix, trace);

    printf("The Levenshtein distance: %d\n", *(levenshteinMatrix + (rows-1)*columns + columns - 1));
    printLevenshteinMatrix(levenshteinMatrix, trace, firstString, secondString);
    traceTransformation(firstString, secondString, trace);
    
    free(firstString); free(secondString);
    free(levenshteinMatrix); free(trace);
    
    return 0;

}

// Read string char by char to avoid memory waste
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
        return max(i, j, levenshteinMatrix, trace, columns);
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

    int minOfNeighborsElement = minOfNeighbors(left, top, topLeft, trace + i*columns + j);

    if ((char)*(firstString + j - 1) == (char)*(secondString + i - 1)) 
        return minOfNeighborsElement;
    else
        return minOfNeighborsElement + 1;

}

void printLevenshteinMatrix(int *levenshteinMatrix, int *trace, char *firstString, char *secondString) {

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

void traceTransformation(char * firstString, char * secondString, int * trace) {

    char *temp = (char *) malloc(strlen(firstString)*sizeof(char));
    int j = strlen(firstString) - 1, i = strlen(secondString) - 1;
    int columns = strlen(firstString), rows = strlen(secondString);

    strcpy(temp, firstString);
    printf("First String: %s\n", temp);

    while (i != 0 && j != 0) {
        int pos = j;
        switch (*(trace + i*columns + j)) {
            case 1: // Insert
                printf("Insert after character %c at position %d\n", (char) *(secondString + i - 1), pos - 1);
                insertCharacter(temp, (char) *(secondString + i - 1), pos);
                printf("New String: %s\n", temp);
                i--;
                break;
            case -1: // Delete
                printf("Delete %c at position %d\n", *(temp + pos - 1), pos - 1);
                delete(temp, pos);
                printf("New String: %s\n", temp);
                j--;
                break;
            case 0:
                printf("Substitute %c at position %d by %c\n", *(temp + pos - 1), pos - 1, (char) *(secondString + i - 1));
                substituteCharacter(temp, (char) *(secondString + i - 1), pos);
                printf("New String: %s\n", temp);
                i--; j--;
                break;
            default: 
                break;
        }
    }

    printf("Final string: %s\n", temp);

}

int minOfNeighbors(int left, int top, int topleft, int * traceElement) {
    int min = 999999;
    if (min >= topleft) {
        min = topleft;
        *traceElement = 0;
    }
    if (min >= top) {
        min = top;
        *traceElement = 1;
    }
    if (min >= left) {
        min = left;
        *traceElement = -1;
    }
    return min;
}

int min(int x, int y) {
    return x > y ? y : x;
}

int max(int i, int j, int * levenshteinMatrix, int * trace, int columns) {
    if (i == 0) {
        *(trace + j) = -1;
        *(levenshteinMatrix + j) = j;
    } else {
        *(trace + i*columns) = 1;
        *(levenshteinMatrix + i*columns) = i;
    }
    return *(levenshteinMatrix + i*columns + j);
}

void insertCharacter(char *string, char character, int pos) {
    char *buf = (char *) malloc((pos)*sizeof(char));
    strncpy(buf, string, pos);
    buf = (char *) realloc(buf, sizeof(char));
    *(buf + pos) = character;
    buf = (char *) realloc(buf, (strlen(string) - pos)*sizeof(char));
    strcpy(buf + pos+1, string + pos);
    string = (char *) realloc(string, sizeof(char));
    strcpy(string, buf);
}

void substituteCharacter(char *string, char character, int pos) {
    char *pointPos = string + pos - 1;
    *pointPos = character;
}

void delete(char *string, int pos) {
    memmove(string + pos - 1, string + pos, (strlen(string) - 1)*sizeof(char));
}
