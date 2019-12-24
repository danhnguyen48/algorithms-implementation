#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define cataluna_map "cataluna.csv"
#define spain_map "spain.csv"

typedef struct {
    unsigned long id;
    char *name;
    double lat, lon;
    unsigned short nsucc;
    unsigned long *successors;
} node;

int get_my_line(FILE *fp) {

    char chunk[128]; // buffer where is read into
    char *line = NULL;
    size_t len_used;
    while (fgets(chunk, sizeof(chunk), fp) != NULL) {
        if (line != NULL) {
            line = realloc(line, sizeof(chunk));
        } else {
            line = malloc(sizeof(chunk));
            len_used = 0;
        }
        size_t len_chunk_used = strlen(chunk);

        memcpy(line + len_used, chunk, len_chunk_used);
        len_used += len_chunk_used;

        if (*(line + len_used - 1) == '\n') {
            int len = strlen(line);
            free(line);
            return len;
        }

    }

    return -1;

}

void readFile(char *file_name, node *nodes) {

    FILE *fp = fopen(file_name, "r");

    if (fp == NULL) {
        printf("Unable to open %s file", file_name);
        exit(1);
    }

    while (get_my_line(fp) != -1) {

    }
    

}

int main(int argc, char **argv) {

    char *file_name = cataluna_map;
    node *nodes;

    if (argc>1 && strcmp(argv[1], "spain") == 0) {
        file_name = spain_map;
    }

    readFile(file_name, nodes);

    return 0;

}
