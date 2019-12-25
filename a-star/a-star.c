#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define cataluna_map "cataluna.csv"
#define spain_map "spain.csv"
#define CHUNK_READING_SIZE 128

/*
### Format: node|@id|@name|@place|@highway|@route|@ref|@oneway|@maxspeed|node_lat|node_lon
### Format: way|@id|@name|@place|@highway|@route|@ref|@oneway|@maxspeed|member nodes|...
### Format: relation|@id|@name|@place|@highway|@route|@ref|@oneway|@maxspeed|relation_type|membertype;@id;@role|...

*/

typedef struct {
    unsigned long id;
    char *name;
    double lat, lon;
    unsigned short nsucc;
    unsigned long *successors;
} node;

int proceed_line(char *line, node **nodes, int *current) {

    if (*line == '#') return 0;

    char *found = NULL;
    int count = 0;
    char *endptr;
    found = strsep(&line, "|");
    if (strcmp(found, "node") != 0) return -1;

    while ((found = strsep(&line, "|")) != NULL) {
        count++;
        switch (count)
        {
        case 1:
            (*nodes + *current)->id = strtoul(found, &endptr, 10);
            break;
        default:
            break;
        }
    }
    (*current)++;

    return 0;

}

// protocol
int get_my_line(FILE *fp, char **line, size_t *max_len); // Gets line by line from file
void readFile(char *file_name, node **nodes, int amount_nodes); // Reads file
void count_nodes(char *file_name, int *amount); // Counts the amount of nodes


int main(int argc, char **argv) {

    char *file_name = cataluna_map;
    node *nodes = NULL;
    int amount_nodes = 0;

    if (argc>1 && strcmp(argv[1], "spain") == 0) {
        file_name = spain_map;
    }

    count_nodes(file_name, &amount_nodes);
    readFile(file_name, &nodes, amount_nodes);

    return 0;

}

int get_my_line(FILE *fp, char **line, size_t *max_len) {

    if(line == NULL || fp == NULL) {
        return -1;
    }

    char chunk[CHUNK_READING_SIZE]; // buffer where is read into

    if (*line == NULL || *max_len < sizeof(chunk)) {
        *line = malloc(sizeof(chunk));
    }
    (*line)[0] = '\0'; // empty string

    while (fgets(chunk, sizeof(chunk), fp) != NULL) {

        size_t len_chunk_used = strlen(chunk);
        size_t len_used = strlen(*line);

        if (*max_len - len_used < len_chunk_used) {
            *max_len += CHUNK_READING_SIZE;
            *line = realloc(*line, *max_len*sizeof(char));
        }

        memcpy(*line + len_used, chunk, len_chunk_used);
        len_used += len_chunk_used;
        *(*line + len_used) = '\0';

        if (*(*line + len_used - 1) == '\n') {
            return len_used;
        }

    }

    return -1;

}

void readFile(char *file_name, node **nodes, int amount_nodes) {

    FILE *fp = fopen(file_name, "r");

    if (fp == NULL) {
        printf("Unable to open %s file", file_name);
        exit(1);
    }

    char *line = NULL;
    size_t max_len = CHUNK_READING_SIZE;
    int current_node = 0;

    *nodes = (node *) malloc(amount_nodes*sizeof(node));

    while (get_my_line(fp, &line, &max_len) != -1) {

        proceed_line(line, nodes, &current_node);

    }

}


void count_nodes(char *file_name, int *amount) {

    FILE *fp = fopen(file_name, "r");
    int begin = 1;

    if (fp == NULL) {
        printf("Unable to open %s file", file_name);
        exit(1);
    }

    char chunk[CHUNK_READING_SIZE]; // buffer where is read into

    while (fgets(chunk, sizeof(chunk), fp) != NULL) {

        if (begin == 1 && chunk[0] == '#') continue;
        size_t len_chunk_used = strlen(chunk);
        if (begin == 1) {
            if (chunk[0] == 'n') (*amount)++;
            else return;
        }

        if (chunk[len_chunk_used-1] == '\n') begin = 1;
        else begin = 0;

    }

}