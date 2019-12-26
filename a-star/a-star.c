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

// PROTOCOL

// Gets line by line from file
int get_my_line(FILE *fp, char **line, size_t *max_len);
// Reads file
void readFile(char *file_name, node **nodes, int amount_nodes);
// Counts the amount of nodes
void count_nodes(char *file_name, int *amount);
// Search node by id
long binari_search_node(unsigned long id, node *nodes, int len_of_node);
// Create a way from first to second
void create_way(int first_node, int second_node, node **nodes);
// Store data from file
int proceed_node(char *line, node **nodes, int *current);
// Store way data from file
int proceed_way(char *line, node **nodes, int amount_nodes);


int proceed_way(char *line, node **nodes, int amount_nodes) {
    
    char *found = NULL;
    found = strsep(&line, "|");
    if (strcmp(found, "way") != 0) return -1;
    int count = 0;
    char *endptr;
    int one_way = 1;
    long first_node = -1, second_node = -1;
    unsigned long node_id;
    unsigned long tempId = 0;

    while ((found = strsep(&line, "|")) != NULL) {
        count++;
        switch (count)
        {
        case 1:
            tempId = strtoul(found, &endptr, 10);
            break;
        case 2:case 3: case 4: case 5: case 6: case 8:
            break;
        case 7:
            if (strcmp(found, "oneway") == 0) one_way = 1;
            else one_way = 0;
            break;
        default:
            node_id = strtoul(found, &endptr, 10);
            if (first_node == -1) {
                first_node = binari_search_node(node_id, *nodes, amount_nodes);
            } else {
                if ((second_node = binari_search_node(node_id, *nodes, amount_nodes)) == -1) {
                    first_node = -1;
                    continue;
                }
            }
            if (first_node == -1 || second_node == -1) continue;

            // Create way from first_node to second_node
            create_way(first_node, second_node, nodes);
            // Create way from second_node to first_node if !oneway
            if (one_way == 0) {
                create_way(second_node, first_node, nodes);
            }
            first_node = second_node;

            break;
        }
    }

    return 0;

}



int main(int argc, char **argv) {

    char *file_name = cataluna_map;
    node *nodes = NULL;
    int amount_nodes = 0;
    unsigned short *nsuccdim;

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

        if (*line == '#') continue;

        if (*line == 'n')
            proceed_node(line, nodes, &current_node);
        else if (*line == 'w') 
            proceed_way(line, nodes, amount_nodes);
        else break;

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

int proceed_node(char *line, node **nodes, int *current) {

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
            (*nodes + *current)->nsucc = 0;
            break;
        case 2:
            (*nodes + *current)->name = (char *) malloc(strlen(found)*sizeof(char));
            strcpy((*nodes + *current)->name, found);
            break;
        case 9:
            (*nodes + *current)->lat = strtod(found, &endptr);
            break;
        case 10:
            (*nodes + *current)->lon = strtod(found, &endptr);
            break;
        default:
            break;
        }
    }
    (*current)++;

    return 0;

}

long binari_search_node(unsigned long id, node *nodes, int len_of_node) {
    int imin = 0;
    int imax = len_of_node - 1;
    while (imax >= imin) {
        int imid = (imin + imax) / 2;
        if (nodes[imid].id == id) return imid;
        else if (nodes[imid].id < id) imin = imid + 1;
        else imax = imid - 1;
    }
    return -1;
}

void create_way(int first_node, int second_node, node **nodes) {

    // Expand one memory block
    if ((*nodes + first_node)->nsucc == 0) {
        (*nodes + first_node)->successors = (unsigned long *) malloc(sizeof(unsigned long));
    } else {
        (*nodes + first_node)->successors = realloc((*nodes + first_node)->successors, sizeof(unsigned long));
    }
    // Create way from first_node to second_node
    (*nodes + first_node)->successors[(*nodes + first_node)->nsucc] = second_node;
    (*nodes + first_node)->nsucc++;

}

