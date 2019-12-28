#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define cataluna_map "cataluna.csv"
#define spain_map "spain.csv"
#define CHUNK_READING_SIZE 128
#define EARTH_RADIUS 6371000
#define MAX_DOUBLE 66.9

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
    long index;
} node;

struct QNode {
    node *key;
    double f;
    double h;
    struct QNode *next;
};

typedef struct {
    struct QNode *front, *rear;
} open;

// PROTOCOL

// Gets line by line from file
int get_my_line(FILE *fp, char **line, size_t *max_len);
// Reads file
void readFile(char *file_name, node **nodes, unsigned int amount_nodes);
// Counts the amount of nodes
void count_nodes(char *file_name, unsigned int *amount);
// Search node by id
long binary_search_node(unsigned long id, node *nodes, int len_of_node);
// Create a way from first to second
void create_way(int first_node, int second_node, node **nodes);
// Store data from file
int proceed_node(char *line, node **nodes, int *current);
// Store way data from file
int proceed_way(char *line, node **nodes, unsigned int amount_nodes);
// Calculate distance between two points
double distance_between_two_points(node *first_point, node *second_point);
// Create open list
open *create_queue();
// Import qnode into queue
int import_queue(open *q, struct QNode *qnode);
// Dequeue queue
struct QNode *de_queue(open *q);
// New QNode
struct QNode *new_qnode(double f, double h, node *key);
// Check if queue is empty => 1 is empty and 0 is no empty
int queue_empty(open *q);
// Convert degree to radians
double convert_radians(double degree);
// Check if node is in open list
struct QNode *is_node_in_open_list(open *q, node n);
// A-star
void a_start(node *source, node *goal, open *q, int **trace, double **g, unsigned int amount_nodes, node **nodes);


int main(int argc, char **argv) {

    char *file_name = cataluna_map;
    node *nodes = NULL;
    unsigned int amount_nodes = 0;
    unsigned short *nsuccdim;
    int *trace = NULL;
    double *g = NULL;
    open *open_list = create_queue();

    if (argc>1 && strcmp(argv[1], "spain") == 0) {
        file_name = spain_map;
    }

    count_nodes(file_name, &amount_nodes);

    // allocate trace and g into memory
    trace = (int *) malloc(amount_nodes*sizeof(int));
    g = (double *) malloc(amount_nodes*sizeof(double));
    // initilize data for trace and g
    memset(trace, -1, amount_nodes*sizeof(int));
    memset(g, MAX_DOUBLE, amount_nodes*sizeof(double));

    readFile(file_name, &nodes, amount_nodes);

    long barcelona_index = binary_search_node(240949599, nodes, amount_nodes);
    long sevilla_index = binary_search_node(195977239, nodes, amount_nodes);

    printf("barcelona: %ld, sevilla: %ld\n", barcelona_index, sevilla_index);

    a_start(&nodes[barcelona_index], &nodes[sevilla_index], open_list, &trace, &g, amount_nodes, &nodes);

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

void readFile(char *file_name, node **nodes, unsigned int amount_nodes) {

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
    fclose(fp);

}


void count_nodes(char *file_name, unsigned int *amount) {

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
            (*nodes + *current)->index = *current;
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

int proceed_way(char *line, node **nodes, unsigned int amount_nodes) {
    
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
                first_node = binary_search_node(node_id, *nodes, amount_nodes);
            } else {
                if ((second_node = binary_search_node(node_id, *nodes, amount_nodes)) == -1) {
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

long binary_search_node(unsigned long id, node *nodes, int len_of_node) {
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

double convert_radians(double degree) {
    return degree*M_PI/180;
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

double distance_between_two_points(node *first_point, node *second_point) {

    double delta_lat = convert_radians(first_point->lat - second_point->lat);
    double delta_lon = convert_radians(first_point->lon - second_point->lon);
    double first_lat = convert_radians(first_point->lat);
    double second_lat = convert_radians(second_point->lat);

    double a = sinf(delta_lat/2)*sinf(delta_lat/2) + cosf(first_lat)*cosf(second_lat)*sinf(delta_lon/2)*sinf(delta_lon/2);
    double c = 2*atan2f(sqrtf(a), sqrtf(1-a));
    double d = EARTH_RADIUS * c;

    return d;

}


open* create_queue() {

    open *q = (open *) malloc(sizeof(open));
    q->front = q->rear = NULL;
    return q;

}

int import_queue(open *q, struct QNode *qnode) {

    if (q->rear == NULL) {
        q->front = q->rear = qnode;
        return 1;
    }

    if (qnode->f <= q->front->f) {
        qnode->next = q->front;
        q->front = qnode;
        return 1;
    }

    struct QNode *temp = q->front;

    while (temp->next != NULL && qnode->f > temp->next->f)
        temp = temp->next;
    
    if (temp->next == NULL) {
        q->rear->next = qnode;
        q->rear = qnode;
    } else {
        qnode->next = temp->next;
        temp->next = qnode;
    }
    return 1;

}

struct QNode *de_queue(open *q) {

    // Queue empty
    if (q->front == NULL) return NULL;

    struct QNode *temp = q->front;
    free(temp);
    q->front = q->front->next;

    // If front NULL, queue will be NULL so change rear to be NULL
    if (q->front == NULL)
        q->rear = NULL;

    return temp;

}

struct QNode *new_qnode(double f, double h, node *key) {

    struct QNode *qnode = (struct QNode *) malloc(sizeof(struct QNode));
    qnode->f = f;
    qnode->h = h;
    qnode->key = key;
    qnode->next = NULL;
    
    return qnode;

}

int queue_empty(open *q) {
    if (q->front == NULL) return 1;
    return 0;
}

struct QNode *is_node_in_open_list(open *q, node n) {

    if (q->rear == NULL || q->front == NULL) return 0;

    struct QNode *temp = q->front;

    while (temp->key->id != n.id) {
        temp = temp->next;
        if (temp == NULL) return NULL;
    }
    return temp;

}


void a_start(node *source, node *goal, open *q, int **trace, double **g, unsigned int amount_nodes, node **nodes) {

    struct QNode *current;
    double estimated_g;
    double distance;
    long index_neighbor;
    double neighbor_h;

    // Create source in queue
    double source_h = distance_between_two_points(source, goal);
    double source_f = source_h;
    struct QNode *source_qnode = new_qnode(source_f, source_h, source);
    *(*g + source->index) = 0;
    import_queue(q, source_qnode);

    while (queue_empty(q) == 0) {

        current = de_queue(q);

        if (current->key->id == goal->id) { // Check lai cho nay
            printf("The lowest cost is %f\n", *(*g + current->key->index));
            return;
        }

        for (int i=0; i<current->key->nsucc; i++) {

            index_neighbor = current->key->successors[i];
            distance = distance_between_two_points(current->key, *nodes + index_neighbor);
            estimated_g = *(*g + current->key->index) + distance;

            if (estimated_g < *(*g + index_neighbor)) { // This path is better than previous one
                
                *(*trace + index_neighbor) = current->key->index;
                *(*g + index_neighbor) = estimated_g;
                
                // Check neighbor is in openlist
                // if yes, update neighbor f
                // if no, insert into queue
                neighbor_h = distance_between_two_points(*nodes + index_neighbor, goal);
                struct QNode *neighbor_qnode_in_q = is_node_in_open_list(q, *(*nodes + index_neighbor));
                if (neighbor_qnode_in_q != NULL) {
                    neighbor_qnode_in_q->f = *(*g + index_neighbor) + neighbor_h;
                } else {
                    double neighbor_f = *(*g + index_neighbor) + neighbor_h;
                    struct QNode *neighbor_qnode = new_qnode(neighbor_f, neighbor_h, *nodes + index_neighbor);
                    import_queue(q, neighbor_qnode);
                }

            }

        }

    }

}