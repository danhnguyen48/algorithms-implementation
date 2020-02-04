#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define cataluna_map "cataluna.csv"
#define spain_map "spain.csv"
#define CHUNK_READING_SIZE 128
#define EARTH_RADIUS 6371000
#define MAX_DOUBLE 66.9
#define OUTPUT_FILE "a-star-output.txt"

/*
### Format: node|@id|@name|@place|@highway|@route|@ref|@oneway|@maxspeed|node_lat|node_lon
### Format: way|@id|@name|@place|@highway|@route|@ref|@oneway|@maxspeed|member nodes|...
### Format: relation|@id|@name|@place|@highway|@route|@ref|@oneway|@maxspeed|relation_type|membertype;@id;@role|...
*/

struct successor {
    unsigned long index;
    struct successor *next;
};

typedef struct {
    struct successor *front, *rear;
} successor_list;

typedef struct {
    unsigned long id;
    char *name;
    double lat, lon;
    unsigned short nsucc;
    successor_list *successor_list;
    long index;
    struct QNode *queue_node;
} node;

struct QNode {
    node *key;
    double g, h;
    struct QNode *trace;
    struct QNode *next;
};

typedef struct {
    struct QNode *front, *rear;
} open;

// Trace
struct TPoint {
    unsigned int index;
    unsigned long id;
    double distance;
    char *name;
    struct TPoint *next;
};

typedef struct {
    struct TPoint *front, *rear;
} TraceQueue;

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
void create_way(unsigned long first_node, unsigned long second_node, node **nodes);
// Store data from file
int proceed_node(char *line, node **nodes, int *current);
// Store way data from file
int proceed_way(char *line, node **nodes, unsigned int amount_nodes);
// Calculate heiristic distance between two points
double heuristic_distance_between_two_points(node *first_point, node *second_point);
// Calculate real distance between two points
double equirectangular_approximation(node *first_point, node *second_point);
// Create open list
open *create_queue();
// Import qnode into queue following sorted f
int import_queue(open *q, struct QNode *qnode);
// Dequeue queue
struct QNode *de_queue(open *q);
// New QNode
struct QNode *new_qnode(double g, double h, node *key);
// Check if queue is empty => 1 is empty and 0 is no empty
int queue_empty(open *q);
// Convert degree to radians
double convert_radians(double degree);
// Check if node is in open list
struct QNode *is_node_in_list(open *q, node *n);
// Create successor linked list
successor_list *create_successor_list();
// Import successor list
int import_successor_list(successor_list *list, struct successor *s);
// New successor
struct successor *new_successor(unsigned long neighbor_index);
// A-star
void a_star(node *source, node *goal, open *q, unsigned int amount_nodes, node **nodes, char *output);
// Track ways
void trace_back(struct QNode *queue_goal, node **nodes, unsigned long source_id, char *output);

int main(int argc, char **argv) {

    char *file_name = cataluna_map;
    char *output_file = NULL;
    node *nodes = NULL;
    unsigned int amount_nodes = 0;
    open *open_list = create_queue();
    unsigned long source_id = 771979683;
    unsigned long goal_id = 429854583;
    clock_t start, end;
    double cpu_time_used;

    start = clock();

    if (argc>1 && strcmp(argv[1], "spain") == 0) {
        file_name = spain_map;
        source_id = 240949599;
        goal_id = 195977239;
    }

    if (argc>2 && strcmp(argv[2], "file") == 0) {
        output_file = OUTPUT_FILE;
    }

    count_nodes(file_name, &amount_nodes);

    readFile(file_name, &nodes, amount_nodes);

    long barcelona_index = binary_search_node(source_id, nodes, amount_nodes);
    long sevilla_index = binary_search_node(goal_id, nodes, amount_nodes);
    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Reading Time: %f\n", cpu_time_used);

    printf("barcelona: %ld, sevilla: %ld\n", barcelona_index, sevilla_index);

    a_star(&nodes[barcelona_index], &nodes[sevilla_index], open_list, amount_nodes, &nodes, output_file);

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("All Time: %f\n", cpu_time_used);

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
            (*nodes + *current)->successor_list = create_successor_list();
            (*nodes + *current)->queue_node = NULL;
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

    while ((found = strsep(&line, "|")) != NULL) {
        count++;
        switch (count)
        {
        case 1:case 2:case 3: case 4: case 5: case 6: case 8:
            break;
        case 7:
            if (strcmp(found, "oneway") == 0) one_way = 1;
            else one_way = 0;
            break;
        default:
            node_id = strtoul(found, &endptr, 10);
            if (first_node == -1) {
                first_node = binary_search_node(node_id, *nodes, amount_nodes);
                continue;
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

void create_way(unsigned long first_node, unsigned long second_node, node **nodes) {

    // Expand one memory block
    struct successor *s = new_successor(second_node);
    import_successor_list((*nodes + first_node)->successor_list, s);

}

double heuristic_distance_between_two_points(node *first_point, node *second_point) {

    double delta_lat = convert_radians(first_point->lat - second_point->lat);
    double delta_lon = convert_radians(first_point->lon - second_point->lon);
    double first_lat = convert_radians(first_point->lat);
    double second_lat = convert_radians(second_point->lat);

    double a = sin(delta_lat/2)*sin(delta_lat/2) + cos(first_lat)*cos(second_lat)*sin(delta_lon/2)*sin(delta_lon/2);
    double c = 2*atan2(sqrt(a), sqrt(1-a));
    double d = EARTH_RADIUS * c; 

    return d;

}

double equirectangular_approximation(node *first_point, node *second_point) {

    double lat1 = convert_radians(first_point->lat);
    double lat2 = convert_radians(second_point->lat);
    double delta_lon = convert_radians(first_point->lon - second_point->lon);
    double x = delta_lon * cos((lat1 + lat2)/2);
    double y = lat1 - lat2;
    double d = sqrt(x*x + y*y) * EARTH_RADIUS;

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

    double qnode_f = qnode->g + qnode->h;
    if (qnode_f <= q->front->g + q->front->h) {
        qnode->next = q->front;
        q->front = qnode;
        return 1;
    }

    struct QNode *temp = q->front;

    while (temp->next != NULL && qnode_f > temp->next->g + temp->next->h)
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
    // free(temp);
    q->front = q->front->next;
    temp->next = NULL;

    // If front NULL, queue will be NULL so change rear to be NULL
    if (q->front == NULL)
        q->rear = NULL;

    return temp;

}

struct QNode *new_qnode(double g, double h, node *key) {

    struct QNode *qnode = (struct QNode *) malloc(sizeof(struct QNode));
    qnode->g = g;
    qnode->h = h;
    qnode->key = key;
    qnode->next = NULL;
    qnode->trace = NULL;
    key->queue_node = qnode;
    
    return qnode;

}

int queue_empty(open *q) {
    if (q->front == NULL) return 1;
    return 0;
}

struct QNode *is_node_in_list(open *q, node *n) {

    if (q->rear == NULL || q->front == NULL || n->queue_node == NULL) return NULL;

    struct QNode *temp = q->front;
    if (temp == n->queue_node) {
        q->front = temp->next;
        temp->next = NULL;
        return temp;
    }

    while (temp->next != NULL) {
        if (temp->next == n->queue_node) {
            struct QNode *qnode = temp->next; // store pointer which is needed to return
            temp->next = qnode->next; // remove n from queue
            qnode->next = NULL;
            if (temp->next == NULL) 
                q->rear = temp; // If the qnode after the needed node is null, the before qnode will be the rear
            return qnode;
        }
        temp = temp->next;
    }

    return NULL;

}

void a_star(node *source, node *goal, open *q, unsigned int amount_nodes, node **nodes, char *output) {

    struct QNode *current;
    double estimated_g;
    double distance;
    long index_neighbor;
    double neighbor_h;

    // Create source in queue
    double source_h = heuristic_distance_between_two_points(source, goal);
    double source_f = source_h;
    struct QNode *source_qnode = new_qnode(0.0, source_h, source);
    import_queue(q, source_qnode);

    while (queue_empty(q) == 0) {

        current = de_queue(q);

        if (current->key->id == goal->id) {
            printf("The lowest cost is %f\n", current->g);
            trace_back(current, nodes, source->id, output);
            return;
        }
        // Generate each state node_successor that come after node_current
        struct successor *succ = current->key->successor_list->front;
        while (succ != NULL) {

            index_neighbor = succ->index;
            distance = equirectangular_approximation(current->key, *nodes + index_neighbor);
            estimated_g = current->g + distance;

            if ((*nodes + index_neighbor)->queue_node == NULL || estimated_g < (*nodes + index_neighbor)->queue_node->g) { // This path is better than previous one
                
                // Check neighbor is in openlist
                // if yes, update neighbor f
                // if no, insert into queue
                struct QNode *neighbor_qnode_in_q = is_node_in_list(q, *nodes + index_neighbor); // Check and dequeue this qnode
                if (neighbor_qnode_in_q != NULL) {
                    neighbor_qnode_in_q->trace = current;
                    neighbor_qnode_in_q->g = estimated_g;
                    // Update position
                    import_queue(q, neighbor_qnode_in_q);
                } else {
                    neighbor_h = heuristic_distance_between_two_points(*nodes + index_neighbor, goal);
                    struct QNode *neighbor_qnode = new_qnode(estimated_g, neighbor_h, *nodes + index_neighbor);
                    neighbor_qnode->trace = current;
                    import_queue(q, neighbor_qnode);
                }

            }

            succ = succ->next;

        }

    }

}

successor_list *create_successor_list() {
    successor_list *sl = (successor_list *) malloc(sizeof(successor_list));
    sl->front = sl->rear = NULL;
    return sl;
}

int import_successor_list(successor_list *list, struct successor *s) {

    if (list->front == NULL || list->rear == NULL) {
        list->front = s;
        list->rear = s;
        return 1;
    }

    list->rear->next = s;
    list->rear = s;

    return 1;

}

struct successor *new_successor(unsigned long neighbor_index) {

    struct successor *s = (struct successor *) malloc(sizeof(struct successor));
    s->index = neighbor_index;
    s->next = NULL;
    
    return s;

}

void trace_back(struct QNode *queue_goal, node **nodes, unsigned long source_id, char *output) {

    TraceQueue *trace_queue = (TraceQueue *) malloc(sizeof(TraceQueue));
    trace_queue->front = trace_queue->rear = NULL;
    struct TPoint *tpoint = (struct TPoint *) malloc(sizeof(struct TPoint));
    tpoint->id = queue_goal->key->id;
    tpoint->index = queue_goal->key->index;
    tpoint->name = (char *) malloc(sizeof(queue_goal->key->name));
    tpoint->name = queue_goal->key->name;
    tpoint->distance = queue_goal->g;
    tpoint->next = NULL;
    trace_queue->front = tpoint;
    trace_queue->rear = tpoint;

    struct QNode *temp = queue_goal;
    while (temp->key->id != source_id) {
        temp = temp->trace;
        struct TPoint *tp = (struct TPoint *) malloc(sizeof(struct TPoint));
        tp->id = temp->key->id;
        tp->index = temp->key->index;
        tp->name = (char *) malloc(sizeof(temp->key->name));
        tp->name = temp->key->name;
        tp->distance = temp->g;
        tp->next = trace_queue->front;
        trace_queue->front = tp;
    }

    struct TPoint *tp = (struct TPoint *) malloc(sizeof(struct TPoint));
    tp = trace_queue->front;

    FILE *fp = NULL;
    if (output != NULL) {
        fp = fopen(output, "w");
    }

    do {
        if (fp == NULL)
            printf("Node id: %lu\t| Distance: %f  \t| Name: %s\n", tp->id, tp->distance, tp->name);
        else
            fprintf(fp, "Node id: %lu\t| Distance: %f  \t| Name: %s\n", tp->id, tp->distance, tp->name);
            // fprintf(fp, "Index: %u | Node id: %lu | Distance: %f | Name: %s\n", tp->index, tp->id, tp->distance, tp->name);
        tp = tp->next;
    } while (tp->id != queue_goal->key->id);
    if (fp == NULL)
        printf("Node id: %lu\t| Distance: %f  \t| Name: %s\n", tp->id, tp->distance, tp->name);
    else
        fprintf(fp, "Node id: %lu\t| Distance: %f  \t| Name: %s\n", tp->id, tp->distance, tp->name);
        // fprintf(fp, "Index: %u | Node id: %lu | Distance: %f | Name: %s\n", tp->index, tp->id, tp->distance, tp->name);

    fclose(fp);

}
