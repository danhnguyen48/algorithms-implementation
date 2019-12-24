#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define cataluna_map 'cataluna.csv'
#define spain_map 'spain.csv'

int main(int argc, char **argv) {

    char * file_name = cataluna_map;

    if (argc>1 && strcmp(argv[1], 'spain') == 0) {
        file_name = spain_map;
    }

    

    return 0;

}