#include <stdio.h>
#include "../include/structures.h"
#include "../include/parameters.h"
#include "../include/logs.h"

int main(int argc, char *argv[]) {
    // Get parameters input file path
    char *parameters_file_path;
    if (argc == 2) {
        parameters_file_path = argv[1];
    } else {
        printf("Usage: %s <parameters_file_path>\n", argv[0]);
        return 1;
    }
    Parameters parameters = read_parameters_file(parameters_file_path);
    log_parameters(&parameters, 0);
    log_parameters(&parameters, 1);
    return 0;
}
