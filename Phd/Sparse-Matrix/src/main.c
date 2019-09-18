#include "config/config_parser.h"
#include "sparse_matrix/sparse_matrix.h"

int main (int argc, char **argv) 
{
    if (argc-1 != 1)
    {
        display_usage(argv);
        exit(EXIT_FAILURE);
    }

    struct csr_matrix *sparse_matrix = new_csr_matrix();

    read_csr_sparse_matrix_from_mtx_file(sparse_matrix,argv[1]);

    print_csr_sparse_matrix(sparse_matrix);

    csr_matrix_vector_multiplication(sparse_matrix);

/*
    struct user_options *options;
    options = new_user_options ();

    struct ode_solver *ode_solver;
    ode_solver = new_ode_solver ();

    // First we have to get the config file path
    get_config_file (argc, argv, options);

    if (options->config_file) 
    {
        // Here we parse the config file
        if (ini_parse (options->config_file, parse_config_file, options) < 0) 
        {
            fprintf (stderr, "Error: Can't load the config file %s\n", options->config_file);
            return EXIT_FAILURE;
        }
    }

    // The command line options always overwrite the config file
    parse_options (argc, argv, options);

    // Create the output dir and the logfile
    if (options->out_dir_name) 
    {
        sds buffer = sdsnew ("");
        create_dir_if_no_exists (options->out_dir_name);
        buffer = sdscatfmt (buffer, "%s/outputlog.txt", options->out_dir_name);
        open_logfile (buffer);
        sdsfree (buffer);
    }

    configure_ode_solver_from_options (ode_solver, options);
    
    solve_celular_model(ode_solver,options);

    free_ode_solver (ode_solver);
    free_user_options (options);

    close_logfile ();
*/

    return EXIT_SUCCESS;
}
