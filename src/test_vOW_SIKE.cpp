
#include <string>
#include <cmath>
#define VOW_SIKE

extern "C"
{
#include "sidh_vow_base.h"
#include "sike_vow_instances.h"
#include "state.h"
#include "bintree.h"
#include "test_extras.h"
}
#include "benchmark_concurrent_ds.h"

#define ATTACK_NAME "vOW_SIKE"
#define BENCH_LOOPS 1

// greek letters
#if (OS_TARGET == OS_WIN)
#define _ALPHA_CHAR "%c", 224
#define _BETA_CHAR "%c", 225
#define _GAMMA_CHAR "%c", 226
#define _DELTA_CHAR "%c", 235
#elif (OS_TARGET == OS_LINUX)
#include <locale.h>
#define _ALPHA_CHAR ("α")
#define _BETA_CHAR ("β")
#define _GAMMA_CHAR ("γ")
#define _DELTA_CHAR ("δ")
#endif

int arg_delta = 0; // delta passed from CLI
int arg_w = 0;     // w passed from CLI

int stats_vow(uint64_t n_cores, bool hansel_gretel, bool collect_stats, uint64_t max_crumbs)
{
    unsigned int i, j;
    uint64_t random_functions, collisions, mem_collisions, dist_points, number_steps_collect, number_steps_locate, number_steps, dist_cols;
    bool success;
    double calendar_time, total_time;
    unsigned long long cycles, cycles1, cycles2;
    shared_state_t S;
#if (OS_TARGET == OS_LINUX)
    // on linux, set utf8 support
    setlocale(LC_ALL, "");
#endif

    printf("\nRunning vOW attack on SIKE");
    printf("\n----------------------------------------------------------------------------------------\n\n");

    for (i = 0; i < NUM_INSTS_STATS; i++)
    {
        success = true;
        random_functions = 0;
        collisions = 0;
        mem_collisions = 0;
        dist_points = 0;
        number_steps_collect = 0;
        number_steps_locate = 0;
        number_steps = 0;
        dist_cols = 0;
        calendar_time = 0;
        total_time = 0;
        cycles = 0;

        // Setting number of cores
        insts_stats[i].N_OF_CORES = (uint16_t)n_cores;
        // Setting Hansel&Gretel optimization
        insts_stats[i].HANSEL_GRETEL = hansel_gretel;
        insts_stats[i].MAX_CRUMBS = max_crumbs;

        if (arg_delta != 0)
        {
            insts_stats[i].delta = arg_delta;
        }
        printf("delta = %ld\n", insts_stats[i].delta);

        init_shared_state(insts_stats + i, &S);
        benchmark_precomputation(insts_stats + i, &S);
        init_shared_state_second_step(insts_stats + i, &S);

        printf("\n----------------------------------------------------------------------------------------\n");
        printf("\nInstance:\t");
        printf("e = %" PRIu64 "\t    ", (uint64_t)insts_stats[i].e);
        if (arg_w != 0)
        {
            insts_stats[i].MEMORY_LOG_SIZE = arg_w;
        }
        printf("w = %" PRIu64 "\t", (uint64_t)insts_stats[i].MEMORY_LOG_SIZE);
        printf(_ALPHA_CHAR);
        printf(" = %.2f\t", insts_stats[i].ALPHA);
        printf(_BETA_CHAR);
        printf(" = %.2f\t", insts_stats[i].BETA);
        printf(_GAMMA_CHAR);
        printf(" = %.2f\t", insts_stats[i].GAMMA);
        printf(_DELTA_CHAR);
        printf(" = %" PRIu64 "\t", (uint64_t)insts_stats[i].delta);
        printf("modulus = %s", insts_stats[0].MODULUS);
        printf("\n\n");
        printf("Number of iterations averaged over: \t\t%" PRIu64 "\n", (uint64_t)BENCH_LOOPS);
        printf("Memory: \t\t\t\t\t");
        printf("RAM\n");
        printf("Number of cores: \t\t\t\t%" PRIu64 "\n", n_cores);
        printf("Hansel & Gretel: \t\t\t\t%s", hansel_gretel ? "Yes" : "No");
        if (hansel_gretel)
            printf(", %" PRIu64 " crumbs\n", insts_stats[i].MAX_CRUMBS);
        else
            printf("\n");
        printf("Statistics only: \t\t\t\t");
        printf("%s\n", collect_stats ? "Yes (only running one function version)" : "No");

        S.MAX_FUNCTION_VERSIONS = 100000;
        S.collect_vow_stats = collect_stats;
        for (j = 0; j < BENCH_LOOPS; j++)
        {
            reset_shared_state(&S);

            if (j > 0)
            {
                S.initial_function_version += 1; /* Maintain across examples, but change */
            }
            S.PRNG_SEED += j; /* Different PRNG seed.. */

            cycles1 = cpucycles();
            vOW(&S);
            cycles2 = cpucycles();
            cycles = cycles + (cycles2 - cycles1);

            success &= S.success;
            random_functions += (uint64_t)S.final_avg_random_functions;
            collisions += S.collisions;
            mem_collisions += S.mem_collisions;
            dist_points += S.dist_points;
            number_steps_collect += S.number_steps_collect;
            number_steps_locate += S.number_steps_locate;
            number_steps += S.number_steps;
            dist_cols += S.dist_cols.size;
            calendar_time += S.wall_time;
            total_time += S.total_time;

#ifdef COLLECT_DATABASE_STATS
            printf("\n");
            printf("avg read  time: %.6f / %.1f = %.6f\n", S.debug_stats[1], S.debug_stats[0], S.debug_stats[1] / S.debug_stats[0]);
            printf("avg write time: %.6f / %.1f = %.6f\n", S.debug_stats[3], S.debug_stats[2], S.debug_stats[3] / S.debug_stats[2]);
            printf("\n");
#endif

            if (!collect_stats)
            {
                printf("\n Iteration %d", j);
                if (S.success)
                    printf(" COMPLETED using %.2f random functions and %.2f seconds", S.final_avg_random_functions, S.wall_time);
                else
                    printf(" INCOMPLETE. Used %.2f random functions and %.2f seconds", S.final_avg_random_functions, S.wall_time);
            }
        }
        free_shared_state(&S);

        if (collect_stats)
        {
            printf("\nNumber of function iterations (i): \t\t%.2f\n", (double)number_steps / (double)random_functions);
            printf("\t For collecting dist. points: \t\t%.2f (%.2f%%)\n",
                   ((double)number_steps_collect / (double)random_functions),
                   100 * ((double)number_steps_collect / (double)number_steps));
            printf("\t For locating collisions: \t\t%.2f (%.2f%%)\n",
                   ((double)number_steps_locate / (double)random_functions),
                   100 * ((double)number_steps_locate / (double)number_steps));
            printf("Number of collisions per function: \t\t%.2f (expected 1.3w = %.2f, ratio = %.2f)\n", ((double)collisions / (double)random_functions),
                   1.3 * pow(2, (double)insts_stats[i].MEMORY_LOG_SIZE),
                   ((double)collisions / (double)random_functions) / (1.3 * pow(2, (double)insts_stats[i].MEMORY_LOG_SIZE)));
            printf("Number of distinct collisions per function (c): %.2f (expected 1.1w = %.2f, ratio = %.2f)\n",
                   ((double)dist_cols / (double)random_functions),
                   1.1 * pow(2, (double)insts_stats[i].MEMORY_LOG_SIZE),
                   ((double)dist_cols / (double)random_functions) / (1.1 * pow(2, (double)insts_stats[i].MEMORY_LOG_SIZE)));
            printf("\n");
            printf("Expected number of function versions (n/(2c)): \t%.2f (expected 0.45n/w = %.2f, ratio = %.2f)\n",
                   pow(2, insts_stats[i].e - 1) / (2 * ((double)dist_cols / (double)random_functions)),
                   0.45 * pow(2, insts_stats[i].e - 1) / pow(2, (double)insts_stats[i].MEMORY_LOG_SIZE),
                   (pow(2, insts_stats[i].e - 1) / (2 * ((double)dist_cols / (double)random_functions))) /
                       (0.45 * pow(2, insts_stats[i].e - 1) / pow(2, (double)insts_stats[i].MEMORY_LOG_SIZE)));
            printf("Expected total run-time (in/(2c)): \t\t%.2f (expected %cn^3/w = %.2f, ratio = %.2f)\n",
                   ((double)number_steps / (double)random_functions) * pow(2, insts_stats[i].e - 1) / (2 * ((double)dist_cols / (double)random_functions)),
                   251,
                   sqrt(pow(pow(2, insts_stats[i].e - 1), 3) / pow(2, insts_stats[i].MEMORY_LOG_SIZE)),
                   (((double)number_steps / (double)random_functions) * pow(2, insts_stats[i].e - 1) / (2 * ((double)dist_cols / (double)random_functions))) / sqrt(pow(pow(2, insts_stats[i].e - 1), 3) / pow(2, insts_stats[i].MEMORY_LOG_SIZE)));
        }
        else
        {
            printf("\nAll tests successful: \t\t\t%s\n\n", success ? "Yes" : "No");
            printf("\n");
            printf("Number of function iterations: \t\t%.2f (expected sqrt(n^3/w) = %.2f, ratio = %.2f)\n",
                   (double)number_steps / (double)BENCH_LOOPS,
                   sqrt(pow(pow(2, insts_stats[i].e - 1), 3) / pow(2, insts_stats[i].MEMORY_LOG_SIZE)),
                   ((double)number_steps / (double)BENCH_LOOPS) / (sqrt(pow(pow(2, insts_stats[i].e - 1), 3) / pow(2, insts_stats[i].MEMORY_LOG_SIZE))));
            printf("\t For collecting dist. points: \t%.2f (%.2f%%)\n",
                   ((double)number_steps_collect / (double)BENCH_LOOPS),
                   100 * ((double)number_steps_collect / (double)number_steps));
            printf("\t For locating collisions: \t%.2f (%.2f%%)\n",
                   ((double)number_steps_locate / (double)BENCH_LOOPS),
                   100 * ((double)number_steps_locate / (double)number_steps));
            printf("Number of function versions: \t\t%.2f (expected 0.45n/w = %.2f, ratio = %.2f)\n",
                   (double)random_functions / (double)BENCH_LOOPS,
                   0.45 * pow(2, insts_stats[i].e - 1) / pow(2, insts_stats[i].MEMORY_LOG_SIZE),
                   ((double)random_functions / (double)BENCH_LOOPS) /
                       (0.45 * pow(2, insts_stats[i].e - 1) / pow(2, insts_stats[i].MEMORY_LOG_SIZE)));
            printf("Number of collisions per function: \t%.2f (expected 1.3w = %.2f, ratio = %.2f)\n",
                   ((double)collisions / (double)random_functions),
                   1.3 * pow(2, insts_stats[i].MEMORY_LOG_SIZE),
                   (((double)collisions / (double)random_functions)) / (1.3 * pow(2, insts_stats[i].MEMORY_LOG_SIZE)));
            printf("Number of mem_cols per func: \t%.2f\n", (double)mem_collisions / (double)random_functions);
        }
        printf("\nTotal time                        : %.2f", total_time);
        printf("\nTotal time    (avg per iteration) : %.2f", total_time / BENCH_LOOPS);
        printf("\nCalendar time                     : %.2f", calendar_time);
        printf("\nCalendar time (avg per iteration) : %.2f", calendar_time / BENCH_LOOPS);
        printf("\nCycle count   (avg per iteration) : %lld\n\n", cycles / BENCH_LOOPS);
    }
    return 0;
}

int main(int argc, char **argv)
{

    int Status = PASSED;
    uint64_t n_cores = 1;       // One core by default
    bool hansel_gretel = false; // Hansel&Gretel optimization is disabled by default
    uint64_t max_crumbs = 10;
    bool collect_stats = false; // Extra collection of stats is disabled by default
    bool help_flag = false;

    // avoid output buffering
    setvbuf(stdout, NULL, _IONBF, 0);

    for (int i = 0; i < argc - 1; i++)
    {
        if (argv[i + 1][0] != '-')
        {
            help_flag = true;
            goto help;
        }
        switch (argv[i + 1][1])
        {
        case 'n':
            n_cores = strtol(argv[i + 2], NULL, 10);
            i++;
            if (n_cores > 1000)
                help_flag = true;
            break;
        case 'H':
            hansel_gretel = true;
            max_crumbs = strtol(argv[i + 2], NULL, 10);
            i++;
            break;
        case 's':
            collect_stats = true;
            break;
        case 'h':
            help_flag = true;
            break;
        case 'd':
            arg_delta = strtol(argv[i + 2], NULL, 10);
            i++;
            break;
        case 'w':
            arg_w = strtol(argv[i + 2], NULL, 10);
            i++;
            break;
        default:
            help_flag = true;
            break;
        }
        if (help_flag)
        {
            goto help;
        }
    }

    if (arg_delta > (insts_stats[0].e - 2))
    {
        printf("Precomputation depth is set too high. Maximal possible value = %ld\n", (insts_stats[0].e - 2));
        return !PASSED;
    }

    omp_set_num_threads(n_cores);
    printf("numthreads: %lu\n", n_cores);
    Status = stats_vow(n_cores, hansel_gretel, collect_stats, max_crumbs); // Testing
    if (Status != PASSED)
    {
        printf("\n\n   Error detected while running attack... \n\n");
        return 1;
    }

help:
    if (help_flag)
    {
        printf("\n Usage:");
        printf("\n test_vOW_SIKE -n [N_CORES] -H [MAX_CRUMBS] -s -h \n");
        printf("\n -n : number of cores (one core by default). Maximum 1000 cores.");
        printf("\n -H : Hansel&Gretel optimization on (off by default).");
        printf("\n -s : collection of attack stats on (off by default).");
        printf("\n -h : this help.\n\n");
    }

    return Status;
}
