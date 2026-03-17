/**
 * @file main.c
 * @brief TDOA Algorithm Demonstration Program
 * @details Demonstrates high-efficiency TDOA positioning algorithms
 *          with Fang's, Chan's, and Taylor series methods.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "tdoa.h"

/* ============================================================================
 * Test Configuration
 * ============================================================================ */

#define NUM_RUNS 100

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

static void print_point(const char *label, const tdoa_point3d_t *point) {
    printf("%s: (%.6f, %.6f, %.6f)\n", label, point->x, point->y, point->z);
}

static double compute_error(const tdoa_point3d_t *estimated, const tdoa_point3d_t *true_pos) {
    double dx = estimated->x - true_pos->x;
    double dy = estimated->y - true_pos->y;
    double dz = estimated->z - true_pos->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

static void print_performance(const char *algorithm, const tdoa_point3d_t *estimated,
                              const tdoa_point3d_t *true_pos, double elapsed_ms) {
    double error = compute_error(estimated, true_pos);
    printf("=== %s Results ===\n", algorithm);
    print_point("  Estimated", estimated);
    print_point("  True Position", true_pos);
    printf("  Error: %.6f m\n", error);
    printf("  Time: %.3f ms\n\n", elapsed_ms);
}

/* ============================================================================
 * Test Cases
 * ============================================================================ */

static int test_basic_2d(void) {
    printf("\n");
    printf("========================================\n");
    printf("Test 1: Basic 2D TDOA (3 Receivers)\n");
    printf("========================================\n");

    tdoa_solver_t solver;
    if (tdoa_solver_init(&solver) != 0) {
        return -1;
    }

    /* Define receiver positions (2D, z=0) */
    tdoa_point3d_t r1 = {0.0, 0.0, 0.0};
    tdoa_point3d_t r2 = {100.0, 0.0, 0.0};
    tdoa_point3d_t r3 = {50.0, 86.6, 0.0};

    tdoa_add_receiver(&solver, &r1, 0.0);
    tdoa_add_receiver(&solver, &r2, 0.0);
    tdoa_add_receiver(&solver, &r3, 0.0);

    /* True emitter position */
    tdoa_point3d_t true_pos = {50.0, 30.0, 0.0};

    /* Generate measurements: 1ns TDOA noise (~0.3m) */
    double noise_std = 1e-9;
    tdoa_generate_measurements(&solver, &true_pos, noise_std);

    /* Solve using Fang's algorithm */
    tdoa_point3d_t result_fang;
    clock_t start = clock();
    int fang_result = tdoa_solve_fang(&solver, &result_fang);
    double elapsed_fang = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    if (fang_result == 0) {
        print_performance("Fang's Algorithm", &result_fang, &true_pos, elapsed_fang);
    }

    /* Solve using Chan's algorithm */
    tdoa_point3d_t result_chan;
    start = clock();
    int chan_result = tdoa_solve_chan(&solver, &result_chan);
    double elapsed_chan = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    if (chan_result == 0) {
        print_performance("Chan's Algorithm", &result_chan, &true_pos, elapsed_chan);
    }

    /* Solve using Taylor series (with centroid as initial guess) */
    tdoa_point3d_t initial_guess = {50.0, 28.0, 0.0};
    tdoa_point3d_t result_taylor;
    start = clock();
    int taylor_result = tdoa_solve_taylor(&solver, &initial_guess, &result_taylor);
    double elapsed_taylor = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    if (taylor_result == 0) {
        print_performance("Taylor Series", &result_taylor, &true_pos, elapsed_taylor);
    }

    tdoa_solver_destroy(&solver);
    return 0;
}

static int test_3d_positioning(void) {
    printf("\n");
    printf("========================================\n");
    printf("Test 2: 3D TDOA Positioning (4 Receivers)\n");
    printf("========================================\n");

    tdoa_solver_t solver;
    if (tdoa_solver_init(&solver) != 0) {
        return -1;
    }

    /* 3D receiver configuration */
    tdoa_point3d_t r1 = {0.0, 0.0, 0.0};
    tdoa_point3d_t r2 = {100.0, 0.0, 0.0};
    tdoa_point3d_t r3 = {50.0, 86.6, 0.0};
    tdoa_point3d_t r4 = {50.0, 28.9, 80.0};

    tdoa_add_receiver(&solver, &r1, 0.0);
    tdoa_add_receiver(&solver, &r2, 0.0);
    tdoa_add_receiver(&solver, &r3, 0.0);
    tdoa_add_receiver(&solver, &r4, 0.0);

    /* True emitter position in 3D */
    tdoa_point3d_t true_pos = {50.0, 30.0, 20.0};

    /* Generate measurements: 1ns noise */
    double noise_std = 1e-9;
    tdoa_generate_measurements(&solver, &true_pos, noise_std);

    /* Solve using Chan's algorithm */
    tdoa_point3d_t result_chan;
    clock_t start = clock();
    int chan_result = tdoa_solve_chan(&solver, &result_chan);
    double elapsed = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    if (chan_result == 0) {
        print_performance("Chan's Algorithm (3D)", &result_chan, &true_pos, elapsed);
    }

    /* Solve using Taylor series */
    tdoa_point3d_t initial_guess = {50.0, 30.0, 25.0};
    tdoa_point3d_t result_taylor;
    start = clock();
    int taylor_result = tdoa_solve_taylor(&solver, &initial_guess, &result_taylor);
    elapsed = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    if (taylor_result == 0) {
        print_performance("Taylor Series (3D)", &result_taylor, &true_pos, elapsed);
    }

    tdoa_solver_destroy(&solver);
    return 0;
}

static int test_monte_carlo(void) {
    printf("\n");
    printf("========================================\n");
    printf("Test 3: Monte Carlo Accuracy Test (%d runs)\n", NUM_RUNS);
    printf("========================================\n");

    srand((unsigned int)time(NULL));

    tdoa_solver_t solver;
    if (tdoa_solver_init(&solver) != 0) {
        return -1;
    }

    /* 4 receivers in square configuration */
    tdoa_point3d_t r1 = {0.0, 0.0, 0.0};
    tdoa_point3d_t r2 = {1000.0, 0.0, 0.0};
    tdoa_point3d_t r3 = {1000.0, 1000.0, 0.0};
    tdoa_point3d_t r4 = {0.0, 1000.0, 0.0};

    tdoa_add_receiver(&solver, &r1, 0.0);
    tdoa_add_receiver(&solver, &r2, 0.0);
    tdoa_add_receiver(&solver, &r3, 0.0);
    tdoa_add_receiver(&solver, &r4, 0.0);

    /* True position */
    tdoa_point3d_t true_pos = {500.0, 500.0, 100.0};

    /* 1ns TDOA noise (~0.3m) */
    double noise_std = 1e-9;

    double total_error = 0.0;
    double max_error = 0.0;
    double min_error = 1e10;
    int success_count = 0;

    clock_t start = clock();

    for (int run = 0; run < NUM_RUNS; run++) {
        tdoa_solver_reset(&solver);

        tdoa_add_receiver(&solver, &r1, 0.0);
        tdoa_add_receiver(&solver, &r2, 0.0);
        tdoa_add_receiver(&solver, &r3, 0.0);
        tdoa_add_receiver(&solver, &r4, 0.0);

        tdoa_generate_measurements(&solver, &true_pos, noise_std);

        /* Use Chan result as initial guess for Taylor */
        tdoa_point3d_t chan_result;
        tdoa_solve_chan(&solver, &chan_result);

        tdoa_point3d_t result;
        if (tdoa_solve_taylor(&solver, &chan_result, &result) == 0) {
            double error = compute_error(&result, &true_pos);
            total_error += error;
            if (error > max_error) max_error = error;
            if (error < min_error) min_error = error;
            success_count++;
        }
    }

    double elapsed = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    printf("Noise level: 1.0 ns TDOA\n");
    printf("Successful runs: %d/%d\n", success_count, NUM_RUNS);
    printf("Results:\n");
    printf("  Average error: %.3f m\n", total_error / success_count);
    printf("  Max error: %.3f m\n", max_error);
    printf("  Min error: %.3f m\n", min_error);
    printf("  Total time: %.2f ms (%.3f ms/run)\n\n", elapsed, elapsed / NUM_RUNS);

    tdoa_solver_destroy(&solver);
    return 0;
}

static int test_acoustic_tdoa(void) {
    printf("\n");
    printf("========================================\n");
    printf("Test 4: Acoustic TDOA (Speed of Sound)\n");
    printf("========================================\n");

    tdoa_solver_t solver;
    if (tdoa_solver_init(&solver) != 0) {
        return -1;
    }

    /* Microphone array (30cm apart) */
    tdoa_point3d_t m1 = {0.0, 0.0, 0.0};
    tdoa_point3d_t m2 = {0.30, 0.0, 0.0};
    tdoa_point3d_t m3 = {0.0, 0.30, 0.0};
    tdoa_point3d_t m4 = {0.0, 0.0, 0.30};

    tdoa_add_receiver(&solver, &m1, 0.0);
    tdoa_add_receiver(&solver, &m2, 0.0);
    tdoa_add_receiver(&solver, &m3, 0.0);
    tdoa_add_receiver(&solver, &m4, 0.0);

    /* True sound source position */
    tdoa_point3d_t true_pos = {1.5, 1.0, 0.5};

    double speed_of_sound = tdoa_get_speed_of_sound(20.0);
    printf("Speed of sound at 20C: %.2f m/s\n", speed_of_sound);

    /* Generate TDOA measurements: 1 microsecond noise */
    double noise_std = 1e-6;
    tdoa_generate_measurements(&solver, &true_pos, noise_std);

    /* Use Chan first then Taylor */
    tdoa_point3d_t chan_result;
    tdoa_solve_chan(&solver, &chan_result);

    tdoa_point3d_t initial_guess = chan_result;
    tdoa_point3d_t result;

    clock_t start = clock();
    int result_code = tdoa_solve_taylor(&solver, &initial_guess, &result);
    double elapsed = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    if (result_code == 0) {
        print_performance("Acoustic TDOA", &result, &true_pos, elapsed);
    }

    tdoa_solver_destroy(&solver);
    return 0;
}

static int test_custom_config(void) {
    printf("\n");
    printf("========================================\n");
    printf("Test 5: Custom Configuration (5 Receivers)\n");
    printf("========================================\n");

    tdoa_solver_t solver;
    if (tdoa_solver_init(&solver) != 0) {
        return -1;
    }

    /* Custom configuration */
    tdoa_config_t config = {
        .tolerance = 1e-6,
        .max_iterations = 100,
        .use_weighting = true
    };
    tdoa_solver_configure(&solver, &config);

    /* 5 receivers in pentagon */
    double radius = 500.0;
    for (int i = 0; i < 5; i++) {
        double angle = 2.0 * M_PI * i / 5.0;
        tdoa_point3d_t pos = {
            radius * cos(angle),
            radius * sin(angle),
            50.0
        };
        tdoa_add_receiver(&solver, &pos, 0.0);
    }

    /* True position */
    tdoa_point3d_t true_pos = {100.0, 200.0, 30.0};

    /* Generate measurements: 1ns noise */
    double noise_std = 1e-9;
    tdoa_generate_measurements(&solver, &true_pos, noise_std);

    /* Solve with Chan's algorithm */
    tdoa_point3d_t result;
    clock_t start = clock();
    int result_code = tdoa_solve_chan(&solver, &result);
    double elapsed = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    if (result_code == 0) {
        print_performance("Custom (Chan)", &result, &true_pos, elapsed);
    }

    /* Solve with Taylor series */
    tdoa_point3d_t initial_guess = result;  /* Use Chan result */
    start = clock();
    result_code = tdoa_solve_taylor(&solver, &initial_guess, &result);
    elapsed = (double)(clock() - start) * 1000.0 / CLOCKS_PER_SEC;

    if (result_code == 0) {
        print_performance("Custom (Taylor)", &result, &true_pos, elapsed);
    }

    tdoa_solver_destroy(&solver);
    return 0;
}

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    printf("========================================\n");
    printf("   High-Efficiency TDOA Demo\n");
    printf("   C Programming Expert Implementation\n");
    printf("========================================\n");

    printf("\nConstants:\n");
    printf("  Speed of light: %.0f m/s\n", tdoa_get_speed_of_light());
    printf("  Speed of sound (20C): %.2f m/s\n", tdoa_get_speed_of_sound(20.0));

    /* Run all tests */
    test_basic_2d();
    test_3d_positioning();
    test_monte_carlo();
    test_acoustic_tdoa();
    test_custom_config();

    printf("========================================\n");
    printf("   All tests completed!\n");
    printf("========================================\n");

    return 0;
}
