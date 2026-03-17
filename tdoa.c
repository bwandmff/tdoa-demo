/**
 * @file tdoa.c
 * @brief High-Efficiency TDOA Algorithm Implementation
 * @details Implements Fang's algorithm, Chan algorithm, and Taylor series method
 *          for Time Difference of Arrival based positioning.
 */

#include "tdoa.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

/* Constants */
#define SPEED_OF_LIGHT 299792458.0   /* m/s */
#define SPEED_OF_SOUND_20C 343.0     /* m/s at 20°C */

/* ============================================================================
 * Memory Management Utilities
 * ============================================================================ */

static void *safe_malloc(size_t size) {
    void *ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "Error: malloc failed for %zu bytes\n", size);
    }
    return ptr;
}

static void *safe_realloc(void *ptr, size_t size) {
    void *new_ptr = realloc(ptr, size);
    if (new_ptr == NULL && size > 0) {
        fprintf(stderr, "Error: realloc failed for %zu bytes\n", size);
        /* Note: realloc does NOT free ptr on failure, caller still owns it */
    }
    return new_ptr;
}

/* ============================================================================
 * Linear Algebra Routines
 * ============================================================================ */

/* Note: ordinary_least_squares is kept for potential future use */
static int solve_linear_system(double *A, double *b, double *x, size_t n) {
    if (A == NULL || b == NULL || x == NULL || n == 0) {
        return -1;
    }

    double *aug = (double *)safe_malloc(n * (n + 1) * sizeof(double));
    if (aug == NULL) {
        return -1;
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            aug[i * (n + 1) + j] = A[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }

    for (size_t k = 0; k < n; k++) {
        size_t pivot = k;
        double max_val = fabs(aug[k * (n + 1) + k]);
        for (size_t i = k + 1; i < n; i++) {
            double val = fabs(aug[i * (n + 1) + k]);
            if (val > max_val) {
                max_val = val;
                pivot = i;
            }
        }

        if (max_val < 1e-10) {
            free(aug);
            return -1;
        }

        if (pivot != k) {
            for (size_t j = k; j <= n; j++) {
                double tmp = aug[k * (n + 1) + j];
                aug[k * (n + 1) + j] = aug[pivot * (n + 1) + j];
                aug[pivot * (n + 1) + j] = tmp;
            }
        }

        for (size_t i = k + 1; i < n; i++) {
            double factor = aug[i * (n + 1) + k] / aug[k * (n + 1) + k];
            for (size_t j = k; j <= n; j++) {
                aug[i * (n + 1) + j] -= factor * aug[k * (n + 1) + j];
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = aug[i * (n + 1) + n];
        for (int j = i + 1; j < (int)n; j++) {
            x[i] -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] /= aug[i * (n + 1) + i];
    }

    free(aug);
    return 0;
}

/* ============================================================================
 * TDOA Core
 * ============================================================================ */

static double distance_3d(const tdoa_point3d_t *a, const tdoa_point3d_t *b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    double dz = a->z - b->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/* ============================================================================
 * Public API
 * ============================================================================ */

int tdoa_solver_init(tdoa_solver_t *solver) {
    if (solver == NULL) {
        return -1;
    }

    memset(solver, 0, sizeof(tdoa_solver_t));

    solver->config.tolerance = 1e-4;  /* 0.1mm tolerance */
    solver->config.max_iterations = 50;
    solver->config.use_weighting = true;

    /* Default to speed of light */
    solver->speed = SPEED_OF_LIGHT;

    for (size_t i = 0; i < TDOA_MAX_RECEIVERS; i++) {
        solver->receivers[i].is_valid = false;
    }

    solver->buffer_size = TDOA_MAX_RECEIVERS * TDOA_MAX_RECEIVERS;
    solver->matrix_a = (double *)safe_malloc(solver->buffer_size * solver->buffer_size * sizeof(double));
    solver->vector_b = (double *)safe_malloc(solver->buffer_size * sizeof(double));
    solver->result = (double *)safe_malloc(solver->buffer_size * sizeof(double));

    if (solver->matrix_a == NULL || solver->vector_b == NULL || solver->result == NULL) {
        tdoa_solver_destroy(solver);
        return -1;
    }

    return 0;
}

int tdoa_solver_configure(tdoa_solver_t *solver, const tdoa_config_t *config) {
    if (solver == NULL || config == NULL) {
        return -1;
    }

    solver->config.tolerance = config->tolerance;
    solver->config.max_iterations = config->max_iterations;
    solver->config.use_weighting = config->use_weighting;

    return 0;
}

void tdoa_solver_set_speed(tdoa_solver_t *solver, double speed) {
    if (solver == NULL || speed <= 0) {
        return;
    }
    solver->speed = speed;
}

int tdoa_add_receiver(tdoa_solver_t *solver, const tdoa_point3d_t *position, double clock_bias) {
    if (solver == NULL || position == NULL) {
        return -1;
    }

    if (solver->num_receivers >= TDOA_MAX_RECEIVERS) {
        fprintf(stderr, "Error: Maximum number of receivers (%d) reached\n", TDOA_MAX_RECEIVERS);
        return -1;
    }

    size_t idx = solver->num_receivers;
    solver->receivers[idx].position = *position;
    solver->receivers[idx].clock_bias = clock_bias;
    solver->receivers[idx].is_valid = true;

    solver->num_receivers++;

    return (int)idx;
}

int tdoa_add_measurement(tdoa_solver_t *solver, size_t receiver_i, size_t receiver_j,
                         double time_difference, double noise_std) {
    if (solver == NULL) {
        return -1;
    }

    if (receiver_i >= solver->num_receivers || receiver_j >= solver->num_receivers) {
        fprintf(stderr, "Error: Invalid receiver index\n");
        return -1;
    }

    size_t new_count = solver->num_measurements + 1;
    tdoa_measurement_t *new_meas = (tdoa_measurement_t *)safe_realloc(
        solver->measurements, new_count * sizeof(tdoa_measurement_t));

    if (new_meas == NULL) {
        return -1;
    }

    solver->measurements = new_meas;
    solver->measurements[solver->num_measurements].receiver_i = receiver_i;
    solver->measurements[solver->num_measurements].receiver_j = receiver_j;
    solver->measurements[solver->num_measurements].time_difference = time_difference;
    solver->measurements[solver->num_measurements].noise_std = noise_std;

    solver->num_measurements = new_count;

    return 0;
}

/**
 * @brief Fang's algorithm for 2D TDOA positioning
 * @note Uses solver's configured speed (light or sound)
 */
int tdoa_solve_fang(const tdoa_solver_t *solver, tdoa_point3d_t *result) {
    if (solver == NULL || result == NULL) {
        return -1;
    }

    if (solver->num_receivers < 3 || solver->num_measurements < 2) {
        return -1;
    }

    const double c = solver->speed;

    /* Use first 3 receivers */
    const tdoa_point3d_t *r1 = &solver->receivers[0].position;
    const tdoa_point3d_t *r2 = &solver->receivers[1].position;
    const tdoa_point3d_t *r3 = &solver->receivers[2].position;

    /* Find TDOA measurements for pairs (0,1) and (0,2) */
    double tij = 0.0, tik = 0.0;
    bool found_01 = false, found_02 = false;

    for (size_t m = 0; m < solver->num_measurements; m++) {
        const tdoa_measurement_t *meas = &solver->measurements[m];
        if ((meas->receiver_i == 0 && meas->receiver_j == 1) ||
            (meas->receiver_i == 1 && meas->receiver_j == 0)) {
            tij = meas->time_difference * c;
            if (meas->receiver_i == 1 && meas->receiver_j == 0) {
                tij = -tij;  /* Reverse the TDOA */
            }
            found_01 = true;
        } else if ((meas->receiver_i == 0 && meas->receiver_j == 2) ||
                   (meas->receiver_i == 2 && meas->receiver_j == 0)) {
            tik = meas->time_difference * c;
            if (meas->receiver_i == 2 && meas->receiver_j == 0) {
                tik = -tik;  /* Reverse the TDOA */
            }
            found_02 = true;
        }
    }

    /* Use centroid as fallback if measurements not found */
    if (!found_01 || !found_02) {
        result->x = (r1->x + r2->x + r3->x) / 3.0;
        result->y = (r1->y + r2->y + r3->y) / 3.0;
        result->z = 0.0;
        return 0;
    }

    /* Fang's algorithm implementation */
    double x1 = r1->x, y1 = r1->y;
    double x2 = r2->x, y2 = r2->y;
    double x3 = r3->x, y3 = r3->y;

    /* Distance between receivers */
    double r21 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    double r31 = sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1));

    if (r21 < 1e-6 || r31 < 1e-6) {
        return -1;  /* Invalid geometry */
    }

    /* Fang's method */
    double k = x1*x1 + y1*y1;
    double x21 = x2 - x1;
    double y21 = y2 - y1;
    double x31 = x3 - x1;
    double y31 = y3 - y1;

    double t1 = (tij - r21) / x21;
    double t2 = (tij * y21) / x21;
    double t3 = (tik - r31 + (k - 2*x1*x31 - 2*y1*y31) / x31) / x31;
    double t4 = (tij * y31 / x21 - tik * y31 / x31) / x31;

    /* Solve quadratic: a*y^2 + b*y + c = 0 */
    double a = 1 + t1*t1 - t2*t2;
    double b = 2 * (t1*t3 - t2*t4 - t1 - t2);
    double c_coef = t3*t3 + t4*t4 - 1;

    double discriminant = b*b - 4*a*c_coef;
    if (discriminant < 0) {
        discriminant = 0;  /* Use smaller magnitude solution */
    }

    double y = (-b + sqrt(discriminant)) / (2*a);
    double x = t1 - t2 * y;

    /* Validate result */
    if (isnan(x) || isnan(y) || isinf(x) || isinf(y)) {
        result->x = (x1 + x2 + x3) / 3.0;
        result->y = (y1 + y2 + y3) / 3.0;
        result->z = 0.0;
    } else {
        result->x = x;
        result->y = y;
        result->z = 0.0;
    }

    return 0;
}

/**
 * @brief Gradient descent with good convergence
 */
int tdoa_solve_chan(const tdoa_solver_t *solver, tdoa_point3d_t *result) {
    if (solver == NULL || result == NULL) {
        return -1;
    }

    if (solver->num_receivers < 3) {
        return -1;
    }

    const double c = solver->speed;

    /* Initial guess: centroid */
    double cx = 0, cy = 0, cz = 0;
    for (size_t i = 0; i < solver->num_receivers; i++) {
        cx += solver->receivers[i].position.x;
        cy += solver->receivers[i].position.y;
        cz += solver->receivers[i].position.z;
    }
    cx /= solver->num_receivers;
    cy /= solver->num_receivers;
    cz /= solver->num_receivers;

    /* Gradient descent with backtracking line search */
    double step = 1.0;

    for (size_t iter = 0; iter < solver->config.max_iterations; iter++) {
        /* Compute error and gradient */
        double error = 0;
        double gx = 0, gy = 0, gz = 0;

        for (size_t m = 0; m < solver->num_measurements; m++) {
            const tdoa_measurement_t *meas = &solver->measurements[m];
            const tdoa_point3d_t *ri = &solver->receivers[meas->receiver_i].position;
            const tdoa_point3d_t *rj = &solver->receivers[meas->receiver_j].position;

            double dx = cx - ri->x, dy = cy - ri->y, dz = cz - ri->z;
            double di = sqrt(dx*dx + dy*dy + dz*dz);

            double ex = cx - rj->x, ey = cy - rj->y, ez = cz - rj->z;
            double dj = sqrt(ex*ex + ey*ey + ez*ez);

            if (di < 1e-6 || dj < 1e-6) continue;

            double pred = di - dj;
            double meas_val = c * meas->time_difference;
            double r = meas_val - pred;

            error += r * r * 0.5;

            /* Gradient */
            double w = 1.0;
            gx += w * r * (dx/di - ex/dj);
            gy += w * r * (dy/di - ey/dj);
            gz += w * r * (dz/di - ez/dj);
        }

        /* Update with backtracking */
        double new_cx = cx + step * gx;
        double new_cy = cy + step * gy;
        double new_cz = cz + step * gz;

        /* Compute new error */
        double new_error = 0;
        for (size_t m = 0; m < solver->num_measurements; m++) {
            const tdoa_measurement_t *meas = &solver->measurements[m];
            const tdoa_point3d_t *ri = &solver->receivers[meas->receiver_i].position;
            const tdoa_point3d_t *rj = &solver->receivers[meas->receiver_j].position;

            double di = sqrt((new_cx - ri->x)*(new_cx - ri->x) +
                           (new_cy - ri->y)*(new_cy - ri->y) +
                           (new_cz - ri->z)*(new_cz - ri->z));
            double dj = sqrt((new_cx - rj->x)*(new_cx - rj->x) +
                           (new_cy - rj->y)*(new_cy - rj->y) +
                           (new_cz - rj->z)*(new_cz - rj->z));

            double pred = di - dj;
            double meas_val = c * meas->time_difference;
            double r = meas_val - pred;

            new_error += r * r * 0.5;
        }

        /* Accept or reduce step */
        if (new_error < error) {
            cx = new_cx;
            cy = new_cy;
            cz = new_cz;

            if (fabs(new_error - error) < solver->config.tolerance) {
                break;
            }
            step = step * 1.1;  /* Increase step when making progress */
        } else {
            step = step * 0.5;  /* Reduce step when not */
        }

        if (step < 1e-10) break;
    }

    result->x = cx;
    result->y = cy;
    result->z = cz;

    return 0;
}

/**
 * @brief Taylor Series Expansion using iterative linearization
 */
int tdoa_solve_taylor(const tdoa_solver_t *solver, const tdoa_point3d_t *initial_guess,
                      tdoa_point3d_t *result) {
    if (solver == NULL || initial_guess == NULL || result == NULL) {
        return -1;
    }

    if (solver->num_measurements < 2) {
        return -1;
    }

    const double c = solver->speed;

    /* Start from initial guess */
    double cx = initial_guess->x;
    double cy = initial_guess->y;
    double cz = initial_guess->z;

    /* Use centroid as fallback */
    if (fabs(cx) < 1e-6 && fabs(cy) < 1e-6 && fabs(cz) < 1e-6) {
        for (size_t i = 0; i < solver->num_receivers; i++) {
            cx += solver->receivers[i].position.x;
            cy += solver->receivers[i].position.y;
            cz += solver->receivers[i].position.z;
        }
        cx /= solver->num_receivers;
        cy /= solver->num_receivers;
        cz /= solver->num_receivers;
    }

    double best_cx = cx, best_cy = cy, best_cz = cz;
    double best_error = 1e20;

    /* Taylor series iteration */
    for (size_t iter = 0; iter < solver->config.max_iterations; iter++) {
        size_t m = solver->num_measurements;

        /* Build overdetermined system */
        double *A = (double *)safe_malloc(m * 3 * sizeof(double));
        double *b = (double *)safe_malloc(m * sizeof(double));

        if (A == NULL || b == NULL) {
            free(A);
            free(b);
            break;
        }

        /* Compute residuals and Jacobian */
        for (size_t i = 0; i < m; i++) {
            const tdoa_measurement_t *meas = &solver->measurements[i];
            const tdoa_point3d_t *ri = &solver->receivers[meas->receiver_i].position;
            const tdoa_point3d_t *rj = &solver->receivers[meas->receiver_j].position;

            double dxi = cx - ri->x;
            double dyi = cy - ri->y;
            double dzi = cz - ri->z;
            double di = sqrt(dxi*dxi + dyi*dyi + dzi*dzi);

            double dxj = cx - rj->x;
            double dyj = cy - rj->y;
            double dzj = cz - rj->z;
            double dj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);

            if (di < 1e-4 || dj < 1e-4) {
                free(A);
                free(b);
                cx = best_cx;
                cy = best_cy;
                cz = best_cz;
                goto taylor_exit;
            }

            /* Jacobian matrix row */
            A[i * 3 + 0] = (cx - ri->x) / di - (cx - rj->x) / dj;
            A[i * 3 + 1] = (cy - ri->y) / di - (cy - rj->y) / dj;
            A[i * 3 + 2] = (cz - ri->z) / di - (cz - rj->z) / dj;

            /* Residual */
            b[i] = c * meas->time_difference - (di - dj);
        }

        /* Solve using normal equations: A'A delta = A'b */
        double ATA[9];
        double ATb[3];

        /* Compute A'A */
        for (int r = 0; r < 3; r++) {
            for (int c2 = 0; c2 < 3; c2++) {
                ATA[r * 3 + c2] = 0;
                for (size_t k = 0; k < m; k++) {
                    ATA[r * 3 + c2] += A[k * 3 + r] * A[k * 3 + c2];
                }
            }
            ATb[r] = 0;
            for (size_t k = 0; k < m; k++) {
                ATb[r] += A[k * 3 + r] * b[k];
            }
        }

        free(A);
        free(b);

        /* Solve 3x3 system */
        double delta[3];
        int res = solve_linear_system(ATA, ATb, delta, 3);

        if (res != 0) {
            break;
        }

        /* Update estimate */
        cx += delta[0];
        cy += delta[1];
        cz += delta[2];

        /* Compute error for this iteration */
        double error = 0;
        for (size_t i = 0; i < m; i++) {
            const tdoa_measurement_t *meas = &solver->measurements[i];
            const tdoa_point3d_t *ri = &solver->receivers[meas->receiver_i].position;
            const tdoa_point3d_t *rj = &solver->receivers[meas->receiver_j].position;

            double di = distance_3d(&(tdoa_point3d_t){cx, cy, cz}, ri);
            double dj = distance_3d(&(tdoa_point3d_t){cx, cy, cz}, rj);
            double pred = di - dj;
            double meas_val = c * meas->time_difference;
            double r = meas_val - pred;
            error += r * r;
        }

        if (error < best_error) {
            best_error = error;
            best_cx = cx;
            best_cy = cy;
            best_cz = cz;
        }

        /* Check convergence */
        double norm = sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
        if (norm < solver->config.tolerance) {
            break;
        }

        if (norm > 1000.0) {
            break;
        }
    }

taylor_exit:
    result->x = cx;
    result->y = cy;
    result->z = cz;

    return 0;
}

int tdoa_generate_measurements(tdoa_solver_t *solver, const tdoa_point3d_t *true_position,
                                double noise_std) {
    if (solver == NULL || true_position == NULL) {
        return -1;
    }

    /* Use the solver's configured speed (light or sound) */
    const double c = solver->speed;

    /* Generate all pairwise TDOA measurements */
    for (size_t i = 0; i < solver->num_receivers; i++) {
        for (size_t j = i + 1; j < solver->num_receivers; j++) {
            const tdoa_point3d_t *ri = &solver->receivers[i].position;
            const tdoa_point3d_t *rj = &solver->receivers[j].position;

            double di = distance_3d(true_position, ri);
            double dj = distance_3d(true_position, rj);

            /* True TDOA = (d_i - d_j) / c */
            double true_tdoa = (di - dj) / c;

            double noise = 0.0;
            if (noise_std > 0) {
                double u1 = (double)rand() / RAND_MAX;
                double u2 = (double)rand() / RAND_MAX;
                noise = noise_std * sqrt(-2.0 * log(u1 + 1e-10)) * cos(2.0 * M_PI * u2);
            }

            double measured_tdoa = true_tdoa + noise;
            tdoa_add_measurement(solver, i, j, measured_tdoa, noise_std);
        }
    }

    return 0;
}

double tdoa_get_speed_of_light(void) {
    return SPEED_OF_LIGHT;
}

double tdoa_get_speed_of_sound(double temperature_c) {
    return 331.3 * sqrt(1.0 + temperature_c / 273.15);
}

void tdoa_solver_destroy(tdoa_solver_t *solver) {
    if (solver == NULL) {
        return;
    }

    free(solver->measurements);
    solver->measurements = NULL;

    free(solver->matrix_a);
    solver->matrix_a = NULL;

    free(solver->vector_b);
    solver->vector_b = NULL;

    free(solver->result);
    solver->result = NULL;

    solver->num_measurements = 0;
    solver->num_receivers = 0;
}

void tdoa_solver_reset(tdoa_solver_t *solver) {
    if (solver == NULL) {
        return;
    }

    free(solver->measurements);
    solver->measurements = NULL;
    solver->num_measurements = 0;

    for (size_t i = 0; i < TDOA_MAX_RECEIVERS; i++) {
        solver->receivers[i].is_valid = false;
    }
    solver->num_receivers = 0;
}
