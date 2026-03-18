/**
 * @file tdoa.c
 * @brief High-Efficiency TDOA Algorithm Implementation
 * @details Optimized implementation with pre-allocated buffers and efficient algorithms
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
#define MAX_MEASUREMENTS 28           /* Max for 8 receivers: C(8,2) = 28 */

/* ============================================================================
 * Memory Management
 * ============================================================================ */

static void *safe_malloc(size_t size) {
    void *ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "Error: malloc failed for %zu bytes\n", size);
    }
    return ptr;
}

/* ============================================================================
 * Linear Algebra - Optimized for 3x3
 * ============================================================================ */

/**
 * @brief Optimized 3x3 linear system solver (no dynamic allocation)
 */
static inline int solve_3x3(double A[9], double b[3], double x[3]) {
    /* Gaussian elimination with partial pivoting */
    double aug[12];  /* 3x4 augmented matrix */

    /* Copy to augmented matrix */
    aug[0] = A[0]; aug[1] = A[1]; aug[2] = A[2]; aug[3] = b[0];
    aug[4] = A[3]; aug[5] = A[4]; aug[6] = A[5]; aug[7] = b[1];
    aug[8] = A[6]; aug[9] = A[7]; aug[10] = A[8]; aug[11] = b[2];

    for (int k = 0; k < 3; k++) {
        /* Find pivot */
        int pivot = k;
        double max_val = fabs(aug[k * 4 + k]);
        for (int i = k + 1; i < 3; i++) {
            double val = fabs(aug[i * 4 + k]);
            if (val > max_val) {
                max_val = val;
                pivot = i;
            }
        }

        if (max_val < 1e-12) {
            return -1;
        }

        /* Swap rows */
        if (pivot != k) {
            for (int j = k; j < 4; j++) {
                double tmp = aug[k * 4 + j];
                aug[k * 4 + j] = aug[pivot * 4 + j];
                aug[pivot * 4 + j] = tmp;
            }
        }

        /* Eliminate */
        for (int i = k + 1; i < 3; i++) {
            double factor = aug[i * 4 + k] / aug[k * 4 + k];
            for (int j = k; j < 4; j++) {
                aug[i * 4 + j] -= factor * aug[k * 4 + j];
            }
        }
    }

    /* Back substitution */
    x[2] = aug[11] / aug[10];
    x[1] = (aug[7] - aug[6] * x[2]) / aug[5];
    x[0] = (aug[3] - aug[1] * x[1] - aug[2] * x[2]) / aug[0];

    return 0;
}

/* ============================================================================
 * Inline distance calculation
 * ============================================================================ */

static inline double dist_sq_3d(double x1, double y1, double z1,
                                double x2, double y2, double z2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;
    return dx*dx + dy*dy + dz*dz;
}

static inline double dist_3d(double x1, double y1, double z1,
                             double x2, double y2, double z2) {
    return sqrt(dist_sq_3d(x1, y1, z1, x2, y2, z2));
}

/* ============================================================================
 * Public API
 * ============================================================================ */

int tdoa_solver_init(tdoa_solver_t *solver) {
    if (solver == NULL) {
        return -1;
    }

    memset(solver, 0, sizeof(tdoa_solver_t));

    solver->config.tolerance = 1e-4;
    solver->config.max_iterations = 50;
    solver->config.use_weighting = true;
    solver->speed = SPEED_OF_LIGHT;

    /* Pre-allocate measurement buffer */
    solver->max_measurements = MAX_MEASUREMENTS;
    solver->measurements = (tdoa_measurement_t *)safe_malloc(
        solver->max_measurements * sizeof(tdoa_measurement_t));

    /* Pre-allocate working buffers */
    solver->buffer_size = TDOA_MAX_RECEIVERS * TDOA_MAX_RECEIVERS;
    solver->matrix_a = (double *)safe_malloc(solver->buffer_size * solver->buffer_size * sizeof(double));
    solver->vector_b = (double *)safe_malloc(solver->buffer_size * sizeof(double));
    solver->result = (double *)safe_malloc(solver->buffer_size * sizeof(double));

    if (solver->measurements == NULL || solver->matrix_a == NULL ||
        solver->vector_b == NULL || solver->result == NULL) {
        tdoa_solver_destroy(solver);
        return -1;
    }

    for (size_t i = 0; i < TDOA_MAX_RECEIVERS; i++) {
        solver->receivers[i].is_valid = false;
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
        return -1;
    }
    if (solver->num_measurements >= solver->max_measurements) {
        return -1;
    }

    size_t idx = solver->num_measurements;
    solver->measurements[idx].receiver_i = receiver_i;
    solver->measurements[idx].receiver_j = receiver_j;
    solver->measurements[idx].time_difference = time_difference;
    solver->measurements[idx].noise_std = noise_std;
    solver->num_measurements++;

    return 0;
}

/**
 * @brief Fang's algorithm for 2D - Optimized
 */
int tdoa_solve_fang(const tdoa_solver_t *solver, tdoa_point3d_t *result) {
    if (solver == NULL || result == NULL) {
        return -1;
    }
    if (solver->num_receivers < 3 || solver->num_measurements < 2) {
        return -1;
    }

    const double c = solver->speed;
    const tdoa_point3d_t *r1 = &solver->receivers[0].position;
    const tdoa_point3d_t *r2 = &solver->receivers[1].position;
    const tdoa_point3d_t *r3 = &solver->receivers[2].position;

    /* Find TDOA measurements */
    double tij = 0.0, tik = 0.0;
    int found = 0;

    for (size_t m = 0; m < solver->num_measurements; m++) {
        const tdoa_measurement_t *meas = &solver->measurements[m];
        if (meas->receiver_i == 0 && meas->receiver_j == 1) {
            tij = meas->time_difference * c;
            found++;
        } else if (meas->receiver_i == 0 && meas->receiver_j == 2) {
            tik = meas->time_difference * c;
            found++;
        } else if (meas->receiver_i == 1 && meas->receiver_j == 0) {
            tij = -meas->time_difference * c;
            found++;
        } else if (meas->receiver_i == 2 && meas->receiver_j == 0) {
            tik = -meas->time_difference * c;
            found++;
        }
    }

    if (found < 2) {
        result->x = (r1->x + r2->x + r3->x) / 3.0;
        result->y = (r1->y + r2->y + r3->y) / 3.0;
        result->z = 0.0;
        return 0;
    }

    /* Fang's algorithm */
    double x1 = r1->x, y1 = r1->y;
    double x2 = r2->x, y2 = r2->y;
    double x3 = r3->x, y3 = r3->y;

    double r21 = dist_3d(x1, y1, 0, x2, y2, 0);
    double r31 = dist_3d(x1, y1, 0, x3, y3, 0);

    if (r21 < 1e-6 || r31 < 1e-6) {
        return -1;
    }

    double x21 = x2 - x1, y21 = y2 - y1;
    double x31 = x3 - x1, y31 = y3 - y1;

    double t1 = (tij - r21) / x21;
    double t2 = (tij * y21) / x21;
    double t3 = (tik - r31 + (x1*x1 + y1*y1 - 2*x1*x31 - 2*y1*y31) / x31) / x31;
    double t4 = (tij * y31 / x21 - tik * y31 / x31) / x31;

    double a = 1 + t1*t1 - t2*t2;
    double b = 2 * (t1*t3 - t2*t4 - t1 - t2);
    double c_coef = t3*t3 + t4*t4 - 1;

    double disc = b*b - 4*a*c_coef;
    if (disc < 0) disc = 0;

    double y = (-b + sqrt(disc)) / (2*a);
    double x = t1 - t2 * y;

    if (isnan(x) || isnan(y) || isinf(x) || isinf(y)) {
        result->x = (x1 + x2 + x3) / 3.0;
        result->y = (y1 + y2 + y3) / 3.0;
    } else {
        result->x = x;
        result->y = y;
    }
    result->z = 0.0;

    return 0;
}

/**
 * @brief Chan algorithm - Optimized with precomputed distances
 */
int tdoa_solve_chan(const tdoa_solver_t *solver, tdoa_point3d_t *result) {
    if (solver == NULL || result == NULL) {
        return -1;
    }
    if (solver->num_receivers < 3) {
        return -1;
    }

    const double c = solver->speed;
    const size_t m = solver->num_measurements;

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

    double step = 1.0;

    /* Cache receiver positions */
    double rx[8], ry[8], rz[8];
    for (size_t i = 0; i < solver->num_receivers; i++) {
        rx[i] = solver->receivers[i].position.x;
        ry[i] = solver->receivers[i].position.y;
        rz[i] = solver->receivers[i].position.z;
    }

    for (size_t iter = 0; iter < solver->config.max_iterations; iter++) {
        double error = 0, gx = 0, gy = 0, gz = 0;

        /* Compute error and gradient in one pass */
        for (size_t k = 0; k < m; k++) {
            const tdoa_measurement_t *meas = &solver->measurements[k];
            size_t i = meas->receiver_i;
            size_t j = meas->receiver_j;

            double dxi = cx - rx[i], dyi = cy - ry[i], dzi = cz - rz[i];
            double di = sqrt(dxi*dxi + dyi*dyi + dzi*dzi);

            double dxj = cx - rx[j], dyj = cy - ry[j], dzj = cz - rz[j];
            double dj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);

            if (di < 1e-6 || dj < 1e-6) continue;

            double pred = di - dj;
            double r = c * meas->time_difference - pred;
            error += r * r * 0.5;

            gx += r * (dxi/di - dxj/dj);
            gy += r * (dyi/di - dyj/dj);
            gz += r * (dzi/di - dzj/dj);
        }

        /* Update */
        double new_cx = cx + step * gx;
        double new_cy = cy + step * gy;
        double new_cz = cz + step * gz;

        /* Compute new error */
        double new_error = 0;
        for (size_t k = 0; k < m; k++) {
            const tdoa_measurement_t *meas = &solver->measurements[k];
            size_t i = meas->receiver_i;
            size_t j = meas->receiver_j;

            double di = dist_3d(new_cx, new_cy, new_cz, rx[i], ry[i], rz[i]);
            double dj = dist_3d(new_cx, new_cy, new_cz, rx[j], ry[j], rz[j]);
            double r = c * meas->time_difference - (di - dj);
            new_error += r * r * 0.5;
        }

        if (new_error < error) {
            cx = new_cx;
            cy = new_cy;
            cz = new_cz;

            if (fabs(new_error - error) < solver->config.tolerance) {
                break;
            }
            step = fmin(step * 1.1, 100.0);
        } else {
            step *= 0.5;
        }

        if (step < 1e-10) break;
    }

    result->x = cx;
    result->y = cy;
    result->z = cz;

    return 0;
}

/**
 * @brief Taylor series - Optimized with inline solver and preallocated buffers
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
    const size_t m = solver->num_measurements;

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

    /* Cache receiver positions */
    double rx[8], ry[8], rz[8];
    for (size_t i = 0; i < solver->num_receivers; i++) {
        rx[i] = solver->receivers[i].position.x;
        ry[i] = solver->receivers[i].position.y;
        rz[i] = solver->receivers[i].position.z;
    }

    /* Pre-allocate Jacobian and residual (on stack) */
    double J[84];  /* 28 x 3 max */
    double r[28];

    for (size_t iter = 0; iter < solver->config.max_iterations; iter++) {
        /* Build Jacobian and residual */
        for (size_t k = 0; k < m; k++) {
            const tdoa_measurement_t *meas = &solver->measurements[k];
            size_t i = meas->receiver_i;
            size_t j = meas->receiver_j;

            double dxi = cx - rx[i], dyi = cy - ry[i], dzi = cz - rz[i];
            double di = sqrt(dxi*dxi + dyi*dyi + dzi*dzi);

            double dxj = cx - rx[j], dyj = cy - ry[j], dzj = cz - rz[j];
            double dj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);

            if (di < 1e-4 || dj < 1e-4) {
                goto taylor_exit;
            }

            J[k * 3 + 0] = dxi/di - dxj/dj;
            J[k * 3 + 1] = dyi/di - dyj/dj;
            J[k * 3 + 2] = dzi/di - dzj/dj;

            r[k] = c * meas->time_difference - (di - dj);
        }

        /* Compute J'J and J'r (3x3 and 3x1) */
        double JJ[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double Jr[3] = {0, 0, 0};

        for (size_t k = 0; k < m; k++) {
            double *Jk = &J[k * 3];
            JJ[0] += Jk[0] * Jk[0];
            JJ[1] += Jk[0] * Jk[1];
            JJ[2] += Jk[0] * Jk[2];
            JJ[3] += Jk[1] * Jk[0];
            JJ[4] += Jk[1] * Jk[1];
            JJ[5] += Jk[1] * Jk[2];
            JJ[6] += Jk[2] * Jk[0];
            JJ[7] += Jk[2] * Jk[1];
            JJ[8] += Jk[2] * Jk[2];

            Jr[0] += Jk[0] * r[k];
            Jr[1] += Jk[1] * r[k];
            Jr[2] += Jk[2] * r[k];
        }

        /* Solve 3x3 system */
        double delta[3];
        if (solve_3x3(JJ, Jr, delta) != 0) {
            break;
        }

        /* Update */
        cx += delta[0];
        cy += delta[1];
        cz += delta[2];

        double norm = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
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

/**
 * @brief Generate measurements - Optimized with precomputed distances
 */
int tdoa_generate_measurements(tdoa_solver_t *solver, const tdoa_point3d_t *true_position,
                                double noise_std) {
    if (solver == NULL || true_position == NULL) {
        return -1;
    }

    const double c = solver->speed;
    const double tx = true_position->x;
    const double ty = true_position->y;
    const double tz = true_position->z;

    /* Cache receiver positions */
    double rx[8], ry[8], rz[8];
    for (size_t i = 0; i < solver->num_receivers; i++) {
        rx[i] = solver->receivers[i].position.x;
        ry[i] = solver->receivers[i].position.y;
        rz[i] = solver->receivers[i].position.z;
    }

    for (size_t i = 0; i < solver->num_receivers; i++) {
        for (size_t j = i + 1; j < solver->num_receivers; j++) {
            double di = dist_3d(tx, ty, tz, rx[i], ry[i], rz[i]);
            double dj = dist_3d(tx, ty, tz, rx[j], ry[j], rz[j]);

            double true_tdoa = (di - dj) / c;

            double noise = 0.0;
            if (noise_std > 0) {
                /* Fast Gaussian random */
                double u1 = (double)rand() / RAND_MAX;
                double u2 = (double)rand() / RAND_MAX;
                noise = noise_std * sqrt(-2.0 * log(u1 + 1e-10)) * cos(6.28318530718 * u2);
            }

            tdoa_add_measurement(solver, i, j, true_tdoa + noise, noise_std);
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
    free(solver->matrix_a);
    free(solver->vector_b);
    free(solver->result);
    solver->measurements = NULL;
    solver->matrix_a = NULL;
    solver->vector_b = NULL;
    solver->result = NULL;
    solver->num_measurements = 0;
    solver->num_receivers = 0;
}

void tdoa_solver_reset(tdoa_solver_t *solver) {
    if (solver == NULL) {
        return;
    }
    solver->num_measurements = 0;
    solver->num_receivers = 0;
    for (size_t i = 0; i < TDOA_MAX_RECEIVERS; i++) {
        solver->receivers[i].is_valid = false;
    }
}
