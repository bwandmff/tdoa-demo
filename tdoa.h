/**
 * @file tdoa.h
 * @brief High-Efficiency TDOA (Time Difference of Arrival) Algorithm Implementation
 * @author C Programming Expert
 * @version 1.0
 * @date 2026-03-17
 *
 * This implementation provides efficient TDOA-based positioning algorithms
 * following C99/C11 standards with proper memory management.
 */

#ifndef TDOA_H_
#define TDOA_H_

#include <stddef.h>   /* for size_t */
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Maximum number of receivers supported
 */
#define TDOA_MAX_RECEIVERS 8

/**
 * @brief Maximum number of iterations for iterative algorithms
 */
#define TDOA_MAX_ITERATIONS 100

/**
 * @brief Default tolerance for convergence
 */
#define TDOA_DEFAULT_TOLERANCE 1e-6

/**
 * @brief 2D point structure
 */
typedef struct {
    double x;
    double y;
} tdoa_point_t;

/**
 * @brief 3D point structure
 */
typedef struct {
    double x;
    double y;
    double z;
} tdoa_point3d_t;

/**
 * @brief Receiver configuration
 */
typedef struct {
    tdoa_point3d_t position;  /**< Receiver position */
    double clock_bias;         /**< Clock bias offset */
    bool is_valid;            /**< Whether receiver is active */
} tdoa_receiver_t;

/**
 * @brief TDOA measurement structure
 */
typedef struct {
    size_t receiver_i;        /**< First receiver index */
    size_t receiver_j;        /**< Second receiver index */
    double time_difference;   /**< TDOA value in seconds */
    double noise_std;         /**< Measurement noise standard deviation */
} tdoa_measurement_t;

/**
 * @brief TDOA solver configuration
 */
typedef struct {
    double tolerance;          /**< Convergence tolerance */
    size_t max_iterations;    /**< Maximum iterations */
    bool use_weighting;       /**< Enable weighted least squares */
} tdoa_config_t;

/**
 * @brief TDOA solver context
 */
typedef struct {
    tdoa_receiver_t receivers[TDOA_MAX_RECEIVERS];
    size_t num_receivers;

    tdoa_measurement_t *measurements;
    size_t num_measurements;
    size_t max_measurements;  /**< Pre-allocated buffer size */

    tdoa_config_t config;

    double speed;  /**< Signal propagation speed (m/s), default: speed of light */

    /* Pre-allocated working buffers */
    double *matrix_a;
    double *vector_b;
    double *result;
    size_t buffer_size;
} tdoa_solver_t;

/**
 * @brief Initialize TDOA solver with default configuration
 * @param solver Pointer to solver context
 * @return 0 on success, -1 on error
 */
int tdoa_solver_init(tdoa_solver_t *solver);

/**
 * @brief Configure TDOA solver
 * @param solver Pointer to solver context
 * @param config Configuration parameters
 * @return 0 on success, -1 on error
 */
int tdoa_solver_configure(tdoa_solver_t *solver, const tdoa_config_t *config);

/**
 * @brief Set signal propagation speed
 * @param solver Pointer to solver context
 * @param speed Propagation speed in m/s (e.g., SPEED_OF_LIGHT or SPEED_OF_SOUND)
 */
void tdoa_solver_set_speed(tdoa_solver_t *solver, double speed);

/**
 * @brief Add a receiver to the solver
 * @param solver Pointer to solver context
 * @param position Receiver position
 * @param clock_bias Clock bias offset
 * @return Receiver index on success, -1 on error
 */
int tdoa_add_receiver(tdoa_solver_t *solver, const tdoa_point3d_t *position, double clock_bias);

/**
 * @brief Add a TDOA measurement
 * @param solver Pointer to solver context
 * @param receiver_i First receiver index
 * @param receiver_j Second receiver index
 * @param time_difference TDOA value
 * @param noise_std Measurement noise
 * @return 0 on success, -1 on error
 */
int tdoa_add_measurement(tdoa_solver_t *solver, size_t receiver_i, size_t receiver_j,
                         double time_difference, double noise_std);

/**
 * @brief Solve TDOA positioning using Fang's algorithm (closed-form, 2D)
 * @param solver Pointer to solver context
 * @param result Output position result
 * @return 0 on success, -1 on error
 */
int tdoa_solve_fang(const tdoa_solver_t *solver, tdoa_point3d_t *result);

/**
 * @brief Solve TDOA positioning using Chan algorithm (iterative WLS, 2D/3D)
 * @param solver Pointer to solver context
 * @param result Output position result
 * @return 0 on success, -1 on error
 */
int tdoa_solve_chan(const tdoa_solver_t *solver, tdoa_point3d_t *result);

/**
 * @brief Solve TDOA using Taylor series expansion (iterative, 2D/3D)
 * @param solver Pointer to solver context
 * @param initial_guess Initial position estimate
 * @param result Output position result
 * @return 0 on success, -1 on error
 */
int tdoa_solve_taylor(const tdoa_solver_t *solver, const tdoa_point3d_t *initial_guess,
                      tdoa_point3d_t *result);

/**
 * @brief Generate simulated TDOA measurements for testing
 * @param solver Pointer to solver context
 * @param true_position True emitter position
 * @param noise_std Add Gaussian noise with this std dev
 * @return 0 on success, -1 on error
 */
int tdoa_generate_measurements(tdoa_solver_t *solver, const tdoa_point3d_t *true_position,
                                double noise_std);

/**
 * @brief Get signal propagation speed (speed of sound or light)
 * @return Propagation speed in m/s
 */
double tdoa_get_speed_of_light(void);

/**
 * @brief Get speed of sound at given temperature
 * @param temperature_c Temperature in Celsius
 * @return Speed of sound in m/s
 */
double tdoa_get_speed_of_sound(double temperature_c);

/**
 * @brief Free solver resources
 * @param solver Pointer to solver context
 */
void tdoa_solver_destroy(tdoa_solver_t *solver);

/**
 * @brief Reset solver state
 * @param solver Pointer to solver context
 */
void tdoa_solver_reset(tdoa_solver_t *solver);

#ifdef __cplusplus
}
#endif

#endif /* TDOA_H_ */
