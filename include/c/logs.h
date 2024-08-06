/**
 * @file logs.h
 * @author Daniel San Martin (dsanmartinreyes@gmail.com)
 * @brief Logging functions for the wildfire simulation.
 * @version 0.1
 * @date 2024-07-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef LOGS_H
#define LOGS_H

#include <stdio.h>
#include "parameters.h"

/**
 * @brief the parameters to a file or standard output.
 *
 * This function logs the given parameters to either a file or standard output,
 * depending on the value of the `to_file` parameter.
 *
 * @param parameters A pointer to the Parameters struct containing the parameters to log.
 * @param to_file    A flag indicating whether to log to a file (1) or standard output (0).
 */
void log_parameters_report(Parameters *parameters, int to_file);

/**
 * @brief Logs the parameters to a file and standard output.
 *
 * This function logs the parameters using the specified parameters.
 *
 * @param parameters The parameters to be logged.
 */
void log_parameters(Parameters *parameters);

/**
 * @brief Logs a message.
 *
 * This function logs a message using the specified parameters.
 * Show and save the message.
 *
 * @param parameters The parameters for logging.
 * @param message The message to be logged.
 */
void log_message(Parameters *parameters, const char *message);

/**
 * @brief Logs the current timestep information.
 *
 * This function logs the information related to the current timestep, including the parameters,
 * the timestep number, the current time, and the time taken for the timestep.
 *
 * @param parameters A pointer to the Parameters struct containing the simulation parameters.
 * @param n The timestep number.
 * @param t_n The current time.
 * @param step_time The time taken for the timestep.
 */
void log_timestep(Parameters *parameters, int n, double t_n, double step_time, double CFL, double T_min, double T_max, double Y_min, double Y_max);

#endif