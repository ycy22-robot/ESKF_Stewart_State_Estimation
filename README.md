# Stewart Platform State Estimation (ESKF) - MATLAB Implementation

This repository contains the official MATLAB implementation for the paper **Accurate and Robust State Estimation for Stewart Platform via Error State Kalman Filter with IMU and Leg Length Fusion**. The code implements a state estimation algorithm for a Stewart platform based on the Error State Kalman Filter (ESKF), fusing IMU and leg length sensor data. The performance—accuracy, robustness, and computational efficiency—is compared against the Extended Kalman Filter (EKF) and Forward Kinematics (FK) methods.

## Contact

If you have any questions, please contact [yuecy22@mails.tsinghua.edu.cn](mailto:yuecy22@mails.tsinghua.edu.cn).

## Getting Started

1. **Clone the repository** 

   ```bash
   git clone https://github.com/ycy22-robot/ESKF_Stewart_State_Estimation.git
   ```

2. **Generate simulation data**

   ```matlab
   generate_robot_sensor_dataset.m
   ```

3. **Run the main script to perform state estimation and save results**

   ```matlab
   run_state_estimation_from_sim_data.m
   ```

4. **Plot comparative results**

   ```matlab
   mian_plot_nominal_sim.mlx
   ```

## Code Structure

| File/Script                            | Description                                                  |
| -------------------------------------- | ------------------------------------------------------------ |
| `generate_robot_sensor_dataset.m`      | Sets motion trajectories and noise parameters to generate simulation data (`rod_data.txt`, `imu_data.txt`), ground truth (`ground_truth.txt`), and bias ground truth (`true_bias_values.mat`) |
| `run_eskf_rod.m`                       | State estimation using ESKF algorithm; takes simulation data as input and outputs estimated states |
| `run_ekf_rod.m`                        | State estimation using EKF algorithm                         |
| `run_fk.m`                             | State estimation using Forward Kinematics                    |
| `run_state_estimation_from_sim_data.m` | Main script; loads data, runs all three algorithms, and saves results in `simulation_results_nominal.mat` |
| `mian_plot_nominal_sim.mlx`            | Plots comparison results of ESKF, EKF, and FK, consistent with experiments in Section 4 of the paper |
| Other utility scripts                  | - `FwKineQ.m`: Forward kinematics calculation<br>- `rotation_matrix_cov_to_euler_cov.m`: Converts ESKF rotation matrix covariance to Euler covariance<br>- `quaternion_cov_to_euler_cov.m`: Converts EKF quaternion covariance to Euler covariance |

**Notes:**

- To obtain initial IMU bias estimates in a static scenario, set all motion amplitudes to zero in `generate_robot_sensor_dataset.m`, run the script, and take the mean of the IMU outputs. See `imu_bias_init.txt`. The paper’s initial values are provided, but you may generate your own.
- The data used in Section 4 of the paper is included; simply run the main script to reproduce those experiments.
- To replicate the robustness tests in Section 5, adjust the noise parameters in `generate_robot_sensor_dataset.m` and repeat the above steps.

## FAQ

- **Q: How can I customize the simulation trajectory or noise parameters?**
   A: Edit the relevant settings in `generate_robot_sensor_dataset.m`.
- **Q: What if I encounter errors when running the code?**
   A: Ensure your MATLAB version and toolboxes meet the requirements. For persistent issues, please contact the author.
