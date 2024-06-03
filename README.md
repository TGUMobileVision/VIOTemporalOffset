# VIO Temporal Offset

***This package is the source code of VIO Temporal Offset Optimization.***

In this paper, an initialization strategy is proposed with respect to monocular visual-inertial localization, which optimizes the temporal offset between a camera and an inertial measurement unit (IMU). Moreover, persistent refinement is involved in front-end processes, and specific states are optimized in back-end processes. Firstly, with assumption that IMU timestamps are valid, influence of the temporal offset between camera and IMU timestamps is analyzed, and the temporal offset is derived from reprojection errors so as to align camera and IMU samplings. Then, taking into account camera samplings, both translation and rotation movement is calculated for projected 3D points, and variation of reprojected image features is analyzed with respect to the corresponding camera keyframe. At last, the temporal offset and specific states are estimated by optimizing reprojection errors, making that influence of the uncertain temporal offset is eliminated for visual-inertial localization. Comparative experiments are conducted to validate the performance of the proposed approach.

**The Related Papers**:

**Online Temporal Calibration for Monocular Visual-Inertial Systems**, Tong Qin, Shaojie Shen, IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS, 2018). [PDF](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8593603)

**For more technical details, please refer to**:

X. Gao, B. Li, W. Shi, and X. Zhang, Visual-inertial odometry systems with online temporal offset optimization, International Journal of Robotics and Automation, to appear, 2024.
