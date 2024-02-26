#include "projection_td_factor.h"
#include "stdio.h"
using namespace std;
using std::ostream;
using std::istream;

Eigen::Matrix2d ProjectionTdFactor::sqrt_info;//sqrt_info为460/1.5的2*2对角矩阵[306.667 0;0 306.667]
double ProjectionTdFactor::sum_t;

ProjectionTdFactor::ProjectionTdFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_pts_j, 
                                       const Eigen::Vector2d &_velocity_i, const Eigen::Vector2d &_velocity_j,
                                       const double _td_i, const double _td_j, const double _row_i, const double _row_j) : 
                                       pts_i(_pts_i), pts_j(_pts_j), 
                                       td_i(_td_i), td_j(_td_j)
{
    velocity_i.x() = _velocity_i.x();
    velocity_i.y() = _velocity_i.y();
    velocity_i.z() = 0;
    velocity_j.x() = _velocity_j.x();
    velocity_j.y() = _velocity_j.y();
    velocity_j.z() = 0;
    row_i = _row_i - ROW / 2;   
    row_j = _row_j - ROW / 2;
    //在parameter.cpp里设置ROW = fsSettings["image_height"]; COL = fsSettings["image_width"];

#ifdef UNIT_SPHERE_ERROR
    Eigen::Vector3d b1, b2;
    Eigen::Vector3d a = pts_j.normalized();
    Eigen::Vector3d tmp(0, 0, 1);
    if(a == tmp)
        tmp << 1, 0, 0;
    b1 = (tmp - a * (a.transpose() * tmp)).normalized();
    b2 = a.cross(b1);
    tangent_base.block<1, 3>(0, 0) = b1.transpose();
    tangent_base.block<1, 3>(1, 0) = b2.transpose();
#endif
};

bool ProjectionTdFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    TicToc tic_toc;
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    double inv_dep_i = parameters[3][0];

    double td = parameters[4][0];

    //std::cout<<"inv_dep_i:  "<<inv_dep_i<<std::endl;
    //std::cout<<"parameters[3][1]： "<<parameters[3][1]<<std::endl;

    Eigen::Quaterniond Qij_c = (Qj * qic).inverse() * (Qi * qic);
    Eigen::Quaterniond Qc_i_(1, Qij_c.x()*td, Qij_c.y()*td, Qij_c.z()*td);
    /*std::cout<<"td_j="<<td_j<<std::endl;
    std::cout<<"td_i="<<td_i<<std::endl;
    std::cout<<"td="<<td<<std::endl;
    std::cout<<"Qij_c.x="<<Qij_c.x()<<std::endl;*/
    //std::cout<<"td - td_i="<<td - td_i<<std::endl;  //td-td_i=0   td=0  Qc_i_.x=-nan  Qc_i_.y=-nan  Qc_i_.z=-nan
    //std::cout<<"td="<<td<<std::endl;
    //std::cout<<"Qc_i_: "<<Qc_i_.w()<<" "<<Qc_i_.x()<<" "<Qc_i_.y()<<Qc_i_.z()<<" "<<std::endl;
    //std::cout<<"Qc_i_.y="<<std::endl;
    //std::cout<<"Qc_i_.z="<<"\n"<<std::endl;





    Eigen::Vector3d pts_i_td, pts_j_td;
    pts_i_td = pts_i - (td - td_i) * velocity_i;
    pts_j_td = pts_j - (td - td_j) * velocity_j;
    //pts_i_td = pts_i - (td - td_i + TR / ROW * row_i) * velocity_i;     //带有td变量
    //pts_j_td = pts_j - (td - td_j + TR / ROW * row_j) * velocity_j;
    Eigen::Vector3d pts_camera_i = pts_i_td / inv_dep_i;
    //std::cout<<"pts_i_td：\n"<<pts_i_td<<std::endl;
    //std::cout<<"pts_j_td：\n"<<pts_j_td<<std::endl;


    //Eigen::Vector3d Pi_ =  Eigen::Matrix3d::Identity() * Pi;
    //Eigen::Quaterniond Qi_ = (Qj* Qj.inverse())*Qi;
    //std::cout<<"sqrt_info：\n"<<sqrt_info<<std::endl;

/*
    //测试head<2>()是输出的什么
    Eigen::Vector3d P1 (1,2,3);
    std::cout<<"P1.head<2>=\n"<<P1.head<2>()<<std::endl;
    //测试完毕，P1.head<2>()输出P1的前两个元素
*/

/*
    //测试四元数乘法哪个是正确的
    Eigen::Quaterniond Q1 (1,2,3,4);
    Eigen::Quaterniond Q2 (5,6,7,8);
    Eigen::Matrix3d R1 = Q1.toRotationMatrix();
    Eigen::Matrix3d R2 = Q2.toRotationMatrix();
    Eigen::Matrix3d R3 = R1 * R2 ;
    Eigen::Quaterniond Q3(R3);      //转换成矩阵的乘法结果不对
    Eigen::Quaterniond Q4=Q1*Q2;    //输出的结果这么算是正确的
    std::cout<<"Q3.w=\n"<<Q3.w()<<std::endl;
    std::cout<<"Q3.x=\n"<<Q3.x()<<std::endl;
    std::cout<<"Q3.y=\n"<<Q3.y()<<std::endl;
    std::cout<<"Q3.z=\n"<<Q3.z()<<std::endl;
    std::cout<<"Q4.w=\n"<<Q4.w()<<std::endl;
    std::cout<<"Q4.x=\n"<<Q4.x()<<std::endl;
    std::cout<<"Q4.y=\n"<<Q4.y()<<std::endl;
    std::cout<<"Q4.z=\n"<<Q4.z()<<std::endl;
    //std::cout<<"R1=\n"<<R1<<std::endl;
    //std::cout<<"R2=\n"<<R2<<std::endl;
    //std::cout<<"R3=\n"<<R3<<std::endl;
 */ 


    //测试矩阵拆开算是否可以   不可以，为什么？？  测试之后这么乘是正确的啊
    //Eigen::Quaterniond Qi_c = Qi * qic;     //Q_ci^w
    //Eigen::Matrix3d Ri_c = Qi_c.toRotationMatrix();
    //Eigen::Vector3d Pi_c = Ri_c * tic;        //P_ci^w (1)
    //Eigen::Vector3d Pi_c = Qi_c * tic;          //       (2)两种表示方法效果一样
    //Eigen::Vector3d pts_w_ = Qi * qic * pts_camera_i + Qi * tic;   //拆开计算直接跑飞,难道是四元数这么相乘是不对的
    //Eigen::Vector3d pts_w_ = Qi_c * pts_camera_i + Pi_c;
    
    //std::cout<<"pts_w.x="<<pts_w.x()<<std::endl;
    //std::cout<<"pts_w.y="<<pts_w.y()<<std::endl;
    //std::cout<<"pts_w.z="<<pts_w.z()<<std::endl;

/* 
    //添加td
    Eigen::Quaterniond Qi_c = Qi * qic;
    Eigen::Quaterniond Qj_c = Qj * qic;
    Eigen::Quaterniond Qij_c = Qj_c.inverse() * Qi_c;

    Eigen::Quaterniond Qc_i_;
    Qc_i_.w()=1;
    Qc_i_.x()=Qij_c.x()/(td_j-td_i)*td;
    Qc_i_.y()=Qij_c.y()/(td_j-td_i)*td;
    Qc_i_.z()=Qij_c.z()/(td_j-td_i)*td; //这样把四元数转换过去

    Eigen::Quaterniond Qi_c_ = Qi_c * Qc_i_;
    Eigen::Vector3d Pi_c = Qi * tic;
    Eigen::Vector3d pts_w = Qi_c_ * pts_camera_i + Pi_c;
    */
    //Eigen::Vector3d pts_w = (Qi * qic) * pts_camera_i + Qi * tic + Pi;
    
       
        //Eigen::Quaterniond Qij_c (1,2,3,4);   //测试
        //std::cout<<"123td:"<<td<<std::endl;
    //Eigen::Quaterniond Qij_c=(Qj * qic).inverse() * (Qj * qic);
        /*std::cout<<"Qij_c.w="<<Qij_c.w()<<std::endl;
        cout<<"Qij_c.x="<<Qij_c.x()<<std::endl;
        cout<<"Qij_c.y="<<Qij_c.y()<<std::endl;
        cout<<"Qij_c.z="<<Qij_c.z()<<std::endl;*/
    //Eigen::Quaterniond Qc_i_(1,Qij_c.x()/(td_j-td_i)*td,Qij_c.y()/(td_j-td_i)*td,Qij_c.z()/(td_j-td_i)*td);
        //std::cout<<"Qc_i_="<<Qc_i_.w()<<Qc_i_.x()<<Qc_i_.y()<<Qc_i_.z()<<std::endl;

    //Eigen::Vector3d pts_w = (Qi * qic) * Qc_i_ * pts_camera_i + Qi * tic + Pi;
    

 

    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;  //源代码
    //std::cout<<"pts_w_拆开计算=\n"<<pts_w_<<std::endl;
    //std::cout<<"pts_w源代码=\n"<<pts_w<<std::endl;
    //Eigen::Vector3d pts_w = Qi * (qic * pts_camera_i + tic) + Pi;     //直接乘到一块这样也可以
    Eigen::Vector3d pts_w_ = Qc_i_ * pts_w;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w_ - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    Eigen::Map<Eigen::Vector2d> residual(residuals);

    //6.18 至此，应该就剩下雅克比矩阵了，先按照拆开算修改雅克比矩阵看看是不是还跑飞

#ifdef UNIT_SPHERE_ERROR 
    residual =  tangent_base * (pts_camera_j.normalized() - pts_j_td.normalized());
#else
    double dep_j = pts_camera_j.z();
    residual = (pts_camera_j / dep_j).head<2>() - pts_j_td.head<2>();
#endif

    residual = sqrt_info * residual;

    if (jacobians)
    {
        Eigen::Matrix3d Ri = Qi.toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();
        Eigen::Matrix3d ric = qic.toRotationMatrix();
        Eigen::Matrix3d Rc_i_ = Qc_i_.toRotationMatrix();   //将Qc_i_转化成矩阵
        Eigen::Matrix<double, 2, 3> reduce(2, 3);
#ifdef UNIT_SPHERE_ERROR
        double norm = pts_camera_j.norm();
        Eigen::Matrix3d norm_jaco;
        double x1, x2, x3;
        x1 = pts_camera_j(0);
        x2 = pts_camera_j(1);
        x3 = pts_camera_j(2);
        norm_jaco << 1.0 / norm - x1 * x1 / pow(norm, 3), - x1 * x2 / pow(norm, 3),            - x1 * x3 / pow(norm, 3),
                     - x1 * x2 / pow(norm, 3),            1.0 / norm - x2 * x2 / pow(norm, 3), - x2 * x3 / pow(norm, 3),
                     - x1 * x3 / pow(norm, 3),            - x2 * x3 / pow(norm, 3),            1.0 / norm - x3 * x3 / pow(norm, 3);
        reduce = tangent_base * norm_jaco;
#else
        reduce << 1. / dep_j, 0, -pts_camera_j(0) / (dep_j * dep_j),
            0, 1. / dep_j, -pts_camera_j(1) / (dep_j * dep_j);
#endif
        reduce = sqrt_info * reduce;

        if (jacobians[0])//i时刻的p和q
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix<double, 3, 6> jaco_i;
            jaco_i.leftCols<3>() = ric.transpose() * Rj.transpose() * Rc_i_;
            jaco_i.rightCols<3>() = ric.transpose() * Rj.transpose() * Rc_i_ * Ri * -Utility::skewSymmetric(pts_imu_i);//前面修改了这里就没有pts_imu_i声明了

            jacobian_pose_i.leftCols<6>() = reduce * jaco_i;
            jacobian_pose_i.rightCols<1>().setZero();
        }

        if (jacobians[1])//j时刻的p和q
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix<double, 3, 6> jaco_j;
            jaco_j.leftCols<3>() = ric.transpose() * -Rj.transpose();
            jaco_j.rightCols<3>() = ric.transpose() * Utility::skewSymmetric(pts_imu_j);

            jacobian_pose_j.leftCols<6>() = reduce * jaco_j;
            jacobian_pose_j.rightCols<1>().setZero();
        }
        if (jacobians[2])//外参数
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex_pose(jacobians[2]);
            Eigen::Matrix<double, 3, 6> jaco_ex;
            jaco_ex.leftCols<3>() = ric.transpose() * (Rj.transpose() * Rc_i_ * Ri - Eigen::Matrix3d::Identity());
            Eigen::Matrix3d tmp_r = ric.transpose() * Rj.transpose() * Rc_i_ * Ri * ric;
            jaco_ex.rightCols<3>() = -tmp_r * Utility::skewSymmetric(pts_camera_i) + Utility::skewSymmetric(tmp_r * pts_camera_i) +
                                     Utility::skewSymmetric(ric.transpose() * (Rj.transpose() * (Ri * tic + Pi - Pj) - tic));
            jacobian_ex_pose.leftCols<6>() = reduce * jaco_ex;
            jacobian_ex_pose.rightCols<1>().setZero();
        }
        if (jacobians[3])//逆深度
        {
            Eigen::Map<Eigen::Vector2d> jacobian_feature(jacobians[3]);
            jacobian_feature = reduce * ric.transpose() * Rj.transpose() * Rc_i_ * Ri * ric * pts_i_td * -1.0 / (inv_dep_i * inv_dep_i);
        }
        if (jacobians[4])
        {
            Eigen::Map<Eigen::Vector2d> jacobian_td(jacobians[4]);
            jacobian_td = reduce * ric.transpose() * Rj.transpose() * Rc_i_ * Ri * ric * velocity_i / inv_dep_i * -1.0  +
                          sqrt_info * velocity_j.head(2);
        }
    }
    sum_t += tic_toc.toc();

    return true;
}

void ProjectionTdFactor::check(double **parameters)
{
    double *res = new double[2];
    double **jaco = new double *[5];
    jaco[0] = new double[2 * 7];
    jaco[1] = new double[2 * 7];
    jaco[2] = new double[2 * 7];
    jaco[3] = new double[2 * 1];
    jaco[4] = new double[2 * 1];
    Evaluate(parameters, res, jaco);
    puts("check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Vector2d>(jaco[3]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Vector2d>(jaco[4]) << std::endl
              << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    double inv_dep_i = parameters[3][0];
    double td = parameters[4][0];

    //Eigen::Quaterniond Qij_c = (Qj * qic).inverse() * (Qi * qic);
    //Eigen::Quaterniond Qc_i_(1, Qij_c.x()/(td_j-td_i)*td, Qij_c.y()/(td_j-td_i)*td, Qij_c.z()/(td_j-td_i)*td);

    Eigen::Vector3d pts_i_td, pts_j_td;
    pts_i_td = pts_i - (td - td_i + TR / ROW * row_i) * velocity_i;
    pts_j_td = pts_j - (td - td_j + TR / ROW * row_j) * velocity_j;
    Eigen::Vector3d pts_camera_i = pts_i_td / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    //Eigen::Vector3d pts_w = (Qi * qic) * Qc_i_ * pts_camera_i + Qi * tic + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    Eigen::Vector2d residual;

#ifdef UNIT_SPHERE_ERROR 
    residual =  tangent_base * (pts_camera_j.normalized() - pts_j_td.normalized());
#else
    double dep_j = pts_camera_j.z();
    residual = (pts_camera_j / dep_j).head<2>() - pts_j_td.head<2>();
#endif
    residual = sqrt_info * residual;

    puts("num");
    std::cout << residual.transpose() << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 2, 20> num_jacobian;
    for (int k = 0; k < 20; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
        double inv_dep_i = parameters[3][0];
        double td = parameters[4][0];


        int a = k / 3, b = k % 3;
        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

        if (a == 0)
            Pi += delta;
        else if (a == 1)
            Qi = Qi * Utility::deltaQ(delta);
        else if (a == 2)
            Pj += delta;
        else if (a == 3)
            Qj = Qj * Utility::deltaQ(delta);
        else if (a == 4)
            tic += delta;
        else if (a == 5)
            qic = qic * Utility::deltaQ(delta);
        else if (a == 6 && b == 0)
            inv_dep_i += delta.x();
        else if (a == 6 && b == 1)
            td += delta.y();

        //Eigen::Quaterniond Qij_c = (Qj * qic).inverse() * (Qi * qic);
        //Eigen::Quaterniond Qc_i_(1, Qij_c.x()/(td_j-td_i)*td, Qij_c.y()/(td_j-td_i)*td, Qij_c.z()/(td_j-td_i)*td);

        Eigen::Vector3d pts_i_td, pts_j_td;
        pts_i_td = pts_i - (td - td_i + TR / ROW * row_i) * velocity_i;
        pts_j_td = pts_j - (td - td_j + TR / ROW * row_j) * velocity_j;
        Eigen::Vector3d pts_camera_i = pts_i_td / inv_dep_i;
        Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
        Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
        Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
        Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);
        Eigen::Vector2d tmp_residual;

#ifdef UNIT_SPHERE_ERROR 
        tmp_residual =  tangent_base * (pts_camera_j.normalized() - pts_j_td.normalized());
#else
        double dep_j = pts_camera_j.z();
        tmp_residual = (pts_camera_j / dep_j).head<2>() - pts_j_td.head<2>();
#endif
        tmp_residual = sqrt_info * tmp_residual;

        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}
