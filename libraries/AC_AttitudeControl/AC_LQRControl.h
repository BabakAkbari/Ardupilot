#pragma once
#define ALLOW_DOUBLE_MATH_FUNCTIONS

/// @file    AC_LQRControl.h
/// @brief   ArduCopter LQR control library


#include <AP_Common/AP_Common.h>
#include <AP_Param/AP_Param.h>
#include <AP_Math/AP_Math.h>
#include <AP_Math/Matrix.h>
#include <AP_Vehicle/AP_Vehicle.h>
#include <AP_AHRS/AP_AHRS_View.h>
#include <AP_Motors/AP_Motors.h>
#include <AP_InertialNav/AP_InertialNav.h>     // Inertial Navigation library
#include "AC_AttitudeControl.h" // Attitude control library
#include <AP_Logger/AP_Logger.h>


class AC_LQRControl {
    public:
    AC_LQRControl( AP_AHRS_View &ahrs,
                        AC_AttitudeControl& attitude_control,
                        const AP_InertialNav& inav,
                        AP_Motors& motors,
                        float dt) :
        _dt(dt),
        _ahrs(ahrs),
        _inav(inav),
        // _aparm(aparm),
        _motors(motors),
        _attitude_control(attitude_control)
        {
            // AP_Param::setup_object_defaults(this, var_info);
            set_Q_matrix();
            set_R_matrix();
        }
    // empty destructor to suppress compiler warning
	// virtual ~AC_LQRControl() {}
    void get_A_matrix(const Matrix<float>& x, const Matrix<float>& u);
    void get_B_matrix(const Matrix<float>& x, const Matrix<float>& u);
    // void lqr(Matrix<float> _Q, Matrix<float> _R, Matrix<float> _A, Matrix<float> _B, Matrix<float> _K);
    bool solveRiccatiIterationC(const Matrix<float>& _A, const Matrix<float>& _B,
                            const Matrix<float>& _Q, const Matrix<float>& _R,
                            Matrix<float> &_P, const double dt = 0.001,
                            const double &tolerance = 1.E-5,
                            const uint32_t iter_max = 100000);
    void linearize();

    void set_Q_matrix();
    void set_R_matrix();

    void set_output();
    void set_refernce(double x, double y, double z);
    // void Write(const char *name, const char *labels, const char *units, const char *mults, const char *fmt, ...);


    // This represents the angular velocity in radians per second in the body frame, used in the angular
    // velocity controller.
    Vector3f            _lqr_rate_target_ang_vel;
    
    Matrix<float> A{10};
    Matrix<float> B{10, 4};
    Matrix<float> Q{10};
    Matrix<float> R{4};
    Matrix<float> P{10};

    Matrix<float> x_ref{1,10};
    Matrix<float> u_ref{1,4};
    Matrix<float> K{4, 10};
    Matrix<float> out {4, 1};
    bool t = false;
    
    // Intersampling period in seconds
    float               _dt;
    size_t e = 0;
    // References to external libraries
    const AP_InertialNav&       _inav;
    const AP_AHRS_View&  _ahrs;
    // const AP_Vehicle::MultiCopter &_aparm;
    AP_Motors&          _motors;
    AC_AttitudeControl&         _attitude_control;
    Quaternion          _attitude_quat;

    
    // User settable parameters
    // static const struct AP_Param::GroupInfo var_info[];

};