#include "AC_LQRControl.h"




void AC_LQRControl::get_A_matrix(const Matrix<float>& x, const Matrix<float>& u)
{
    float wx = u.get(0, 0);
    float wy = u.get(0, 1);
    float wz = u.get(0, 2);
    float norm_thrust = u.get(0, 3);
    // float wx = _ahrs.get_gyro().x;
    // float wy = _ahrs.get_gyro().y;
    // float wz = _ahrs.get_gyro().z;
    // float norm_thrust = -AP::ins().get_accel().z;
    // printf("accel : %.4f \n", norm_thrust);
    Quaternion q{x.get(0, 3), x.get(0, 4), x.get(0, 5), x.get(0, 6)};
    Matrix<float> q_partial_correction{4};
    Matrix<float> dqdot_dq{4};
    Matrix<float> dvdot_dq{3, 4};
    Matrix<float> q_vec{4, 1};

    q_vec.set(0, 0, q[0]);
    q_vec.set(1, 0, q[1]);
    q_vec.set(2, 0, q[2]);
    q_vec.set(3, 0, q[3]);

    A.zero();

    //Position
    A.set(0, 7, 1);
    A.set(1, 8, 1);
    A.set(2, 9, 1);
    
    Matrix<float> Identity{4};

    //Orientation
    q_partial_correction = powf(q.length(), -1.0)*(Identity.identity() - powf(q.length(), -2.0)*(q_vec * q_vec.transpose()));

    dqdot_dq << 0, -wx, -wy, -wz,
                wx, 0, wz, -wy,
                wy, -wz, 0, wx,
                wz, wy, -wx, 0;
    dqdot_dq = 0.5 * dqdot_dq * q_partial_correction;

    A.set(3, 3, dqdot_dq.get(0, 0));
    A.set(3, 4, dqdot_dq.get(0, 1));
    A.set(3, 5, dqdot_dq.get(0, 2));
    A.set(3, 6, dqdot_dq.get(0, 3));

    A.set(4, 3,dqdot_dq.get(1, 0));
    A.set(4, 4,dqdot_dq.get(1, 1));
    A.set(4, 5,dqdot_dq.get(1, 2));
    A.set(4, 6,dqdot_dq.get(1, 3));

    A.set(5, 3,dqdot_dq.get(2, 0));
    A.set(5, 4,dqdot_dq.get(2, 1));
    A.set(5, 5,dqdot_dq.get(2, 2));
    A.set(5, 6,dqdot_dq.get(2, 3));

    A.set(6, 3,dqdot_dq.get(3, 0));
    A.set(6, 4,dqdot_dq.get(3, 1));
    A.set(6, 5,dqdot_dq.get(3, 2));
    A.set(6, 6,dqdot_dq.get(3, 3));


    //Velocity
    dvdot_dq << q[2],  q[3],  q[0], q[1],
               -q[1], -q[0],  q[3], q[2],
                q[0], -q[1], -q[2], q[3];

    dvdot_dq = 2*norm_thrust*dvdot_dq*q_partial_correction;

    A.set(7, 3, dvdot_dq.get(0, 0));
    A.set(7, 4, dvdot_dq.get(0, 1));
    A.set(7, 5, dvdot_dq.get(0, 2));
    A.set(7, 6, dvdot_dq.get(0, 3));

    A.set(8, 3, dvdot_dq.get(1, 0));
    A.set(8, 4, dvdot_dq.get(1, 1));
    A.set(8, 5, dvdot_dq.get(1, 2));
    A.set(8, 6, dvdot_dq.get(1, 3));

    A.set(9, 3, dvdot_dq.get(2, 0));
    A.set(9, 4, dvdot_dq.get(2, 1));
    A.set(9, 5, dvdot_dq.get(2, 2));
    A.set(9, 6, dvdot_dq.get(2, 3));



}

void AC_LQRControl::get_B_matrix(const Matrix<float>& x, const Matrix<float>& u)
{
   Quaternion q{x.get(0, 3), x.get(0, 4), x.get(0, 5), x.get(0, 6)};
   Matrix<float> dvdot_dc{3, 1};
   Matrix<float> dqdot_dw{4, 3};

   B.zero();

   dvdot_dc << (q[0]*q[2] + q[1]*q[3]),
               (q[2]*q[3] - q[0]*q[1]),
               powf(q[0],2) - powf(q[1],2) - powf(q[2],2) + powf(q[3],2);

   B.set(7, 3, dvdot_dc.get(0, 0));
   B.set(8, 3, dvdot_dc.get(1, 0));
   B.set(9, 3, dvdot_dc.get(2, 0));

   dqdot_dw << -q[1], -q[2], -q[3],
                q[0], -q[3], -q[2],
                q[3],  q[0],  q[1],
               -q[2],  q[1],  q[0];

   dqdot_dw = 0.5 * dqdot_dw;

   B.set(3, 0, dqdot_dw.get(0, 0));
   B.set(3, 1, dqdot_dw.get(0, 1));
   B.set(3, 2, dqdot_dw.get(0, 2));

   B.set(4, 0, dqdot_dw.get(1, 0));
   B.set(4, 1, dqdot_dw.get(1, 1));
   B.set(4, 2, dqdot_dw.get(1, 2));

   B.set(5, 0, dqdot_dw.get(2, 0));
   B.set(5, 1, dqdot_dw.get(2, 1));
   B.set(5, 2, dqdot_dw.get(2, 2));

   B.set(6, 0, dqdot_dw.get(3, 0));
   B.set(6, 1, dqdot_dw.get(3, 1));
   B.set(6, 2, dqdot_dw.get(3, 2));

}

bool AC_LQRControl::solveRiccatiIterationC(const Matrix<float>& _A, const Matrix<float>& _B,
                            const Matrix<float>& _Q, const Matrix<float>& _R,
                            Matrix<float> &_P, const double dt,
                            const double &tolerance,
                            const uint32_t iter_max) {
  _P = _Q; // initialize

  Matrix<float> P_next{10};

  Matrix<float> AT = _A.transpose();
  Matrix<float> BT = _B.transpose();
  Matrix<float> Rinv = _R.inv();

  double diff;
  for (uint32_t i = 0; i < iter_max; ++i) {
    P_next = _P + (_P * _A + AT * _P - _P * _B * Rinv * BT * _P + _Q) * dt;
    diff = fabs((P_next - _P).maxCoeff());
    _P = P_next;
    if (diff < tolerance) {
      return true;
    }
  }
  return false; // over iteration limit
}

void AC_LQRControl::linearize()
{
    Vector3f curr_pos = _inav.get_position() / 100.0f;
    Vector3f curr_vel = _inav.get_velocity() / 100.0f;
    // _ahrs.get_quat_body_to_ned(_attitude_quat);
    _attitude_quat.from_euler(_ahrs.roll, _ahrs.pitch,  _ahrs.yaw);
    Matrix<float> x{1, 10};
    // printf("pos: %.2f, %.2f, %.2f \n", curr_pos.x, curr_pos.y, curr_pos.z);
    x << curr_pos.x, curr_pos.y, curr_pos.z,
        _attitude_quat[0], _attitude_quat[1], _attitude_quat[2], _attitude_quat[3],
        curr_vel.x, curr_vel.y, curr_vel.z;

    
    
    get_A_matrix(x, u_ref);
    get_B_matrix(x, u_ref);
    // printf("A:\n");
    // A.printmat();
    // printf(solveRiccatiIterationC(A, B, Q, R, P,0.001, 1.E-4, 1000000) ? "true \n" : "false \n");
    // solveRiccatiIterationC(A, B, Q, R, P,0.001, 1.E-4, 1000000);
    solveRiccatiIterationC(A, B, Q, R, P,0.001, 1.E-3, 100000);
    // solveRiccatiIterationC(A, B, Q, R, P,0.001, 1.E-1, 1000);
    // solveRiccatiIterationC(A, B, Q, R, P,0.01, 1.E-3, 10000);
    // solveRiccatiIterationC(A, B, Q, R, P,0.001, 1.E-2, 10000);
    Matrix<float> K_new{4, 10};
    // A.printmat();
    // B.printmat();

    // printf("p:\n");
    // P.printmat();
    // printf("BT:\n");
    // (B.transpose()).printmat();
    // printf("BT*P:\n");
    // (B.transpose() * P).printmat();
    K_new = -1 * R.inv() * B.transpose() * P;

    // sock.sendto(&data_to_send, sizeof(data_to_send), "127.0.0.1", 9002);
    // printf("Knew:\n");
    // K_new.printmat();
    // K_new.printmat();
    // printf( K_new.hasNaN() ? "true" : "false");
    // if( !K_new.hasNaN() )
    K = K_new;
    // if (t == false)
    // {
    //   K = K_new;
    //   t = true;
    // }

    
    
    // printf("R : %.2f, %.2f, %.2f, %.2f,\n", R.inv().get(0, 0), R.inv().get(1, 1), R.inv().get(2, 2), R.inv().get(3, 3));
    // printf("Q : %.2f, %.2f, %.2f, %.2f,\n", Q.get(0, 0), Q.get(1, 1), Q.get(2, 2), Q.get(3, 3));
    // Matrix<float> v {3,2};
    // Matrix<float> n {3,3};
    // v << 2, 3,
    // 0, 5,
    // 5, 2;
    // n << 3, 0, 1,
    //  0, 1, 0, 
    //  18, 3, 0;

}

void AC_LQRControl::set_Q_matrix()
{
  Q.set(0, 0, 1.0f);
  Q.set(1, 1, 1.0f);


  // Q.set(1, 2, 100.0f);
  // Q.set(2, 2, 30000.0f);
  // Q.set(2, 2, 1000000.0f);
  Q.set(2, 2, 10.0f);
  Q.set(3, 3, 0.1f);
  // Q.set(4, 4, 4.0f);
  // Q.set(5, 5, 4.0f);
  // Q.set(6, 6, 4.0f);
  Q.set(4, 4, 0.1f);
  Q.set(5, 5, 0.1f);
  Q.set(6, 6, 0.1f);
 
  // Q.set(7, 7, 1000.0f);
  // Q.set(8, 8, 0.0f);

  Q.set(7, 7, 20.0f);
  Q.set(8, 8, 20.0f);
  Q.set(9, 9, 20.0f);
}

void AC_LQRControl::set_R_matrix()
{
  R << 5.0f, 0, 0, 0,
       0, 5.0f, 0, 0,
       0, 0, 15.0f ,0,
      //  0, 0, 0, 0.05f;
      0, 0, 0, 0.1f;
}

void AC_LQRControl::set_output()
{
  Vector3f curr_pos = _inav.get_position() / 100.0f;
  Vector3f curr_vel = _inav.get_velocity() / 100.0f;
  _ahrs.get_quat_body_to_ned(_attitude_quat);
  
  Matrix<float> x{1, 10};
  Matrix<float> x_error{1, 10};

  // printf("pos: %.2f, %.2f, %.2f \n", curr_pos.x, curr_pos.y, curr_pos.z);
  x << curr_pos.x, curr_pos.y, curr_pos.z,
      _attitude_quat[0], _attitude_quat[1], _attitude_quat[2], _attitude_quat[3],
      curr_vel.x, curr_vel.y, curr_vel.z;

  x_error = x - x_ref;
  x_error.set(0, 0, x_ref.get(0, 0) - x.get(0, 0));
  x_error.set(0, 1, x_ref.get(0, 1) - x.get(0, 1));
  x_error.set(0, 7, x_ref.get(0, 7) - x.get(0, 7));
  x_error.set(0, 8, x_ref.get(0, 8) - x.get(0, 8));
  Quaternion q{x.get(0, 3), x.get(0, 4), x.get(0, 5), x.get(0, 6)};
  Quaternion qref{x_ref.get(0, 3), x_ref.get(0, 4), x_ref.get(0, 5), x_ref.get(0, 6)};
  Quaternion qerror =qref.inverse() * q;
  // Quaternion qerror;
  // float roll_ref, pitch_ref, yaw_ref;
  // error.to_euler(roll_ref, pitch_ref, yaw_ref);
  // printf("%.2f, %.2f , %.2f \n", degrees(roll_ref), degrees(pitch_ref), degrees(yaw_ref));
  // printf("%.2f, %.2f , %.2f \n", degrees(_ahrs.roll), degrees(_ahrs.pitch), degrees(_ahrs.yaw));
  // qerror.from_euler(_ahrs.roll - roll_ref, _ahrs.pitch - pitch_ref, _ahrs.yaw - yaw_ref);

  // Quaternion qerror;
  // qerror.from_euler(_ahrs.roll - radians(0), _ahrs.pitch - radians(0),  _ahrs.yaw- radians(0));
  // qerror.to_euler(roll_ref, pitch_ref, yaw_ref);
  // printf("%.2f, %.2f , %.2f \n", degrees(roll_ref), degrees(pitch_ref), degrees(yaw_ref));

  x_error.set(0, 3, qerror[0]);
  x_error.set(0, 4, qerror[1]);
  x_error.set(0, 5, qerror[2]);
  x_error.set(0, 6, qerror[3]);
  
  AP::logger().Write("LQR", "TimeUS,X,Y,Z,TX,TY,TZ,VX,VY,VZ,TVX,TVY,TVZ",
                   "smmmmmmnnnnnn", // units: seconds, meters
                   "F000000000000", // mult: 1e-6, 1e-2
                   "Qffffffffffff", // format: uint64_t, float
                   AP_HAL::micros64(),
                   curr_pos.x,
                   curr_pos.y,
                   curr_pos.z,
                   x_ref.get(0, 0),
                   x_ref.get(0, 1),
                   x_ref.get(0, 2),
                   curr_vel.x,
                   curr_vel.y,
                   curr_vel.z,
                   x_ref.get(0, 7),
                   x_ref.get(0, 8),
                   x_ref.get(0, 9));
  Matrix<float> out {4, 1};
  // qref.to_euler(roll_ref, pitch_ref, yaw_ref);
  // printf("eref :%.4f , %.4f , %.4f , %.4f\n", error[0], error[1], error[2], error[3]);
  // printf("error:\n");
  // x_error.printmat();
  // K.printmat();
  out = u_ref.transpose() + K * x_error.transpose();
  // (K * x_error.transpose()).printmat();
  K.printmat();
  // out.printmat();
  float rate_x = constrain_float(out.get(0, 0), -2.0f, 2.0f);
  _attitude_control._rate_target_ang_vel.x = rate_x;
  // _attitude_control._rate_target_ang_vel.x = out.get(0, 0); 

  float rate_y = constrain_float(out.get(1, 0), -2.0f, 2.0f);
  _attitude_control._rate_target_ang_vel.y = rate_y;
  // _attitude_control._rate_target_ang_vel.y = out.get(1, 0);

  float rate_z = constrain_float(out.get(2, 0), -2.0f, 2.0f);
  _attitude_control._rate_target_ang_vel.z = rate_z;
  // _attitude_control._rate_target_ang_vel.z = out.get(2, 0);


  // printf("wx: %.2f \n",out.get(0, 0));
  // printf("wy: %.2f \n",out.get(1, 0));
  // printf("wz: %.2f \n",out.get(2, 0));
  // printf("throttle : %.2f \n",out.get(3, 0));
  // out.printmat();
  // float thr = out.get(3, 0)/ 8 - (5/4);
  // thr = constrain_float(thr, -1.0f, 1.0f);
  // printf("thr : %.4f \n", _motors.get_throttle_hover());
  float j = 9.8 / _motors.get_throttle_hover();
  // _attitude_control.set_throttle_out((thr/16) - 1/8 , true, 2.0f);
  float thr = constrain_float(out.get(3, 0)/j, 0.0f, 1.0f);
  _attitude_control.set_throttle_out(thr, false, 2.0f);
  // _motors.set_throttle(thr);

}

void AC_LQRControl::set_refernce(double x, double y, double z)
{
  
  x_ref.set(0, 0, x);
  x_ref.set(0, 1, y);
  x_ref.set(0, 2, z);
  
  Quaternion q{};
  q.from_euler(radians(0.0f), radians(0.0f), radians(0.0f));
  
  x_ref.set(0, 3 , q[0]);
  x_ref.set(0, 4 , q[1]);
  x_ref.set(0, 5 , q[2]);
  x_ref.set(0, 6 , q[3]);

  x_ref.set(0, 7, 0);
  x_ref.set(0, 8, 0);
  x_ref.set(0, 9, 0);

  u_ref.set(0, 0, 0);
  u_ref.set(0, 1, 0);
  u_ref.set(0, 2, 0);
  u_ref.set(0, 3, 9.8);
}