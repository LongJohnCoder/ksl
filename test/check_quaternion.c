#include <config.h>
#include <stdlib.h>

#include <assert.h>
#include <check.h>
#include <float.h>
#include <stdio.h>
#include <time.h>

//#include "linalg.h"
#include "matrix.h"
#include "quaternion.h"

START_TEST(test_quaternion_create) {
  double sf = sqrt(1.2);
  ksl_quaternion_t q = ksl_quaternion(0.2 / sf, 0.4 / sf, 0.6 / sf, 0.8 / sf);
  ck_assert(q.x == q.at[0]);
  ck_assert(q.y == q.at[1]);
  ck_assert(q.z == q.at[2]);
  ck_assert(q.w == q.at[3]);
  ck_assert(q.x == q.r.x);
  ck_assert(q.y == q.r.y);
  ck_assert(q.z == q.r.z);
  ck_assert_double_eq_tol(q.x, 0.2 / sf, 1e-9);
  ck_assert_double_eq_tol(q.y, 0.4 / sf, 1e-9);
  ck_assert_double_eq_tol(q.z, 0.6 / sf, 1e-9);
  ck_assert_double_eq_tol(q.w, 0.8 / sf, 1e-9);
}
END_TEST

START_TEST(test_quaternionf_create) {
  float sf = sqrt(1.2);
  ksl_quaternionf_t q = ksl_quaternionf(0.2 / sf, 0.4 / sf, 0.6 / sf, 0.8 / sf);
  ck_assert(q.x == q.at[0]);
  ck_assert(q.y == q.at[1]);
  ck_assert(q.z == q.at[2]);
  ck_assert(q.w == q.at[3]);
  ck_assert(q.x == q.r.x);
  ck_assert(q.y == q.r.y);
  ck_assert(q.z == q.r.z);
  ck_assert_float_eq_tol(q.x, 0.2 / sf, 1e-5);
  ck_assert_float_eq_tol(q.y, 0.4 / sf, 1e-5);
  ck_assert_float_eq_tol(q.z, 0.6 / sf, 1e-5);
  ck_assert_float_eq_tol(q.w, 0.8 / sf, 1e-5);
}
END_TEST

START_TEST(test_quaternion_unit) {
  double sf = sqrt(1.2);
  ksl_quaternion_t q = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ck_assert_double_eq_tol(q.x, 0.2 / sf, 1e-9);
  ck_assert_double_eq_tol(q.y, 0.4 / sf, 1e-9);
  ck_assert_double_eq_tol(q.z, 0.6 / sf, 1e-9);
  ck_assert_double_eq_tol(q.w, 0.8 / sf, 1e-9);
}
END_TEST

START_TEST(test_quaternionf_unit) {
  float sf = sqrt(1.2);
  ksl_quaternionf_t q = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ck_assert_float_eq_tol(q.x, 0.2 / sf, 1e-5);
  ck_assert_float_eq_tol(q.y, 0.4 / sf, 1e-5);
  ck_assert_float_eq_tol(q.z, 0.6 / sf, 1e-5);
  ck_assert_float_eq_tol(q.w, 0.8 / sf, 1e-5);
}
END_TEST

START_TEST(test_quaternion_alloc) {
  ksl_quaternion_t* q3 = ksl_quaternion_alloc(3);
  ksl_quaternion_t q = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  for(int i = 0; i < 3; i++) {
    q3[i] = q;
    ck_assert_double_eq_tol(q3[i].x, q.x, 1e-9);
    ck_assert_double_eq_tol(q3[i].y, q.y, 1e-9);
    ck_assert_double_eq_tol(q3[i].z, q.z, 1e-9);
    ck_assert_double_eq_tol(q3[i].w, q.w, 1e-9);
  }
}
END_TEST

START_TEST(test_quaternionf_alloc) {
  ksl_quaternionf_t* q3 = ksl_quaternionf_alloc(3);
  ksl_quaternionf_t q = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  for(int i = 0; i < 3; i++) {
    q3[i] = q;
    ck_assert_float_eq_tol(q3[i].x, q.x, 1e-5);
    ck_assert_float_eq_tol(q3[i].y, q.y, 1e-5);
    ck_assert_float_eq_tol(q3[i].z, q.z, 1e-5);
    ck_assert_float_eq_tol(q3[i].w, q.w, 1e-5);
  }
}
END_TEST

START_TEST(test_quaternion_maxIndex) {
  ksl_quaternion_t q = ksl_quaternion(0.2, -0.2, -0.6, 0.6);
  int i = ksl_quaternion_maxIndex(&q);
  ck_assert(i == 2);
  q = ksl_quaternion(0.6, -0.2, -0.6, 0.2);
  i = ksl_quaternion_maxIndex(&q);
  ck_assert(i == 0);
}
END_TEST

START_TEST(test_quaternionf_maxIndex) {
  ksl_quaternionf_t q = ksl_quaternionf(0.2, -0.2, -0.6, 0.6);
  int i = ksl_quaternionf_maxIndex(&q);
  ck_assert(i == 2);
  q = ksl_quaternionf(0.6, -0.2, -0.6, 0.2);
  i = ksl_quaternionf_maxIndex(&q);
  ck_assert(i == 0);
}
END_TEST

START_TEST(test_dot_qq) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2 = ksl_quaternion_unit(0.8, 0.6, 0.4, 0.2);
  ck_assert_double_eq_tol(ksl_dot_qq(&q1, &q2), 2.0 / 3.0, 1e-9);
}
END_TEST

START_TEST(test_dot_qqf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2 = ksl_quaternionf_unit(0.8, 0.6, 0.4, 0.2);
  ck_assert_float_eq_tol(ksl_dot_qqf(&q1, &q2), 2.0 / 3.0, 1e-5);
}
END_TEST

START_TEST(test_quaternion_normalize) {
  double sf = sqrt(1.2);
  ksl_quaternion_t q = ksl_quaternion(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_normalize(&q);
  ck_assert_double_eq_tol(q.x, 0.2 / sf, 1e-9);
  ck_assert_double_eq_tol(q.y, 0.4 / sf, 1e-9);
  ck_assert_double_eq_tol(q.z, 0.6 / sf, 1e-9);
  ck_assert_double_eq_tol(q.w, 0.8 / sf, 1e-9);
}
END_TEST

START_TEST(test_quaternionf_normalize) {
  double sf = sqrt(1.2);
  ksl_quaternionf_t q = ksl_quaternionf(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_normalize(&q);
  ck_assert_float_eq_tol(q.x, 0.2 / sf, 1e-5);
  ck_assert_float_eq_tol(q.y, 0.4 / sf, 1e-5);
  ck_assert_float_eq_tol(q.z, 0.6 / sf, 1e-5);
  ck_assert_float_eq_tol(q.w, 0.8 / sf, 1e-5);
}
END_TEST

START_TEST(test_axpy_qq) {
  double a = 3.0;
  ksl_quaternion_t q1 = ksl_quaternion(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2 = ksl_quaternion(0.7, 0.5, 0.3, 0.1);
  ksl_axpy_qq(a, &q1, &q2);
  ck_assert_double_eq_tol(q2.x, 1.3, 1e-9);
  ck_assert_double_eq_tol(q2.y, 1.7, 1e-9);
  ck_assert_double_eq_tol(q2.z, 2.1, 1e-9);
  ck_assert_double_eq_tol(q2.w, 2.5, 1e-9);
}
END_TEST

START_TEST(test_axpy_qqf) {
  float a = 3.0;
  ksl_quaternionf_t q1 = ksl_quaternionf(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2 = ksl_quaternionf(0.7, 0.5, 0.3, 0.1);
  ksl_axpy_qqf(a, &q1, &q2);
  ck_assert_float_eq_tol(q2.x, 1.3, 1e-5);
  ck_assert_float_eq_tol(q2.y, 1.7, 1e-5);
  ck_assert_float_eq_tol(q2.z, 2.1, 1e-5);
  ck_assert_float_eq_tol(q2.w, 2.5, 1e-5);
}
END_TEST

START_TEST(test_xpy_qq) {
  ksl_quaternion_t q1 = ksl_quaternion(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2 = ksl_quaternion(0.7, 0.5, 0.3, 0.1);
  ksl_xpy_qq(&q1, &q2);
  ck_assert_double_eq_tol(q2.x, 0.9, 1e-9);
  ck_assert_double_eq_tol(q2.y, 0.9, 1e-9);
  ck_assert_double_eq_tol(q2.z, 0.9, 1e-9);
  ck_assert_double_eq_tol(q2.w, 0.9, 1e-9);
}
END_TEST

START_TEST(test_xpy_qqf) {
  ksl_quaternionf_t q1 = ksl_quaternionf(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2 = ksl_quaternionf(0.7, 0.5, 0.3, 0.1);
  ksl_xpy_qqf(&q1, &q2);
  ck_assert_float_eq_tol(q2.x, 0.9, 1e-5);
  ck_assert_float_eq_tol(q2.y, 0.9, 1e-5);
  ck_assert_float_eq_tol(q2.z, 0.9, 1e-5);
  ck_assert_float_eq_tol(q2.w, 0.9, 1e-5);
}
END_TEST

START_TEST(test_nxpy_qq) {
  ksl_quaternion_t q1 = ksl_quaternion(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2 = ksl_quaternion(0.7, 0.5, 0.3, 0.1);
  ksl_nxpy_qq(&q1, &q2);
  ck_assert_double_eq_tol(q2.x, 0.5, 1e-9);
  ck_assert_double_eq_tol(q2.y, 0.1, 1e-9);
  ck_assert_double_eq_tol(q2.z, -0.3, 1e-9);
  ck_assert_double_eq_tol(q2.w, -0.7, 1e-9);
}
END_TEST

START_TEST(test_nxpy_qqf) {
  ksl_quaternionf_t q1 = ksl_quaternionf(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2 = ksl_quaternionf(0.7, 0.5, 0.3, 0.1);
  ksl_nxpy_qqf(&q1, &q2);
  ck_assert_float_eq_tol(q2.x, 0.5, 1e-5);
  ck_assert_float_eq_tol(q2.y, 0.1, 1e-5);
  ck_assert_float_eq_tol(q2.z, -0.3, 1e-5);
  ck_assert_float_eq_tol(q2.w, -0.7, 1e-5);
}
END_TEST

START_TEST(test_mat3x3_to_quaternion) {
  ksl_mat3x3_t* R = ksl_mat3x3_alloc(24);
  R[0] = ksl_mat3x3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
  R[1] = ksl_mat3x3(1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
  R[2] = ksl_mat3x3(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0);
  R[3] = ksl_mat3x3(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0);
  R[4] = ksl_mat3x3(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0);
  R[5] = ksl_mat3x3(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
  R[6] = ksl_mat3x3(-1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0);
  R[7] = ksl_mat3x3(-1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0);
  R[8] = ksl_mat3x3(0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0);
  R[9] = ksl_mat3x3(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  R[10] = ksl_mat3x3(0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
  R[11] = ksl_mat3x3(0.0, 1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0);
  R[12] = ksl_mat3x3(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  R[13] = ksl_mat3x3(0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0);
  R[14] = ksl_mat3x3(0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0);
  R[15] = ksl_mat3x3(0.0, -1.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0);
  R[16] = ksl_mat3x3(0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0);
  R[17] = ksl_mat3x3(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0);
  R[18] = ksl_mat3x3(0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0);
  R[19] = ksl_mat3x3(0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  R[20] = ksl_mat3x3(0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0);
  R[21] = ksl_mat3x3(0.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0);
  R[22] = ksl_mat3x3(0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0);
  R[23] = ksl_mat3x3(0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  for(int i = 0; i < 24; i++) {
    double det = 0.0;
    det += R[i].m00 * (R[i].m11 * R[i].m22 - R[i].m21 * R[i].m12);
    det -= R[i].m01 * (R[i].m10 * R[i].m22 - R[i].m20 * R[i].m12);
    det += R[i].m02 * (R[i].m10 * R[i].m21 - R[i].m20 * R[i].m11);
    ck_assert_double_eq_tol(det, 1.0, 1e-9);
    ksl_quaternion_t q;
    ksl_mat3x3_to_quaternion(&R[i], &q);
    ksl_mat3x3_t R1;
    ksl_quaternion_to_mat3x3(&q, &R1);
    for(int j = 0; j < 9; j++) {
      ck_assert_double_eq_tol(R[i].as_array[j], R1.as_array[j], 1e-9);
    }
  }
  for(int i = 0; i < 100; i++) {
    ksl_quaternion_t q1;
    for(int j = 0; j < 4; j++) {
      q1.at[j] = 2.0 * (rand()) - 1.0;
    }
    ksl_quaternion_normalize(&q1);
    ksl_mat3x3_t R1;
    ksl_quaternion_to_mat3x3(&q1, &R1);
    ksl_quaternion_t q2;
    ksl_mat3x3_to_quaternion(&R1, &q2);
    int k = ksl_quaternion_maxIndex(&q1);
    if(fabs(q1.at[k] - q2.at[k]) < 1e-2) {
      ck_assert_double_eq_tol(q1.x, q2.x, 1e-9);
      ck_assert_double_eq_tol(q1.y, q2.y, 1e-9);
      ck_assert_double_eq_tol(q1.z, q2.z, 1e-9);
      ck_assert_double_eq_tol(q1.w, q2.w, 1e-9);
    } else {
      ck_assert_double_eq_tol(q1.x, -q2.x, 1e-9);
      ck_assert_double_eq_tol(q1.y, -q2.y, 1e-9);
      ck_assert_double_eq_tol(q1.z, -q2.z, 1e-9);
      ck_assert_double_eq_tol(q1.w, -q2.w, 1e-9);
    }
  }
}
END_TEST

START_TEST(test_mat3x3f_to_quaternionf) {
  ksl_mat3x3f_t* R = ksl_mat3x3f_alloc(24);
  R[0] = ksl_mat3x3f(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
  R[1] = ksl_mat3x3f(1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
  R[2] = ksl_mat3x3f(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0);
  R[3] = ksl_mat3x3f(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0);
  R[4] = ksl_mat3x3f(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0);
  R[5] = ksl_mat3x3f(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
  R[6] = ksl_mat3x3f(-1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0);
  R[7] = ksl_mat3x3f(-1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0);
  R[8] = ksl_mat3x3f(0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0);
  R[9] = ksl_mat3x3f(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  R[10] = ksl_mat3x3f(0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
  R[11] = ksl_mat3x3f(0.0, 1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0);
  R[12] = ksl_mat3x3f(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  R[13] = ksl_mat3x3f(0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0);
  R[14] = ksl_mat3x3f(0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0);
  R[15] = ksl_mat3x3f(0.0, -1.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0);
  R[16] = ksl_mat3x3f(0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0);
  R[17] = ksl_mat3x3f(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0);
  R[18] = ksl_mat3x3f(0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0);
  R[19] = ksl_mat3x3f(0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  R[20] = ksl_mat3x3f(0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0);
  R[21] = ksl_mat3x3f(0.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0);
  R[22] = ksl_mat3x3f(0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0);
  R[23] = ksl_mat3x3f(0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  for(int i = 0; i < 24; i++) {
    float det = 0.0;
    det += R[i].m00 * (R[i].m11 * R[i].m22 - R[i].m21 * R[i].m12);
    det -= R[i].m01 * (R[i].m10 * R[i].m22 - R[i].m20 * R[i].m12);
    det += R[i].m02 * (R[i].m10 * R[i].m21 - R[i].m20 * R[i].m11);
    ck_assert_float_eq_tol(det, 1.0, 1e-5);
    ksl_quaternionf_t q;
    ksl_mat3x3f_to_quaternionf(&R[i], &q);
    ksl_mat3x3f_t R1;
    ksl_quaternionf_to_mat3x3f(&q, &R1);
    for(int j = 0; j < 9; j++) {
      ck_assert_float_eq_tol(R[i].as_array[j], R1.as_array[j], 1e-5);
    }
  }
  for(int i = 0; i < 100; i++) {
    ksl_quaternionf_t q1;
    for(int j = 0; j < 4; j++) {
      q1.at[j] = 2.0 * (rand()) - 1.0;
    }
    ksl_quaternionf_normalize(&q1);
    ksl_mat3x3f_t R1;
    ksl_quaternionf_to_mat3x3f(&q1, &R1);
    ksl_quaternionf_t q2;
    ksl_mat3x3f_to_quaternionf(&R1, &q2);
    int k = ksl_quaternionf_maxIndex(&q1);
    if(fabs(q1.at[k] - q2.at[k]) < 1e-2) {
      ck_assert_float_eq_tol(q1.x, q2.x, 1e-5);
      ck_assert_float_eq_tol(q1.y, q2.y, 1e-5);
      ck_assert_float_eq_tol(q1.z, q2.z, 1e-5);
      ck_assert_float_eq_tol(q1.w, q2.w, 1e-5);
    } else {
      ck_assert_float_eq_tol(q1.x, -q2.x, 1e-5);
      ck_assert_float_eq_tol(q1.y, -q2.y, 1e-5);
      ck_assert_float_eq_tol(q1.z, -q2.z, 1e-5);
      ck_assert_float_eq_tol(q1.w, -q2.w, 1e-5);
    }
  }
}
END_TEST

START_TEST(test_quaternion_to_mat3x3) {
  for(int i = 0; i < 100; i++) {
    ksl_quaternion_t q1;
    for(int j = 0; j < 4; j++) {
      q1.at[j] = 2.0 * (rand()) - 1.0;
    }
    ksl_quaternion_normalize(&q1);
    ksl_mat3x3_t R1;
    ksl_quaternion_to_mat3x3(&q1, &R1);
    ksl_quaternion_t q2;
    ksl_mat3x3_to_quaternion(&R1, &q2);
    int k = ksl_quaternion_maxIndex(&q1);
    if(fabs(q1.at[k] - q2.at[k]) < 1e-2) {
      ck_assert_double_eq_tol(q1.x, q2.x, 1e-9);
      ck_assert_double_eq_tol(q1.y, q2.y, 1e-9);
      ck_assert_double_eq_tol(q1.z, q2.z, 1e-9);
      ck_assert_double_eq_tol(q1.w, q2.w, 1e-9);
    } else {
      ck_assert_double_eq_tol(q1.x, -q2.x, 1e-9);
      ck_assert_double_eq_tol(q1.y, -q2.y, 1e-9);
      ck_assert_double_eq_tol(q1.z, -q2.z, 1e-9);
      ck_assert_double_eq_tol(q1.w, -q2.w, 1e-9);
    }
  }
}
END_TEST

START_TEST(test_quaternionf_to_mat3x3f) {
  for(int i = 0; i < 100; i++) {
    ksl_quaternionf_t q1;
    for(int j = 0; j < 4; j++) {
      q1.at[j] = 2.0 * (rand()) - 1.0;
    }
    ksl_quaternionf_normalize(&q1);
    ksl_mat3x3f_t R1;
    ksl_quaternionf_to_mat3x3f(&q1, &R1);
    ksl_quaternionf_t q2;
    ksl_mat3x3f_to_quaternionf(&R1, &q2);
    int k = ksl_quaternionf_maxIndex(&q1);
    if(fabs(q1.at[k] - q2.at[k]) < 1e-2) {
      ck_assert_float_eq_tol(q1.x, q2.x, 1e-5);
      ck_assert_float_eq_tol(q1.y, q2.y, 1e-5);
      ck_assert_float_eq_tol(q1.z, q2.z, 1e-5);
      ck_assert_float_eq_tol(q1.w, q2.w, 1e-5);
    } else {
      ck_assert_float_eq_tol(q1.x, -q2.x, 1e-5);
      ck_assert_float_eq_tol(q1.y, -q2.y, 1e-5);
      ck_assert_float_eq_tol(q1.z, -q2.z, 1e-5);
      ck_assert_float_eq_tol(q1.w, -q2.w, 1e-5);
    }
  }
}
END_TEST

START_TEST(test_dpi_p_to_dq) {
  ksl_vec3_t dpi_p = ksl_vec3(1.0, 2.0, 3.0);
  ksl_quaternion_t q = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t dq;
  ksl_dpi_p_to_dq(&dpi_p, &q, &dq);
  ksl_quaternion_t dqv;
  dqv.x = 0.5 * (q.w * dpi_p.x + q.z * dpi_p.y - q.y * dpi_p.z);
  dqv.y = 0.5 * (-q.z * dpi_p.x + q.w * dpi_p.y + q.x * dpi_p.z);
  dqv.z = 0.5 * (q.y * dpi_p.x - q.x * dpi_p.y + q.w * dpi_p.z);
  dqv.w = 0.5 * (-q.x * dpi_p.x - q.y * dpi_p.y - q.z * dpi_p.z);
  ck_assert_double_eq_tol(dq.x, dqv.x, 1e-9);
  ck_assert_double_eq_tol(dq.y, dqv.y, 1e-9);
  ck_assert_double_eq_tol(dq.z, dqv.z, 1e-9);
  ck_assert_double_eq_tol(dq.w, dqv.w, 1e-9);
}
END_TEST

START_TEST(test_dpi_pf_to_dqf) {
  ksl_vec3f_t dpi_p = ksl_vec3f(1.0, 2.0, 3.0);
  ksl_quaternionf_t q = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t dq;
  ksl_dpi_pf_to_dqf(&dpi_p, &q, &dq);
  ksl_quaternionf_t dqv;
  dqv.x = 0.5 * (q.w * dpi_p.x + q.z * dpi_p.y - q.y * dpi_p.z);
  dqv.y = 0.5 * (-q.z * dpi_p.x + q.w * dpi_p.y + q.x * dpi_p.z);
  dqv.z = 0.5 * (q.y * dpi_p.x - q.x * dpi_p.y + q.w * dpi_p.z);
  dqv.w = 0.5 * (-q.x * dpi_p.x - q.y * dpi_p.y - q.z * dpi_p.z);
  ck_assert_float_eq_tol(dq.x, dqv.x, 1e-5);
  ck_assert_float_eq_tol(dq.y, dqv.y, 1e-5);
  ck_assert_float_eq_tol(dq.z, dqv.z, 1e-5);
  ck_assert_float_eq_tol(dq.w, dqv.w, 1e-5);
}
END_TEST

START_TEST(test_dpi_c_to_dq) {
  ksl_vec3_t dpi_p = ksl_vec3(1.0, 2.0, 3.0);
  ksl_quaternion_t q = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t dq;
  ksl_dpi_c_to_dq(&dpi_p, &q, &dq);
  ksl_quaternion_t dqv;
  dqv.x = 0.5 * (q.w * dpi_p.x - q.z * dpi_p.y + q.y * dpi_p.z);
  dqv.y = 0.5 * (q.z * dpi_p.x + q.w * dpi_p.y - q.x * dpi_p.z);
  dqv.z = 0.5 * (-q.y * dpi_p.x + q.x * dpi_p.y + q.w * dpi_p.z);
  dqv.w = 0.5 * (-q.x * dpi_p.x - q.y * dpi_p.y - q.z * dpi_p.z);
  ck_assert_double_eq_tol(dq.x, dqv.x, 1e-9);
  ck_assert_double_eq_tol(dq.y, dqv.y, 1e-9);
  ck_assert_double_eq_tol(dq.z, dqv.z, 1e-9);
  ck_assert_double_eq_tol(dq.w, dqv.w, 1e-9);
}
END_TEST

START_TEST(test_dpi_cf_to_dqf) {
  ksl_vec3f_t dpi_p = ksl_vec3f(1.0, 2.0, 3.0);
  ksl_quaternionf_t q = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t dq;
  ksl_dpi_cf_to_dqf(&dpi_p, &q, &dq);
  ksl_quaternionf_t dqv;
  dqv.x = 0.5 * (q.w * dpi_p.x - q.z * dpi_p.y + q.y * dpi_p.z);
  dqv.y = 0.5 * (q.z * dpi_p.x + q.w * dpi_p.y - q.x * dpi_p.z);
  dqv.z = 0.5 * (-q.y * dpi_p.x + q.x * dpi_p.y + q.w * dpi_p.z);
  dqv.w = 0.5 * (-q.x * dpi_p.x - q.y * dpi_p.y - q.z * dpi_p.z);
  ck_assert_float_eq_tol(dq.x, dqv.x, 1e-5);
  ck_assert_float_eq_tol(dq.y, dqv.y, 1e-5);
  ck_assert_float_eq_tol(dq.z, dqv.z, 1e-5);
  ck_assert_float_eq_tol(dq.w, dqv.w, 1e-5);
}
END_TEST

START_TEST(test_dq_to_dpi_p) {
  ksl_quaternion_t q = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t dq = ksl_quaternion(1.0, 2.0, 3.0, 4.0);
  ksl_vec3_t dpi_p;
  ksl_dq_to_dpi_p(&q, &dq, &dpi_p);
  ksl_vec3_t dpi_pv;
  dpi_pv.x = 2.0 * (q.w * dq.x - q.z * dq.y + q.y * dq.z - q.x * dq.w);
  dpi_pv.x = 2.0 * (q.z * dq.x + q.w * dq.y - q.x * dq.z - q.y * dq.w);
  dpi_pv.x = 2.0 * (-q.y * dq.x + q.x * dq.y + q.w * dq.z - q.z * dq.w);
  ck_assert_double_eq_tol(dpi_p.x, dpi_pv.x, 1e-9);
  ck_assert_double_eq_tol(dpi_p.y, dpi_pv.y, 1e-9);
  ck_assert_double_eq_tol(dpi_p.z, dpi_pv.z, 1e-9);
}
END_TEST

START_TEST(test_dqf_to_dpi_pf) {
  ksl_quaternionf_t q = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t dq = ksl_quaternionf(1.0, 2.0, 3.0, 4.0);
  ksl_vec3f_t dpi_p;
  ksl_dqf_to_dpi_pf(&q, &dq, &dpi_p);
  ksl_vec3f_t dpi_pv;
  dpi_pv.x = 2.0 * (q.w * dq.x - q.z * dq.y + q.y * dq.z - q.x * dq.w);
  dpi_pv.x = 2.0 * (q.z * dq.x + q.w * dq.y - q.x * dq.z - q.y * dq.w);
  dpi_pv.x = 2.0 * (-q.y * dq.x + q.x * dq.y + q.w * dq.z - q.z * dq.w);
  ck_assert_float_eq_tol(dpi_p.x, dpi_pv.x, 1e-5);
  ck_assert_float_eq_tol(dpi_p.y, dpi_pv.y, 1e-5);
  ck_assert_float_eq_tol(dpi_p.z, dpi_pv.z, 1e-5);
}
END_TEST

START_TEST(test_dq_to_dpi_c) {
  ksl_quaternion_t q = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t dq = ksl_quaternion(1.0, 2.0, 3.0, 4.0);
  ksl_vec3_t dpi_c;
  ksl_dq_to_dpi_c(&q, &dq, &dpi_c);
  ksl_vec3_t dpi_cv;
  dpi_cv.x = 2.0 * (q.w * dq.x + q.z * dq.y - q.y * dq.z - q.x * dq.w);
  dpi_cv.x = 2.0 * (-q.z * dq.x + q.w * dq.y + q.x * dq.z - q.y * dq.w);
  dpi_cv.x = 2.0 * (q.y * dq.x - q.x * dq.y + q.w * dq.z - q.z * dq.w);
  ck_assert_double_eq_tol(dpi_c.x, dpi_cv.x, 1e-9);
  ck_assert_double_eq_tol(dpi_c.y, dpi_cv.y, 1e-9);
  ck_assert_double_eq_tol(dpi_c.z, dpi_cv.z, 1e-9);
}
END_TEST

START_TEST(test_dqf_to_dpi_cf) {
  ksl_quaternionf_t q = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t dq = ksl_quaternionf(1.0, 2.0, 3.0, 4.0);
  ksl_vec3f_t dpi_c;
  ksl_dqf_to_dpi_cf(&q, &dq, &dpi_c);
  ksl_vec3f_t dpi_cv;
  dpi_cv.x = 2.0 * (q.w * dq.x + q.z * dq.y - q.y * dq.z - q.x * dq.w);
  dpi_cv.x = 2.0 * (-q.z * dq.x + q.w * dq.y + q.x * dq.z - q.y * dq.w);
  dpi_cv.x = 2.0 * (q.y * dq.x - q.x * dq.y + q.w * dq.z - q.z * dq.w);
  ck_assert_float_eq_tol(dpi_c.x, dpi_cv.x, 1e-5);
  ck_assert_float_eq_tol(dpi_c.y, dpi_cv.y, 1e-5);
  ck_assert_float_eq_tol(dpi_c.z, dpi_cv.z, 1e-5);
}
END_TEST

START_TEST(test_product_qq) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2 = ksl_quaternion_unit(0.1, 0.3, 0.5, 0.7);
  ksl_quaternion_t q3;
  ksl_product_qq(&q1, &q2, &q3);
  ksl_quaternion_t q3v;
  q3v.x = q1.w * q2.x - q1.z * q2.y + q1.y * q2.z + q1.x * q2.w;
  q3v.y = q1.z * q2.x + q1.w * q2.y - q1.x * q2.z + q1.y * q2.w;
  q3v.z = -q1.y * q2.x + q1.x * q2.y + q1.w * q2.z + q1.z * q2.w;
  q3v.w = -q1.x * q2.x - q1.y * q2.y - q1.z * q2.z + q1.w * q2.w;
  ck_assert_double_eq_tol(q3.x, q3v.x, 1e-9);
  ck_assert_double_eq_tol(q3.y, q3v.y, 1e-9);
  ck_assert_double_eq_tol(q3.z, q3v.z, 1e-9);
  ck_assert_double_eq_tol(q3.w, q3v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qqf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2 = ksl_quaternionf_unit(0.1, 0.3, 0.5, 0.7);
  ksl_quaternionf_t q3;
  ksl_product_qqf(&q1, &q2, &q3);
  ksl_quaternionf_t q3v;
  q3v.x = q1.w * q2.x - q1.z * q2.y + q1.y * q2.z + q1.x * q2.w;
  q3v.y = q1.z * q2.x + q1.w * q2.y - q1.x * q2.z + q1.y * q2.w;
  q3v.z = -q1.y * q2.x + q1.x * q2.y + q1.w * q2.z + q1.z * q2.w;
  q3v.w = -q1.x * q2.x - q1.y * q2.y - q1.z * q2.z + q1.w * q2.w;
  ck_assert_float_eq_tol(q3.x, q3v.x, 1e-5);
  ck_assert_float_eq_tol(q3.y, q3v.y, 1e-5);
  ck_assert_float_eq_tol(q3.z, q3v.z, 1e-5);
  ck_assert_float_eq_tol(q3.w, q3v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjq) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2 = ksl_quaternion_unit(0.1, 0.3, 0.5, 0.7);
  ksl_quaternion_t q3;
  ksl_product_qconjq(&q1, &q2, &q3);
  ksl_quaternion_t q3v;
  q3v.x = q1.w * q2.x + q1.z * q2.y - q1.y * q2.z - q1.x * q2.w;
  q3v.y = -q1.z * q2.x + q1.w * q2.y + q1.x * q2.z - q1.y * q2.w;
  q3v.z = q1.y * q2.x - q1.x * q2.y + q1.w * q2.z - q1.z * q2.w;
  q3v.w = q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
  ck_assert_double_eq_tol(q3.x, q3v.x, 1e-9);
  ck_assert_double_eq_tol(q3.y, q3v.y, 1e-9);
  ck_assert_double_eq_tol(q3.z, q3v.z, 1e-9);
  ck_assert_double_eq_tol(q3.w, q3v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjqf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2 = ksl_quaternionf_unit(0.1, 0.3, 0.5, 0.7);
  ksl_quaternionf_t q3;
  ksl_product_qconjqf(&q1, &q2, &q3);
  ksl_quaternionf_t q3v;
  q3v.x = q1.w * q2.x + q1.z * q2.y - q1.y * q2.z - q1.x * q2.w;
  q3v.y = -q1.z * q2.x + q1.w * q2.y + q1.x * q2.z - q1.y * q2.w;
  q3v.z = q1.y * q2.x - q1.x * q2.y + q1.w * q2.z - q1.z * q2.w;
  q3v.w = q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
  ck_assert_float_eq_tol(q3.x, q3v.x, 1e-5);
  ck_assert_float_eq_tol(q3.y, q3v.y, 1e-5);
  ck_assert_float_eq_tol(q3.z, q3v.z, 1e-5);
  ck_assert_float_eq_tol(q3.w, q3v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qqconj) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2 = ksl_quaternion_unit(0.1, 0.3, 0.5, 0.7);
  ksl_quaternion_t q3;
  ksl_product_qqconj(&q1, &q2, &q3);
  ksl_quaternion_t q3v;
  q3v.x = -q1.w * q2.x + q1.z * q2.y - q1.y * q2.z + q1.x * q2.w;
  q3v.y = -q1.z * q2.x - q1.w * q2.y + q1.x * q2.z + q1.y * q2.w;
  q3v.z = q1.y * q2.x - q1.x * q2.y - q1.w * q2.z + q1.z * q2.w;
  q3v.w = q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
  ck_assert_double_eq_tol(q3.x, q3v.x, 1e-9);
  ck_assert_double_eq_tol(q3.y, q3v.y, 1e-9);
  ck_assert_double_eq_tol(q3.z, q3v.z, 1e-9);
  ck_assert_double_eq_tol(q3.w, q3v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qqconjf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2 = ksl_quaternionf_unit(0.1, 0.3, 0.5, 0.7);
  ksl_quaternionf_t q3;
  ksl_product_qqconjf(&q1, &q2, &q3);
  ksl_quaternionf_t q3v;
  q3v.x = -q1.w * q2.x + q1.z * q2.y - q1.y * q2.z + q1.x * q2.w;
  q3v.y = -q1.z * q2.x - q1.w * q2.y + q1.x * q2.z + q1.y * q2.w;
  q3v.z = q1.y * q2.x - q1.x * q2.y - q1.w * q2.z + q1.z * q2.w;
  q3v.w = q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
  ck_assert_float_eq_tol(q3.x, q3v.x, 1e-5);
  ck_assert_float_eq_tol(q3.y, q3v.y, 1e-5);
  ck_assert_float_eq_tol(q3.z, q3v.z, 1e-5);
  ck_assert_float_eq_tol(q3.w, q3v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjqconj) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2 = ksl_quaternion_unit(0.1, 0.3, 0.5, 0.7);
  ksl_quaternion_t q3;
  ksl_product_qconjqconj(&q1, &q2, &q3);
  ksl_quaternion_t q3v;
  q3v.x = -q1.w * q2.x - q1.z * q2.y + q1.y * q2.z - q1.x * q2.w;
  q3v.y = q1.z * q2.x - q1.w * q2.y - q1.x * q2.z - q1.y * q2.w;
  q3v.z = -q1.y * q2.x + q1.x * q2.y - q1.w * q2.z - q1.z * q2.w;
  q3v.w = -q1.x * q2.x - q1.y * q2.y - q1.z * q2.z + q1.w * q2.w;
  ck_assert_double_eq_tol(q3.x, q3v.x, 1e-9);
  ck_assert_double_eq_tol(q3.y, q3v.y, 1e-9);
  ck_assert_double_eq_tol(q3.z, q3v.z, 1e-9);
  ck_assert_double_eq_tol(q3.w, q3v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjqconjf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2 = ksl_quaternionf_unit(0.1, 0.3, 0.5, 0.7);
  ksl_quaternionf_t q3;
  ksl_product_qconjqconjf(&q1, &q2, &q3);
  ksl_quaternionf_t q3v;
  q3v.x = -q1.w * q2.x - q1.z * q2.y + q1.y * q2.z - q1.x * q2.w;
  q3v.y = q1.z * q2.x - q1.w * q2.y - q1.x * q2.z - q1.y * q2.w;
  q3v.z = -q1.y * q2.x + q1.x * q2.y - q1.w * q2.z - q1.z * q2.w;
  q3v.w = -q1.x * q2.x - q1.y * q2.y - q1.z * q2.z + q1.w * q2.w;
  ck_assert_float_eq_tol(q3.x, q3v.x, 1e-5);
  ck_assert_float_eq_tol(q3.y, q3v.y, 1e-5);
  ck_assert_float_eq_tol(q3.z, q3v.z, 1e-5);
  ck_assert_float_eq_tol(q3.w, q3v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qxq) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qxq(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.w;
  q2v.y = dc[0] * q1.y - dc[1] * q1.z;
  q2v.z = dc[1] * q1.y + dc[0] * q1.z;
  q2v.w = -dc[1] * q1.x + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qxqf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qxqf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.w;
  q2v.y = dc[0] * q1.y - dc[1] * q1.z;
  q2v.z = dc[1] * q1.y + dc[0] * q1.z;
  q2v.w = -dc[1] * q1.x + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qqx) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qqx(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.w;
  q2v.y = dc[0] * q1.y + dc[1] * q1.z;
  q2v.z = dc[0] * q1.z - dc[1] * q1.y;
  q2v.w = dc[0] * q1.w - dc[1] * q1.x;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qqxf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qqxf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.w;
  q2v.y = dc[0] * q1.y + dc[1] * q1.z;
  q2v.z = dc[0] * q1.z - dc[1] * q1.y;
  q2v.w = dc[0] * q1.w - dc[1] * q1.x;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qxconjq) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qxconjq(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.w;
  q2v.y = dc[0] * q1.y + dc[1] * q1.z;
  q2v.z = -dc[1] * q1.y + dc[0] * q1.z;
  q2v.w = dc[1] * q1.x + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qxconjqf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qxconjqf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.w;
  q2v.y = dc[0] * q1.y + dc[1] * q1.z;
  q2v.z = -dc[1] * q1.y + dc[0] * q1.z;
  q2v.w = +dc[1] * q1.x + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjqx) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qconjqx(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.w;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.z;
  q2v.z = -dc[0] * q1.z + dc[1] * q1.y;
  q2v.w = dc[0] * q1.w + dc[1] * q1.x;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjqxf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qconjqxf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.w;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.z;
  q2v.z = -dc[0] * q1.z + dc[1] * q1.y;
  q2v.w = dc[0] * q1.w + dc[1] * q1.x;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qxqconj) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qxqconj(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.w;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.z;
  q2v.z = -dc[1] * q1.y - dc[0] * q1.z;
  q2v.w = dc[1] * q1.x + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qxqconjf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qxqconjf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.w;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.z;
  q2v.z = -dc[1] * q1.y - dc[0] * q1.z;
  q2v.w = dc[1] * q1.x + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qqxconj) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qqxconj(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.w;
  q2v.y = dc[0] * q1.y - dc[1] * q1.z;
  q2v.z = dc[0] * q1.z + dc[1] * q1.y;
  q2v.w = dc[0] * q1.w + dc[1] * q1.x;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qqxconjf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qqxconjf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.w;
  q2v.y = dc[0] * q1.y - dc[1] * q1.z;
  q2v.z = dc[0] * q1.z + dc[1] * q1.y;
  q2v.w = dc[0] * q1.w + dc[1] * q1.x;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qxconjqconj) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qxconjqconj(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.w;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.z;
  q2v.z = dc[1] * q1.y - dc[0] * q1.z;
  q2v.w = -dc[1] * q1.x + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qxconjqconjf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qxconjqconjf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.w;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.z;
  q2v.z = dc[1] * q1.y - dc[0] * q1.z;
  q2v.w = -dc[1] * q1.x + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjqxconj) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qconjqxconj(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.w;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.z;
  q2v.z = -dc[0] * q1.z - dc[1] * q1.y;
  q2v.w = dc[0] * q1.w - dc[1] * q1.x;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjqxconjf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qconjqxconjf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.w;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.z;
  q2v.z = -dc[0] * q1.z - dc[1] * q1.y;
  q2v.w = dc[0] * q1.w - dc[1] * q1.x;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST
//~~~~~

START_TEST(test_product_qyq) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qyq(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.z;
  q2v.y = dc[0] * q1.y + dc[1] * q1.w;
  q2v.z = -dc[1] * q1.x + dc[0] * q1.z;
  q2v.w = -dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qyqf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qyqf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.z;
  q2v.y = dc[0] * q1.y + dc[1] * q1.w;
  q2v.z = -dc[1] * q1.x + dc[0] * q1.z;
  q2v.w = -dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qqy) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qqy(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.z;
  q2v.y = dc[0] * q1.y + dc[1] * q1.w;
  q2v.z = dc[0] * q1.z + dc[1] * q1.x;
  q2v.w = dc[0] * q1.w - dc[1] * q1.y;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qqyf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qqyf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.z;
  q2v.y = dc[0] * q1.y + dc[1] * q1.w;
  q2v.z = dc[0] * q1.z + dc[1] * q1.x;
  q2v.w = dc[0] * q1.w - dc[1] * q1.y;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qyconjq) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qyconjq(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.z;
  q2v.y = dc[0] * q1.y - dc[1] * q1.w;
  q2v.z = dc[1] * q1.x + dc[0] * q1.z;
  q2v.w = dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qyconjqf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qyconjqf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.z;
  q2v.y = dc[0] * q1.y - dc[1] * q1.w;
  q2v.z = dc[1] * q1.x + dc[0] * q1.z;
  q2v.w = dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjqy) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qconjqy(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.z;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.w;
  q2v.z = -dc[0] * q1.z - dc[1] * q1.x;
  q2v.w = dc[0] * q1.w + dc[1] * q1.y;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjqyf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qconjqyf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.z;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.w;
  q2v.z = -dc[0] * q1.z - dc[1] * q1.x;
  q2v.w = dc[0] * q1.w + dc[1] * q1.y;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST
//
START_TEST(test_product_qyqconj) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qyqconj(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.z;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.w;
  q2v.z = dc[1] * q1.x - dc[0] * q1.z;
  q2v.w = dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qyqconjf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qyqconjf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.z;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.w;
  q2v.z = dc[1] * q1.x - dc[0] * q1.z;
  q2v.w = dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qqyconj) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qqyconj(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.z;
  q2v.y = dc[0] * q1.y - dc[1] * q1.w;
  q2v.z = dc[0] * q1.z - dc[1] * q1.x;
  q2v.w = dc[0] * q1.w + dc[1] * q1.y;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qqyconjf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qqyconjf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.z;
  q2v.y = dc[0] * q1.y - dc[1] * q1.w;
  q2v.z = dc[0] * q1.z - dc[1] * q1.x;
  q2v.w = dc[0] * q1.w + dc[1] * q1.y;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qyconjqconj) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qyconjqconj(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.z;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.w;
  q2v.z = -dc[1] * q1.x - dc[0] * q1.z;
  q2v.w = -dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qyconjqconjf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qyconjqconjf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.z;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.w;
  q2v.z = -dc[1] * q1.x - dc[0] * q1.z;
  q2v.w = -dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjqyconj) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qconjqyconj(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.z;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.w;
  q2v.z = +dc[1] * q1.x - dc[0] * q1.z;
  q2v.w = -dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjqyconjf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qconjqyconjf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.z;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.w;
  q2v.z = +dc[1] * q1.x - dc[0] * q1.z;
  q2v.w = -dc[1] * q1.y + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST
//~~~~~~~~~~~~~~~~~~~

START_TEST(test_product_qzq) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qzq(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.y;
  q2v.y = dc[0] * q1.y + dc[1] * q1.x;
  q2v.z = dc[1] * q1.w + dc[0] * q1.z;
  q2v.w = -dc[1] * q1.z + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qzqf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qzqf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.y;
  q2v.y = dc[0] * q1.y + dc[1] * q1.x;
  q2v.z = dc[1] * q1.w + dc[0] * q1.z;
  q2v.w = -dc[1] * q1.z + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qqz) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qqz(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.y;
  q2v.y = dc[0] * q1.y - dc[1] * q1.x;
  q2v.z = dc[0] * q1.z + dc[1] * q1.w;
  q2v.w = dc[0] * q1.w - dc[1] * q1.z;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qqzf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qqzf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.y;
  q2v.y = dc[0] * q1.y - dc[1] * q1.x;
  q2v.z = dc[0] * q1.z + dc[1] * q1.w;
  q2v.w = dc[0] * q1.w - dc[1] * q1.z;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qzconjq) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qzconjq(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.y;
  q2v.y = dc[0] * q1.y - dc[1] * q1.x;
  q2v.z = -dc[1] * q1.w + dc[0] * q1.z;
  q2v.w = dc[1] * q1.z + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qzconjqf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qzconjqf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x + dc[1] * q1.y;
  q2v.y = dc[0] * q1.y - dc[1] * q1.x;
  q2v.z = -dc[1] * q1.w + dc[0] * q1.z;
  q2v.w = dc[1] * q1.z + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjqz) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qconjqz(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.y;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.x;
  q2v.z = -dc[0] * q1.z + dc[1] * q1.w;
  q2v.w = dc[0] * q1.w + dc[1] * q1.z;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjqzf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qconjqzf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.y;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.x;
  q2v.z = -dc[0] * q1.z + dc[1] * q1.w;
  q2v.w = dc[0] * q1.w + dc[1] * q1.z;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qzqconj) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qzqconj(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.y;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.x;
  q2v.z = dc[1] * q1.w - dc[0] * q1.z;
  q2v.w = dc[1] * q1.z + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qzqconjf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qzqconjf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.y;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.x;
  q2v.z = dc[1] * q1.w - dc[0] * q1.z;
  q2v.w = dc[1] * q1.z + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qqzconj) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qqzconj(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.y;
  q2v.y = dc[0] * q1.y + dc[1] * q1.x;
  q2v.z = dc[0] * q1.z - dc[1] * q1.w;
  q2v.w = dc[0] * q1.w + dc[1] * q1.z;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qqzconjf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qqzconjf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = dc[0] * q1.x - dc[1] * q1.y;
  q2v.y = dc[0] * q1.y + dc[1] * q1.x;
  q2v.z = dc[0] * q1.z - dc[1] * q1.w;
  q2v.w = dc[0] * q1.w + dc[1] * q1.z;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qzconjqconj) {
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_qzconjqconj(dc, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.y;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.x;
  q2v.z = -dc[1] * q1.w - dc[0] * q1.z;
  q2v.w = -dc[1] * q1.z + dc[0] * q1.w;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qzconjqconjf) {
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_qzconjqconjf(dc, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x - dc[1] * q1.y;
  q2v.y = -dc[0] * q1.y + dc[1] * q1.x;
  q2v.z = -dc[1] * q1.w - dc[0] * q1.z;
  q2v.w = -dc[1] * q1.z + dc[0] * q1.w;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjqzconj) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  double dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternion_t q2;
  ksl_product_qconjqzconj(&q1, dc, &q2);
  ksl_quaternion_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.y;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.x;
  q2v.z = -dc[0] * q1.z - dc[1] * q1.w;
  q2v.w = dc[0] * q1.w - dc[1] * q1.z;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjqzconjf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  float dc[2] = {cos(1.0), sin(1.0)};
  ksl_quaternionf_t q2;
  ksl_product_qconjqzconjf(&q1, dc, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = -dc[0] * q1.x + dc[1] * q1.y;
  q2v.y = -dc[0] * q1.y - dc[1] * q1.x;
  q2v.z = -dc[0] * q1.z - dc[1] * q1.w;
  q2v.w = dc[0] * q1.w - dc[1] * q1.z;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qv) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_vec3_t v = ksl_vec3(1.0, 2.0, 3.0);
  ksl_quaternion_t q2;
  ksl_product_qv(&q1, &v, &q2);
  ksl_quaternion_t q2v;
  q2v.x = q1.w * v.x - q1.z * v.y + q1.y * v.z;
  q2v.y = q1.z * v.x + q1.w * v.y - q1.x * v.z;
  q2v.z = -q1.y * v.x + q1.x * v.y + q1.w * v.z;
  q2v.w = -q1.x * v.x - q1.y * v.y - q1.z * v.z;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qvf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_vec3f_t v = ksl_vec3f(1.0, 2.0, 3.0);
  ksl_quaternionf_t q2;
  ksl_product_qvf(&q1, &v, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = q1.w * v.x - q1.z * v.y + q1.y * v.z;
  q2v.y = q1.z * v.x + q1.w * v.y - q1.x * v.z;
  q2v.z = -q1.y * v.x + q1.x * v.y + q1.w * v.z;
  q2v.w = -q1.x * v.x - q1.y * v.y - q1.z * v.z;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_qconjv) {
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_vec3_t v = ksl_vec3(1.0, 2.0, 3.0);
  ksl_quaternion_t q2;
  ksl_product_qconjv(&q1, &v, &q2);
  ksl_quaternion_t q2v;
  q2v.x = q1.w * v.x + q1.z * v.y - q1.y * v.z;
  q2v.y = -q1.z * v.x + q1.w * v.y + q1.x * v.z;
  q2v.z = q1.y * v.x - q1.x * v.y + q1.w * v.z;
  q2v.w = q1.x * v.x + q1.y * v.y + q1.z * v.z;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_qconjvf) {
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_vec3f_t v = ksl_vec3f(1.0, 2.0, 3.0);
  ksl_quaternionf_t q2;
  ksl_product_qconjvf(&q1, &v, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = q1.w * v.x + q1.z * v.y - q1.y * v.z;
  q2v.y = -q1.z * v.x + q1.w * v.y + q1.x * v.z;
  q2v.z = q1.y * v.x - q1.x * v.y + q1.w * v.z;
  q2v.w = q1.x * v.x + q1.y * v.y + q1.z * v.z;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_vq) {
  ksl_vec3_t v = ksl_vec3(1.0, 2.0, 3.0);
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_vq(&v, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = q1.w * v.x + q1.z * v.y - q1.y * v.z;
  q2v.y = -q1.z * v.x + q1.w * v.y + q1.x * v.z;
  q2v.z = q1.y * v.x - q1.x * v.y + q1.w * v.z;
  q2v.w = -q1.x * v.x - q1.y * v.y - q1.z * v.z;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_vqf) {
  ksl_vec3f_t v = ksl_vec3f(1.0, 2.0, 3.0);
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_vqf(&v, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = q1.w * v.x + q1.z * v.y - q1.y * v.z;
  q2v.y = -q1.z * v.x + q1.w * v.y + q1.x * v.z;
  q2v.z = q1.y * v.x - q1.x * v.y + q1.w * v.z;
  q2v.w = -q1.x * v.x - q1.y * v.y - q1.z * v.z;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

START_TEST(test_product_vqconj) {
  ksl_vec3_t v = ksl_vec3(1.0, 2.0, 3.0);
  ksl_quaternion_t q1 = ksl_quaternion_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternion_t q2;
  ksl_product_vqconj(&v, &q1, &q2);
  ksl_quaternion_t q2v;
  q2v.x = q1.w * v.x - q1.z * v.y + q1.y * v.z;
  q2v.y = q1.z * v.x + q1.w * v.y - q1.x * v.z;
  q2v.z = -q1.y * v.x + q1.x * v.y + q1.w * v.z;
  q2v.w = q1.x * v.x + q1.y * v.y + q1.z * v.z;
  ck_assert_double_eq_tol(q2.x, q2v.x, 1e-9);
  ck_assert_double_eq_tol(q2.y, q2v.y, 1e-9);
  ck_assert_double_eq_tol(q2.z, q2v.z, 1e-9);
  ck_assert_double_eq_tol(q2.w, q2v.w, 1e-9);
}
END_TEST

START_TEST(test_product_vqconjf) {
  ksl_vec3f_t v = ksl_vec3f(1.0, 2.0, 3.0);
  ksl_quaternionf_t q1 = ksl_quaternionf_unit(0.2, 0.4, 0.6, 0.8);
  ksl_quaternionf_t q2;
  ksl_product_vqconjf(&v, &q1, &q2);
  ksl_quaternionf_t q2v;
  q2v.x = q1.w * v.x - q1.z * v.y + q1.y * v.z;
  q2v.y = q1.z * v.x + q1.w * v.y - q1.x * v.z;
  q2v.z = -q1.y * v.x + q1.x * v.y + q1.w * v.z;
  q2v.w = q1.x * v.x + q1.y * v.y + q1.z * v.z;
  ck_assert_float_eq_tol(q2.x, q2v.x, 1e-5);
  ck_assert_float_eq_tol(q2.y, q2v.y, 1e-5);
  ck_assert_float_eq_tol(q2.z, q2v.z, 1e-5);
  ck_assert_float_eq_tol(q2.w, q2v.w, 1e-5);
}
END_TEST

Suite* quaternion_suite(void) {
  Suite* s = suite_create("quaternion");

  /* Core test case */
  TCase* tc_core = tcase_create("core");

  tcase_add_test(tc_core, test_quaternion_create);
  tcase_add_test(tc_core, test_quaternionf_create);
  tcase_add_test(tc_core, test_quaternion_unit);
  tcase_add_test(tc_core, test_quaternionf_unit);
  tcase_add_test(tc_core, test_quaternion_alloc);
  tcase_add_test(tc_core, test_quaternionf_alloc);
  tcase_add_test(tc_core, test_quaternion_maxIndex);
  tcase_add_test(tc_core, test_quaternionf_maxIndex);
  tcase_add_test(tc_core, test_dot_qq);
  tcase_add_test(tc_core, test_dot_qqf);
  tcase_add_test(tc_core, test_quaternion_normalize);
  tcase_add_test(tc_core, test_quaternionf_normalize);
  tcase_add_test(tc_core, test_axpy_qq);
  tcase_add_test(tc_core, test_axpy_qqf);
  tcase_add_test(tc_core, test_xpy_qq);
  tcase_add_test(tc_core, test_xpy_qqf);
  tcase_add_test(tc_core, test_nxpy_qq);
  tcase_add_test(tc_core, test_nxpy_qqf);
  tcase_add_test(tc_core, test_mat3x3_to_quaternion);
  tcase_add_test(tc_core, test_mat3x3f_to_quaternionf);
  tcase_add_test(tc_core, test_quaternion_to_mat3x3);
  tcase_add_test(tc_core, test_quaternionf_to_mat3x3f);
  tcase_add_test(tc_core, test_dpi_p_to_dq);
  tcase_add_test(tc_core, test_dpi_pf_to_dqf);
  tcase_add_test(tc_core, test_dpi_c_to_dq);
  tcase_add_test(tc_core, test_dpi_cf_to_dqf);
  tcase_add_test(tc_core, test_dq_to_dpi_p);
  tcase_add_test(tc_core, test_dqf_to_dpi_pf);
  tcase_add_test(tc_core, test_dq_to_dpi_c);
  tcase_add_test(tc_core, test_dqf_to_dpi_cf);
  tcase_add_test(tc_core, test_product_qq);
  tcase_add_test(tc_core, test_product_qqf);
  tcase_add_test(tc_core, test_product_qconjq);
  tcase_add_test(tc_core, test_product_qconjqf);
  tcase_add_test(tc_core, test_product_qqconj);
  tcase_add_test(tc_core, test_product_qqconjf);
  tcase_add_test(tc_core, test_product_qconjqconj);
  tcase_add_test(tc_core, test_product_qconjqconjf);
  tcase_add_test(tc_core, test_product_qxq);
  tcase_add_test(tc_core, test_product_qxqf);
  tcase_add_test(tc_core, test_product_qqx);
  tcase_add_test(tc_core, test_product_qqxf);
  tcase_add_test(tc_core, test_product_qxconjq);
  tcase_add_test(tc_core, test_product_qxconjqf);
  tcase_add_test(tc_core, test_product_qconjqx);
  tcase_add_test(tc_core, test_product_qconjqxf);
  tcase_add_test(tc_core, test_product_qxqconj);
  tcase_add_test(tc_core, test_product_qxqconjf);
  tcase_add_test(tc_core, test_product_qqxconj);
  tcase_add_test(tc_core, test_product_qqxconjf);
  tcase_add_test(tc_core, test_product_qxconjqconj);
  tcase_add_test(tc_core, test_product_qxconjqconjf);
  tcase_add_test(tc_core, test_product_qconjqxconj);
  tcase_add_test(tc_core, test_product_qconjqxconjf);
  tcase_add_test(tc_core, test_product_qyq);
  tcase_add_test(tc_core, test_product_qyqf);
  tcase_add_test(tc_core, test_product_qqy);
  tcase_add_test(tc_core, test_product_qqyf);
  tcase_add_test(tc_core, test_product_qyconjq);
  tcase_add_test(tc_core, test_product_qyconjqf);
  tcase_add_test(tc_core, test_product_qconjqy);
  tcase_add_test(tc_core, test_product_qconjqyf);
  tcase_add_test(tc_core, test_product_qyqconj);
  tcase_add_test(tc_core, test_product_qyqconjf);
  tcase_add_test(tc_core, test_product_qqyconj);
  tcase_add_test(tc_core, test_product_qqyconjf);
  tcase_add_test(tc_core, test_product_qyconjqconj);
  tcase_add_test(tc_core, test_product_qyconjqconjf);
  tcase_add_test(tc_core, test_product_qconjqyconj);
  tcase_add_test(tc_core, test_product_qconjqyconjf);
  tcase_add_test(tc_core, test_product_qyq);
  tcase_add_test(tc_core, test_product_qyqf);
  tcase_add_test(tc_core, test_product_qqy);
  tcase_add_test(tc_core, test_product_qqyf);
  tcase_add_test(tc_core, test_product_qyconjq);
  tcase_add_test(tc_core, test_product_qyconjqf);
  tcase_add_test(tc_core, test_product_qconjqy);
  tcase_add_test(tc_core, test_product_qconjqyf);
  tcase_add_test(tc_core, test_product_qyqconj);
  tcase_add_test(tc_core, test_product_qyqconjf);
  tcase_add_test(tc_core, test_product_qqyconj);
  tcase_add_test(tc_core, test_product_qqyconjf);
  tcase_add_test(tc_core, test_product_qyconjqconj);
  tcase_add_test(tc_core, test_product_qyconjqconjf);
  tcase_add_test(tc_core, test_product_qconjqyconj);
  tcase_add_test(tc_core, test_product_qconjqyconjf);
  tcase_add_test(tc_core, test_product_qzq);
  tcase_add_test(tc_core, test_product_qzqf);
  tcase_add_test(tc_core, test_product_qqz);
  tcase_add_test(tc_core, test_product_qqzf);
  tcase_add_test(tc_core, test_product_qzconjq);
  tcase_add_test(tc_core, test_product_qzconjqf);
  tcase_add_test(tc_core, test_product_qconjqz);
  tcase_add_test(tc_core, test_product_qconjqzf);
  tcase_add_test(tc_core, test_product_qzqconj);
  tcase_add_test(tc_core, test_product_qzqconjf);
  tcase_add_test(tc_core, test_product_qqzconj);
  tcase_add_test(tc_core, test_product_qqzconjf);
  tcase_add_test(tc_core, test_product_qzconjqconj);
  tcase_add_test(tc_core, test_product_qzconjqconjf);
  tcase_add_test(tc_core, test_product_qconjqzconj);
  tcase_add_test(tc_core, test_product_qconjqzconjf);
  tcase_add_test(tc_core, test_product_qv);
  tcase_add_test(tc_core, test_product_qvf);
  tcase_add_test(tc_core, test_product_qconjv);
  tcase_add_test(tc_core, test_product_qconjvf);
  tcase_add_test(tc_core, test_product_vq);
  tcase_add_test(tc_core, test_product_vqf);
  tcase_add_test(tc_core, test_product_vqconj);
  tcase_add_test(tc_core, test_product_vqconjf);

  suite_add_tcase(s, tc_core);

  return s;
}

int main(void) {
  int number_failed;
  Suite* s = quaternion_suite();
  SRunner* sr = srunner_create(s);
  srunner_set_log(sr, "check_quaternion.log");
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
