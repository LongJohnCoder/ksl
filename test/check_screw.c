#include <config.h>

#include <check.h>
#include <stdio.h>
#include <stdlib.h>

#include "coscrew.h"
#include "screw.h"

static FILE* dbg;

START_TEST(test_screw_create) {
  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  int i;
  double v;
  for(i = 0, v = 1.0; i < 6; i++, v += 1.0) {
    ck_assert(s.at[i] == v);
  }
  ck_assert((s.lin.x == 1.0) && (s.m0 == 1.0));
  ck_assert((s.lin.y == 2.0) && (s.m1 == 2.0));
  ck_assert((s.lin.z == 3.0) && (s.m2 == 3.0));
  ck_assert((s.ang.x == 4.0) && (s.m3 == 4.0));
  ck_assert((s.ang.y == 5.0) && (s.m4 == 5.0));
  ck_assert((s.ang.z == 6.0) && (s.m5 == 6.0));
}
END_TEST

START_TEST(test_screwf_create) {
  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  float v;
  int i;
  for(i = 0, v = 1.0; i < 6; i++, v += 1.0) {
    ck_assert(s.at[i] == v);
  }
  ck_assert((s.lin.x == 1.0) && (s.m0 == 1.0));
  ck_assert((s.lin.y == 2.0) && (s.m1 == 2.0));
  ck_assert((s.lin.z == 3.0) && (s.m2 == 3.0));
  ck_assert((s.ang.x == 4.0) && (s.m3 == 4.0));
  ck_assert((s.ang.y == 5.0) && (s.m4 == 5.0));
  ck_assert((s.ang.z == 6.0) && (s.m5 == 6.0));
}
END_TEST

START_TEST(test_screw_alloc) {
  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t* s3 = ksl_screw_alloc(3);
  for(int i = 0; i < 3; i++) {
    ksl_screw_copy(&s, &s3[i]);
    int j;
    double v;
    for(j = 0, v = 1.0; j < 6; j++, v += 1.0) {
      ck_assert(s3[i].at[j] == v);
    }
  }
}
END_TEST

START_TEST(test_screwf_alloc) {
  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t* s3 = ksl_screwf_alloc(3);
  for(int i = 0; i < 3; i++) {
    ksl_screwf_copy(&s, &s3[i]);
    int j;
    float v;
    for(j = 0, v = 1.0; j < 6; j++, v += 1.0) {
      ck_assert(s3[i].at[j] == v);
    }
  }
}
END_TEST

START_TEST(test_screw_norm) {
  ksl_screw_t s1 = {{7.0, 8.0, -3.0, 2.0, 1.0, -2.0}};
  ck_assert_double_eq_tol(ksl_screw_norm(&s1), 3.0, 1.e-9);
  ksl_screw_t s2 = {{2.0, 1.0, -2.0, 0.0, 0.0, 0.0}};
  ck_assert_double_eq_tol(ksl_screw_norm(&s2), 3.0, 1.e-9);
}
END_TEST

START_TEST(test_screwf_norm) {
  ksl_screwf_t s1 = {{7.0, 8.0, -3.0, 2.0, 1.0, -2.0}};
  ck_assert_float_eq_tol(ksl_screwf_norm(&s1), 3.0, 1.e-6);
  ksl_screwf_t s2 = {{2.0, 1.0, -2.0, 0.0, 0.0, 0.0}};
  ck_assert_float_eq_tol(ksl_screwf_norm(&s2), 3.0, 1.e-6);
}
END_TEST

START_TEST(test_screw_normalize) {
  ksl_screw_t s1 = {{7.0, 8.0, -3.0, 2.0, 1.0, -2.0}};
  ksl_screw_normalize(&s1);
  ck_assert_double_eq_tol(s1.m0, 7.0 / 3.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m1, 8.0 / 3.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m2, -3.0 / 3.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m3, 2.0 / 3.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m4, 1.0 / 3.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m5, -2.0 / 3.0, 1.e-9);
  ksl_screw_t s2 = {{2.0, 1.0, -2.0, 0.0, 0.0, 0.0}};
  ksl_screw_normalize(&s2);
  ck_assert_double_eq_tol(s2.m0, 2.0 / 3.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 1.0 / 3.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, -2.0 / 3.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 0.0, 1.e-9);
}
END_TEST

START_TEST(test_screwf_normalize) {
  ksl_screwf_t s1 = {{7.0, 8.0, -3.0, 2.0, 1.0, -2.0}};
  ksl_screwf_normalize(&s1);
  ck_assert_float_eq_tol(s1.m0, 7.0 / 3.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m1, 8.0 / 3.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m2, -3.0 / 3.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m3, 2.0 / 3.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m4, 1.0 / 3.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m5, -2.0 / 3.0, 1.e-6);
  ksl_screwf_t s2 = {{2.0, 1.0, -2.0, 0.0, 0.0, 0.0}};
  ksl_screwf_normalize(&s2);
  ck_assert_float_eq_tol(s2.m0, 2.0 / 3.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 1.0 / 3.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, -2.0 / 3.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 0.0, 1.e-6);
}
END_TEST

START_TEST(test_dot_cs) {
  ksl_screw_t s = {{1.0, 2.0, 3.0, -3.0, -2.0, -1.0}};
  ksl_coscrew_t c = {{-1.0, -2.0, -3.0, 3.0, 2.0, 1.0}};
  ck_assert_double_eq_tol(ksl_dot_cs(&c, &s), -28.0, 1.e-9);
}
END_TEST

START_TEST(test_dot_csf) {
  ksl_screwf_t s = {{1.0, 2.0, 3.0, -3.0, -2.0, -1.0}};
  ksl_coscrewf_t c = {{-1.0, -2.0, -3.0, 3.0, 2.0, 1.0}};
  ck_assert_float_eq_tol(ksl_dot_csf(&c, &s), -28.0, 1.e-6);
}
END_TEST

START_TEST(test_axpy_ss) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{-2.0, -4.0, -6.0, -8.0, -10.0, -12.0}};
  ksl_axpy_ss(5.0, &s1, &s2);
  ck_assert_double_eq_tol(s2.m0, 3.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 6.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 9.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 12.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 15.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 18.0, 1.e-9);
}
END_TEST

START_TEST(test_axpy_ssf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{-2.0, -4.0, -6.0, -8.0, -10.0, -12.0}};
  ksl_axpy_ssf(5.0, &s1, &s2);
  ck_assert_float_eq_tol(s2.m0, 3.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 6.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 9.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 12.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 15.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 18.0, 1.e-6);
}
END_TEST

START_TEST(test_screw_copy) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0}};
  ksl_screw_copy(&s1, &s2);
  ck_assert_double_eq_tol(s2.m0, 1.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 2.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 3.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 4.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 5.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 6.0, 1.e-9);
}
END_TEST

START_TEST(test_screwf_copy) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0}};
  ksl_screwf_copy(&s1, &s2);
  ck_assert_float_eq_tol(s2.m0, 1.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 2.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 3.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 4.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 5.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 6.0, 1.e-6);
}
END_TEST

START_TEST(test_screw_inverted) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_inverted(&s1, &s2);
  ck_assert_double_eq_tol(s2.m0, -1.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, -2.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, -3.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, -4.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, -5.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, -6.0, 1.e-9);
}
END_TEST

START_TEST(test_screwf_inverted) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_inverted(&s1, &s2);
  ck_assert_float_eq_tol(s2.m0, -1.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, -2.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, -3.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, -4.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, -5.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, -6.0, 1.e-6);
}
END_TEST

START_TEST(test_screw_invert) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_invert(&s1);
  ck_assert_double_eq_tol(s1.m0, -1.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m1, -2.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m2, -3.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m3, -4.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m4, -5.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m5, -6.0, 1.e-9);
}
END_TEST

START_TEST(test_screwf_invert) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_invert(&s1);
  ck_assert_float_eq_tol(s1.m0, -1.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m1, -2.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m2, -3.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m3, -4.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m4, -5.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m5, -6.0, 1.e-6);
}
END_TEST

START_TEST(test_product_as) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0}};
  ksl_product_as(2.0, &s1, &s2);
  ck_assert_double_eq_tol(s2.m0, 2.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 4.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 6.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 8.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 10.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 12.0, 1.e-9);
}
END_TEST

START_TEST(test_product_asf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0}};
  ksl_product_asf(2.0, &s1, &s2);
  ck_assert_float_eq_tol(s2.m0, 2.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 4.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 6.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 8.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 10.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 12.0, 1.e-6);
}
END_TEST

START_TEST(test_add_ss) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_add_ss(&s1, &s2, &s2);
  ck_assert_double_eq_tol(s2.m0, 11.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 22.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 33.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 44.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 55.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 66.0, 1.e-9);
}
END_TEST

START_TEST(test_add_ssf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_add_ssf(&s1, &s2, &s2);
  ck_assert_float_eq_tol(s2.m0, 11.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 22.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 33.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 44.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 55.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 66.0, 1.e-6);
}
END_TEST

START_TEST(test_xpy_ss) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_xpy_ss(&s1, &s2);
  ck_assert_double_eq_tol(s2.m0, 11.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 22.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 33.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 44.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 55.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 66.0, 1.e-9);
}
END_TEST

START_TEST(test_xpy_ssf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_xpy_ssf(&s1, &s2);
  ck_assert_float_eq_tol(s2.m0, 11.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 22.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 33.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 44.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 55.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 66.0, 1.e-6);
}
END_TEST

START_TEST(test_subtract_ss) {
  ksl_screw_t s1 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_screw_t s2 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_subtract_ss(&s1, &s2, &s2);
  ck_assert_double_eq_tol(s2.m0, 9.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 18.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 27.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 36.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 45.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 54.0, 1.e-9);
}
END_TEST

START_TEST(test_subtract_ssf) {
  ksl_screwf_t s1 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_screwf_t s2 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_subtract_ssf(&s1, &s2, &s2);
  ck_assert_float_eq_tol(s2.m0, 9.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 18.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 27.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 36.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 45.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 54.0, 1.e-6);
}
END_TEST

START_TEST(test_nxpy_ss) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_nxpy_ss(&s1, &s2);
  ck_assert_double_eq_tol(s2.m0, 9.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 18.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 27.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 36.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 45.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 54.0, 1.e-9);
}
END_TEST

START_TEST(test_nxpy_ssf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_nxpy_ssf(&s1, &s2);
  ck_assert_float_eq_tol(s2.m0, 9.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 18.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 27.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 36.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 45.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 54.0, 1.e-6);
}
END_TEST

START_TEST(test_add_sst) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{10.0, 20.0, 30.0, 0.0, 0.0, 0.0}};
  ksl_add_sst(&s1, &s2, &s2);
  ck_assert_double_eq_tol(s2.m0, 11.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 22.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 33.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 4.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 5.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 6.0, 1.e-9);
}
END_TEST

START_TEST(test_add_sstf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{10.0, 20.0, 30.0, 0.0, 0.0, 0.0}};
  ksl_add_sstf(&s1, &s2, &s2);
  ck_assert_float_eq_tol(s2.m0, 11.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 22.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 33.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 4.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 5.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 6.0, 1.e-6);
}
END_TEST

START_TEST(test_htx) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htx(&D, &s);
  ck_assert_double_eq_tol(s.m0, 1.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, 2.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, 3.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 0.0, 1.e-9);
}
END_TEST

START_TEST(test_htxf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htxf(&D, &s);
  ck_assert_float_eq_tol(s.m0, 1.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, 2.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, 3.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 0.0, 1.e-6);
}
END_TEST

START_TEST(test_hty) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hty(&D, &s);
  ck_assert_double_eq_tol(s.m0, 4.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, 5.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, 6.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 0.0, 1.e-9);
}
END_TEST

START_TEST(test_htyf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htyf(&D, &s);
  ck_assert_float_eq_tol(s.m0, 4.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, 5.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, 6.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 0.0, 1.e-6);
}
END_TEST

START_TEST(test_htz) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htz(&D, &s);
  ck_assert_double_eq_tol(s.m0, 7.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, 8.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, 9.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 0.0, 1.e-9);
}
END_TEST

START_TEST(test_htzf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htzf(&D, &s);
  ck_assert_float_eq_tol(s.m0, 7.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, 8.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, 9.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 0.0, 1.e-6);
}
END_TEST

START_TEST(test_hrx) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hrx(&D, &s);
  ck_assert_double_eq_tol(s.m0, 9.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, -18.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, 9.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 1.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 2.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 3.0, 1.e-9);
}
END_TEST

START_TEST(test_hrxf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hrxf(&D, &s);
  ck_assert_float_eq_tol(s.m0, 9.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, -18.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, 9.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 1.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 2.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 3.0, 1.e-6);
}
END_TEST

START_TEST(test_hry) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hry(&D, &s);
  ck_assert_double_eq_tol(s.m0, 6.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, -12.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, 6.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 4.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 5.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 6.0, 1.e-9);
}
END_TEST

START_TEST(test_hryf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hryf(&D, &s);
  ck_assert_float_eq_tol(s.m0, 6.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, -12.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, 6.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 4.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 5.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 6.0, 1.e-6);
}
END_TEST

START_TEST(test_hrz) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hrz(&D, &s);
  ck_assert_double_eq_tol(s.m0, 3.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, -6.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, 3.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 7.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 8.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 9.0, 1.e-9);
}
END_TEST

START_TEST(test_hrzf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hrzf(&D, &s);
  ck_assert_float_eq_tol(s.m0, 3.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, -6.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, 3.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 7.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 8.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 9.0, 1.e-6);
}
END_TEST

START_TEST(test_htxinv) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htxinv(&D, &s);
  ck_assert_double_eq_tol(s.m0, -1.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, -2.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, -3.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 0.0, 1.e-9);
}
END_TEST

START_TEST(test_htxinvf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htxinvf(&D, &s);
  ck_assert_float_eq_tol(s.m0, -1.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, -2.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, -3.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 0.0, 1.e-6);
}
END_TEST

START_TEST(test_htyinv) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htyinv(&D, &s);
  ck_assert_double_eq_tol(s.m0, -4.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, -5.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, -6.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 0.0, 1.e-9);
}
END_TEST

START_TEST(test_htyinvf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htyinvf(&D, &s);
  ck_assert_float_eq_tol(s.m0, -4.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, -5.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, -6.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 0.0, 1.e-6);
}
END_TEST

START_TEST(test_htzinv) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htzinv(&D, &s);
  ck_assert_double_eq_tol(s.m0, -7.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, -8.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, -9.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 0.0, 1.e-9);
}
END_TEST

START_TEST(test_htzinvf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htzinvf(&D, &s);
  ck_assert_float_eq_tol(s.m0, -7.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, -8.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, -9.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, 0.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, 0.0, 1.e-6);
}
END_TEST

START_TEST(test_hrxinv) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hrxinv(&D, &s);
  ck_assert_double_eq_tol(s.m0, -9.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, 18.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, -9.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, -1.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, -2.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, -3.0, 1.e-9);
}
END_TEST

START_TEST(test_hrxinvf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hrxinvf(&D, &s);
  ck_assert_float_eq_tol(s.m0, -9.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, 18.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, -9.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, -1.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, -2.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, -3.0, 1.e-6);
}
END_TEST

START_TEST(test_hryinv) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hryinv(&D, &s);
  ck_assert_double_eq_tol(s.m0, -6.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, 12.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, -6.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, -4.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, -5.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, -6.0, 1.e-9);
}
END_TEST

START_TEST(test_hryinvf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hryinvf(&D, &s);
  ck_assert_float_eq_tol(s.m0, -6.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, 12.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, -6.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, -4.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, -5.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, -6.0, 1.e-6);
}
END_TEST

START_TEST(test_hrzinv) {
  ksl_SE3_t D;
  D =
    ksl_SE3_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hrzinv(&D, &s);
  ck_assert_double_eq_tol(s.m0, -3.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, 6.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, -3.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, -7.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, -8.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, -9.0, 1.e-9);
}
END_TEST

START_TEST(test_hrzinvf) {
  ksl_SE3f_t D;
  D =
    ksl_SE3f_cmo(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  ksl_screwf_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_hrzinvf(&D, &s);
  ck_assert_float_eq_tol(s.m0, -3.0, 1.e-6);
  ck_assert_float_eq_tol(s.m1, 6.0, 1.e-6);
  ck_assert_float_eq_tol(s.m2, -3.0, 1.e-6);
  ck_assert_float_eq_tol(s.m3, -7.0, 1.e-6);
  ck_assert_float_eq_tol(s.m4, -8.0, 1.e-6);
  ck_assert_float_eq_tol(s.m5, -9.0, 1.e-6);
}
END_TEST

START_TEST(test_cross_ss) {
  ksl_screw_t s1 = {{4.0, 5.0, 6.0, 1.0, 2.0, 3.0}};
  ksl_screw_t s2 = {{6.0, 5.0, 4.0, 3.0, 2.0, 1.0}};
  fprintf(dbg, "\nIn cross_ss\n");
  fprintf(dbg, "s1: [%g, %g, %g, %g, %g, %g]\n", s1.m0, s1.m1, s1.m2, s1.m3,
          s1.m4, s1.m5);
  fprintf(dbg, "s2: [%g, %g, %g, %g, %g, %g]\n", s2.m0, s2.m1, s2.m2, s2.m3,
          s2.m4, s2.m5);
  ksl_cross_ss(&s1, &s2, &s2);
  fprintf(dbg, "s1: [%g, %g, %g, %g, %g, %g]\n", s1.m0, s1.m1, s1.m2, s1.m3,
          s1.m4, s1.m5);
  fprintf(dbg, "s2: [%g, %g, %g, %g, %g, %g]\n", s2.m0, s2.m1, s2.m2, s2.m3,
          s2.m4, s2.m5);
  ck_assert_double_eq_tol(s2.m0, -14.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 28.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, -14.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, -4.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 8.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, -4.0, 1.e-9);
}
END_TEST

START_TEST(test_cross_ssf) {
  ksl_screwf_t s1 = {{4.0, 5.0, 6.0, 1.0, 2.0, 3.0}};
  ksl_screwf_t s2 = {{6.0, 5.0, 4.0, 3.0, 2.0, 1.0}};
  fprintf(dbg, "\nIn cross_ssf\n");
  fprintf(dbg, "s1: [%g, %g, %g, %g, %g, %g]\n", s1.m0, s1.m1, s1.m2, s1.m3,
          s1.m4, s1.m5);
  fprintf(dbg, "s2: [%g, %g, %g, %g, %g, %g]\n", s2.m0, s2.m1, s2.m2, s2.m3,
          s2.m4, s2.m5);
  ksl_cross_ssf(&s1, &s2, &s2);
  fprintf(dbg, "s1: [%g, %g, %g, %g, %g, %g]\n", s1.m0, s1.m1, s1.m2, s1.m3,
          s1.m4, s1.m5);
  fprintf(dbg, "s2: [%g, %g, %g, %g, %g, %g]\n", s2.m0, s2.m1, s2.m2, s2.m3,
          s2.m4, s2.m5);
  // ck_assert_float_eq_tol(s2.m0, -14.0, 1.e-6);
  // ck_assert_float_eq_tol(s2.m1, 28.0, 1.e-6);
  // ck_assert_float_eq_tol(s2.m2, -14.0, 1.e-6);
  // ck_assert_float_eq_tol(s2.m3, -4.0, 1.e-6);
  // ck_assert_float_eq_tol(s2.m4, 8.0, 1.e-6);
  // ck_assert_float_eq_tol(s2.m5, -4.0, 1.e-6);
}
END_TEST

START_TEST(test_screw_invert) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_invert(&s1);
  ck_assert_double_eq_tol(s1.m0, -1.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m1, -2.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m2, -3.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m3, -4.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m4, -5.0, 1.e-9);
  ck_assert_double_eq_tol(s1.m5, -6.0, 1.e-9);
}
END_TEST

START_TEST(test_screwf_invert) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_invert(&s1);
  ck_assert_float_eq_tol(s1.m0, -1.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m1, -2.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m2, -3.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m3, -4.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m4, -5.0, 1.e-6);
  ck_assert_float_eq_tol(s1.m5, -6.0, 1.e-6);
}
END_TEST

START_TEST(test_product_as) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0}};
  ksl_product_as(2.0, &s1, &s2);
  ck_assert_double_eq_tol(s2.m0, 2.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 4.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 6.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 8.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 10.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 12.0, 1.e-9);
}
END_TEST

START_TEST(test_product_asf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0}};
  ksl_product_asf(2.0, &s1, &s2);
  ck_assert_float_eq_tol(s2.m0, 2.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 4.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 6.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 8.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 10.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 12.0, 1.e-6);
}
END_TEST

START_TEST(test_add_ss) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_add_ss(&s1, &s2, &s2);
  ck_assert_double_eq_tol(s2.m0, 11.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 22.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 33.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 44.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 55.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 66.0, 1.e-9);
}
END_TEST

START_TEST(test_add_ssf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_add_ssf(&s1, &s2, &s2);
  ck_assert_float_eq_tol(s2.m0, 11.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 22.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 33.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 44.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 55.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 66.0, 1.e-6);
}
END_TEST

START_TEST(test_xpy_ss) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_xpy_ss(&s1, &s2);
  ck_assert_double_eq_tol(s2.m0, 11.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 22.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 33.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 44.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 55.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 66.0, 1.e-9);
}
END_TEST

START_TEST(test_xpy_ssf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_xpy_ssf(&s1, &s2);
  ck_assert_float_eq_tol(s2.m0, 11.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 22.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 33.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 44.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 55.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 66.0, 1.e-6);
}
END_TEST

START_TEST(test_subtract_ss) {
  ksl_screw_t s1 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_screw_t s2 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_subtract_ss(&s1, &s2, &s2);
  ck_assert_double_eq_tol(s2.m0, 9.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 18.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 27.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 36.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 45.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 54.0, 1.e-9);
}
END_TEST

START_TEST(test_subtract_ssf) {
  ksl_screwf_t s1 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_screwf_t s2 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_subtract_ssf(&s1, &s2, &s2);
  ck_assert_float_eq_tol(s2.m0, 9.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 18.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 27.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 36.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 45.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 54.0, 1.e-6);
}
END_TEST

START_TEST(test_nxpy_ss) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_nxpy_ss(&s1, &s2);
  ck_assert_double_eq_tol(s2.m0, 9.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 18.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 27.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 36.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 45.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 54.0, 1.e-9);
}
END_TEST

START_TEST(test_nxpy_ssf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
  ksl_nxpy_ssf(&s1, &s2);
  ck_assert_float_eq_tol(s2.m0, 9.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 18.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 27.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 36.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 45.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 54.0, 1.e-6);
}
END_TEST

START_TEST(test_add_sst) {
  ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screw_t s2 = {{10.0, 20.0, 30.0, 0.0, 0.0, 0.0}};
  ksl_add_sst(&s1, &s2, &s2);
  ck_assert_double_eq_tol(s2.m0, 11.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m1, 22.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m2, 33.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m3, 4.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m4, 5.0, 1.e-9);
  ck_assert_double_eq_tol(s2.m5, 6.0, 1.e-9);
}
END_TEST

START_TEST(test_add_sstf) {
  ksl_screwf_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_screwf_t s2 = {{10.0, 20.0, 30.0, 0.0, 0.0, 0.0}};
  ksl_add_sstf(&s1, &s2, &s2);
  ck_assert_float_eq_tol(s2.m0, 11.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m1, 22.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m2, 33.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m3, 4.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m4, 5.0, 1.e-6);
  ck_assert_float_eq_tol(s2.m5, 6.0, 1.e-6);
}
END_TEST

START_TEST(test_htx) {
  ksl_SE3_t D;
  D = ksl_SE3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);
  ksl_screw_t s = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
  ksl_htx(&D, &s);
  ck_assert_double_eq_tol(s.m0, 1.0, 1.e-9);
  ck_assert_double_eq_tol(s.m1, 2.0, 1.e-9);
  ck_assert_double_eq_tol(s.m2, 3.0, 1.e-9);
  ck_assert_double_eq_tol(s.m3, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m4, 0.0, 1.e-9);
  ck_assert_double_eq_tol(s.m5, 0.0, 1.e-9);
}
END_TEST

// START_TEST(test_htx) {
//   ksl_SE3_t Di;
//
//   ksl_screw_t s1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}};
//   ksl_screw_t s2 = {{10.0, 20.0, 30.0, 40.0, 50.0, 60.0}};
//   ksl_nxpy_ss(&s1, &s2);
//   ck_assert_double_eq_tol(s2.m0, 9.0, 1.e-9);
//   ck_assert_double_eq_tol(s2.m1, 18.0, 1.e-9);
//   ck_assert_double_eq_tol(s2.m2, 27.0, 1.e-9);
//   ck_assert_double_eq_tol(s2.m3, 36.0, 1.e-9);
//   ck_assert_double_eq_tol(s2.m4, 45.0, 1.e-9);
//   ck_assert_double_eq_tol(s2.m5, 54.0, 1.e-9);
// }
// END_TEST

Suite* screw_suite(void) {
  Suite* s = suite_create("screw");
  TCase* tc_core = tcase_create("core");
  tcase_add_test(tc_core, test_screw_create);
  tcase_add_test(tc_core, test_screwf_create);
  tcase_add_test(tc_core, test_screw_alloc);
  tcase_add_test(tc_core, test_screwf_alloc);
  tcase_add_test(tc_core, test_screw_norm);
  tcase_add_test(tc_core, test_screwf_norm);
  tcase_add_test(tc_core, test_screw_normalize);
  tcase_add_test(tc_core, test_screwf_normalize);
  tcase_add_test(tc_core, test_dot_cs);
  tcase_add_test(tc_core, test_dot_csf);
  tcase_add_test(tc_core, test_axpy_ss);
  tcase_add_test(tc_core, test_axpy_ssf);
  tcase_add_test(tc_core, test_screw_copy);
  tcase_add_test(tc_core, test_screwf_copy);
  tcase_add_test(tc_core, test_screw_inverted);
  tcase_add_test(tc_core, test_screwf_inverted);
  tcase_add_test(tc_core, test_screw_invert);
  tcase_add_test(tc_core, test_screwf_invert);
  tcase_add_test(tc_core, test_product_as);
  tcase_add_test(tc_core, test_product_asf);
  tcase_add_test(tc_core, test_add_ss);
  tcase_add_test(tc_core, test_add_ssf);
  tcase_add_test(tc_core, test_xpy_ss);
  tcase_add_test(tc_core, test_xpy_ssf);
  tcase_add_test(tc_core, test_subtract_ss);
  tcase_add_test(tc_core, test_subtract_ssf);
  tcase_add_test(tc_core, test_nxpy_ss);
  tcase_add_test(tc_core, test_nxpy_ssf);
  tcase_add_test(tc_core, test_add_sst);
  tcase_add_test(tc_core, test_add_sstf);
  tcase_add_test(tc_core, test_htx);
  tcase_add_test(tc_core, test_htxf);
  tcase_add_test(tc_core, test_hty);
  tcase_add_test(tc_core, test_htyf);
  tcase_add_test(tc_core, test_htz);
  tcase_add_test(tc_core, test_htzf);
  tcase_add_test(tc_core, test_hrx);
  tcase_add_test(tc_core, test_hrxf);
  tcase_add_test(tc_core, test_hry);
  tcase_add_test(tc_core, test_hryf);
  tcase_add_test(tc_core, test_hrz);
  tcase_add_test(tc_core, test_hrzf);
  tcase_add_test(tc_core, test_htxinv);
  tcase_add_test(tc_core, test_htxinvf);
  tcase_add_test(tc_core, test_htyinv);
  tcase_add_test(tc_core, test_htyinvf);
  tcase_add_test(tc_core, test_htzinv);
  tcase_add_test(tc_core, test_htzinvf);
  tcase_add_test(tc_core, test_hrxinv);
  tcase_add_test(tc_core, test_hrxinvf);
  tcase_add_test(tc_core, test_hryinv);
  tcase_add_test(tc_core, test_hryinvf);
  tcase_add_test(tc_core, test_hrzinv);
  tcase_add_test(tc_core, test_hrzinvf);
  tcase_add_test(tc_core, test_cross_ss);
  tcase_add_test(tc_core, test_cross_ssf);

  suite_add_tcase(s, tc_core);
  return s;
}

int main(void) {
  dbg = fopen("test_debug.txt", "wa");
  int number_failed;
  Suite* s = screw_suite();
  SRunner* sr = srunner_create(s);
  srunner_set_log(sr, "check_screw.log");
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  fclose(dbg);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
