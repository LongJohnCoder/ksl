#include <config.h>

#include <check.h>
#include <stdio.h>
#include <stdlib.h>

#include "linalg.h"
#include "util.h"

static FILE* dbg;

START_TEST(test_swap_unit_stride) {
  double a[20] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  double b[20] = {11, 12, 13, 14, 15, 6,  7,  8,  9,  10,
                  1,  2,  3,  4,  5,  16, 17, 18, 19, 20};
  ksl_swapArray(5, a, 1, a + 10, 1);
  ck_assert(ksl_allclose(20, a, b));
}
END_TEST

START_TEST(test_swap_nonunit_stride) {
  double a[20] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  double b[20] = {4,  2,  3,  1,  5,  9,  7,  8,  6,  10,
                  14, 12, 13, 11, 15, 19, 17, 18, 16, 20};
  ksl_swapArray(4, a, 5, a + 3, 5);
  ck_assert(ksl_allclose(20, a, b));
}
END_TEST

START_TEST(test_maxIndex) {
  double a[10] = {-21, 45, -622, 431, 622, 0, -1, -633, 1, -99};
  ck_assert(ksl_maxIndex(10, a) == 7);
}
END_TEST

START_TEST(test_linalg_lu_full_rmo) {
  double AO[30] = {3,  1, -3, 3, 0,  2,  2, 2, -6, 3, 1, 7, -4, 0,  0,
                   -3, 1, 3,  8, -4, 12, 0, 4, 4,  6, 1, 3, -3, -2, -2};
  double LU[30] = {12,    4,         8,    0,          4,  -4,        -0.5, 9,
                   6,     3,         3,    0,          0,  1.0 / 3.0, -6,   -4,
                   0,     0,         0.25, -1.0 / 3.0, -1, -6,        -2,   2,
                   -0.25, 1.0 / 3.0, -0.5, 0,          0,  0};

  int pr[5];
  int prB[5] = {3, 1, 2, 4, 0};
  int pc[6];
  int pcB[6] = {2, 5, 0, 3, 4, 1};
  double eps = 1e-12;

  ksl_linalg_lu_full_rmo(5, 6, AO, eps, pr, pc);

  ck_assert(ksl_allclose(30, AO, LU));
  ck_assert(ksl_allclosei(5, pr, prB));
  ck_assert(ksl_allclosei(6, pc, pcB));
}
END_TEST

START_TEST(test_linalg_lu_full_specified_rmo) {
  double AO[30] = {3,  1, -3, 3, 0,  2,  2, 2, -6, 3, 1, 7, -4, 0,  0,
                   -3, 1, 3,  8, -4, 12, 0, 4, 4,  6, 1, 3, -3, -2, -2};
  double LU[30] = {12,    4,         8,    0,          -4, 4,         -0.5, 9,
                   6,     3,         0,    3,          0,  1.0 / 3.0, -6,   -4,
                   0,     0,         0.25, -1.0 / 3.0, -1, -6,        2,    -2,
                   -0.25, 1.0 / 3.0, -0.5, -0,         0,  0};
  int pr[5];
  int prB[5] = {3, 1, 2, 4, 0};
  int pc[6];
  int pcB[6] = {2, 5, 0, 3, 1, 4};
  double eps = 1e-12;

  ksl_linalg_lu_full_specified_rmo(5, 6, AO, eps, pr, pc, 4);

  ck_assert(ksl_allclose(30, AO, LU));
  ck_assert(ksl_allclosei(5, pr, prB));
  ck_assert(ksl_allclosei(6, pc, pcB));
}
END_TEST

START_TEST(test_linalg_lu_rmo) {
  double AO[24] = {12, 4, 8,  0,  4, -4, -6, 7,  2, 3,  1,  2,
                   0,  3, -4, -3, 1, 0,  3,  -2, 6, -3, -2, 1};
  double LU[24] = {12, 4, 8,    0,          4,  -4,        -0.5, 9,
                   6,  3, 3,    0,          0,  1.0 / 3.0, -6,   -4,
                   0,  0, 0.25, -1.0 / 3.0, -1, -6,        -2,   2};

  ksl_linalg_lu_rmo(4, 6, AO);

  ck_assert(ksl_allclose(24, AO, LU));
}
END_TEST

START_TEST(test_linalg_lu_setBMatrix_rmo) {
  double LU[30] = {12,    4,         8,    0,          4,  -4,        -0.5, 9,
                   6,     3,         3,    0,          0,  1.0 / 3.0, -6,   -4,
                   0,     0,         0.25, -1.0 / 3.0, -1, -6,        -2,   2,
                   -0.25, 1.0 / 3.0, -0.5, 0,          0,  0};
  double LUB[30] = {12,    4,          8,    0,  29.0 / 81.0, -38.0 / 81.0,
                    -0.5,  9,          6,    3,  10.0 / 27.0, -1.0 / 27.0,
                    0,     1.0 / 3.0,  -6,   -4, -2.0 / 9.0,  2.0 / 9.0,
                    0.25,  -1.0 / 3.0, -1,   -6, 1.0 / 3.0,   -1.0 / 3.0,
                    -0.25, 1.0 / 3.0,  -0.5, 0,  0,           0};

  ksl_linalg_lu_setBMatrix_rmo(5, 6, 4, LU);

  ck_assert(ksl_allclose(30, LU, LUB));
}
END_TEST

START_TEST(test_linalg_lu_setCMatrix_rmo) {
  double LU[30] = {12,    4,         8,    0,          4,  -4,        -0.5, 9,
                   6,     3,         3,    0,          0,  1.0 / 3.0, -6,   -4,
                   0,     0,         0.25, -1.0 / 3.0, -1, -6,        -2,   2,
                   -0.25, 1.0 / 3.0, -0.5, 0,          0,  0};
  double LUC[30] = {12, 4,  8,  0,         4,  -4,  -0.5, 9, 6,    3,
                    3,  0,  0,  1.0 / 3.0, -6, -4,  0,    0, 0.25, -1.0 / 3.0,
                    -1, -6, -2, 2,         0,  0.5, -0.5, 0, 0,    0};

  ksl_linalg_lu_setCMatrix_rmo(5, 6, 4, LU);

  ck_assert(ksl_allclose(30, LU, LUC));
}
END_TEST

START_TEST(test_linalg_ldlt_rmo) {
  double A[16] = {4, 3,          2,          1,
                  3, 21.0 / 4.0, 7.0 / 2.0,  7.0 / 4.0,
                  2, 7.0 / 2.0,  13.0 / 3.0, 13.0 / 6.0,
                  1, 7.0 / 4.0,  13.0 / 6.0, 25.0 / 12.0};
  double LDLT[16] = {4,         0.75,      0.5, 0.25,      0.75, 3,
                     2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 2,    0.5,
                     0.25,      1.0 / 3.0, 0.5, 1};
  ksl_linalg_ldlt_rmo(A, 4);

  ck_assert(ksl_allclose(16, A, LDLT));
}
END_TEST

START_TEST(test_linalg_cholesky_rmo) {
  double A[16] = {4, 3,          2,          1,
                  3, 21.0 / 4.0, 7.0 / 2.0,  7.0 / 4.0,
                  2, 7.0 / 2.0,  13.0 / 3.0, 13.0 / 6.0,
                  1, 7.0 / 4.0,  13.0 / 6.0, 25.0 / 12.0};
  double UTU[16] = {2,
                    1.5,
                    1,
                    0.5,
                    1.5,
                    sqrt(3),
                    2.0 * sqrt(3) / 3.0,
                    sqrt(3) / 3.0,
                    1,
                    2.0 * sqrt(3) / 3.0,
                    sqrt(2),
                    sqrt(2) / 2.0,
                    0.5,
                    sqrt(3) / 3.0,
                    sqrt(2) / 2.0,
                    1};

  ksl_linalg_cholesky_rmo(A, 4);

  ck_assert(ksl_allclose(16, A, UTU));
}
END_TEST

START_TEST(test_linalg_ldlt_forwardElimination_rmo) {
  double LDLT[16] = {4,         0.75,      0.5, 0.25,      0.75, 3,
                     2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 2,    0.5,
                     0.25,      1.0 / 3.0, 0.5, 1};
  double b[4] = {8, 18, 12, 2};
  double y[4];
  double yB[4] = {8, 12, 0, -4};

  ksl_linalg_ldlt_forwardElimination_rmo(LDLT, b, y, 4);

  ck_assert(ksl_allclose(4, y, yB));
}

END_TEST

START_TEST(test_linalg_ldlt_backwardSubstitution_rmo) {
  double LDLT[16] = {4,         0.75,      0.5, 0.25,      0.75, 3,
                     2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 2,    0.5,
                     0.25,      1.0 / 3.0, 0.5, 1};
  double y[4] = {2, 4, 0, -4};
  double x[4];
  double xB[4] = {-1, 4, 2, -4};

  ksl_linalg_ldlt_backwardSubstitution_rmo(LDLT, y, x, 4);

  ck_assert(ksl_allclose(4, x, xB));
}
END_TEST

START_TEST(test_linalg_cholesky_forwardElimination_rmo) {
  double UTU[16] = {2,
                    1.5,
                    1,
                    0.5,
                    1.5,
                    sqrt(3),
                    2.0 * sqrt(3) / 3.0,
                    sqrt(3) / 3.0,
                    1,
                    2.0 * sqrt(3) / 3.0,
                    sqrt(2),
                    sqrt(2) / 2.0,
                    0.5,
                    sqrt(3) / 3.0,
                    sqrt(2) / 2.0,
                    1};
  double b[4] = {8, 18, 12, 2};
  double y[4];
  double yB[4] = {4, 4 * sqrt(3), 0, -4};

  ksl_linalg_cholesky_forwardElimination_rmo(UTU, b, y, 4);

  ck_assert(ksl_allclose(4, y, yB));
}
END_TEST

START_TEST(test_linalg_cholesky_backwardSubstitution_rmo) {
  double UTU[16] = {2,
                    1.5,
                    1,
                    0.5,
                    1.5,
                    sqrt(3),
                    2.0 * sqrt(3) / 3.0,
                    sqrt(3) / 3.0,
                    1,
                    2.0 * sqrt(3) / 3.0,
                    sqrt(2),
                    sqrt(2) / 2.0,
                    0.5,
                    sqrt(3) / 3.0,
                    sqrt(2) / 2.0,
                    1};
  double y[4] = {4, 4 * sqrt(3), 0, -4};
  double x[4];
  double xB[4] = {-1, 4, 2, -4};

  ksl_linalg_cholesky_backwardSubstitution_rmo(UTU, y, x, 4);

  ck_assert(ksl_allclose(4, x, xB));
}
END_TEST

START_TEST(test_linalg_ldlt_solve_rmo) {
  double LDLT[16] = {4,         0.75,      0.5, 0.25,      0.75, 3,
                     2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 2,    0.5,
                     0.25,      1.0 / 3.0, 0.5, 1};
  double b[4] = {8, 18, 12, 2};
  double x[4];
  double xB[4] = {-1, 4, 2, -4};

  ksl_linalg_ldlt_solve_rmo(LDLT, b, x, 4);

  ck_assert(ksl_allclose(4, x, xB));
}
END_TEST

START_TEST(test_linalg_cholesky_solve_rmo) {
  double UTU[16] = {2,
                    1.5,
                    1,
                    0.5,
                    1.5,
                    sqrt(3),
                    2.0 * sqrt(3) / 3.0,
                    sqrt(3) / 3.0,
                    1,
                    2.0 * sqrt(3) / 3.0,
                    sqrt(2),
                    sqrt(2) / 2.0,
                    0.5,
                    sqrt(3) / 3.0,
                    sqrt(2) / 2.0,
                    1};
  double b[4] = {8, 18, 12, 2};
  double x[4];
  double xB[4] = {-1, 4, 2, -4};

  ksl_linalg_cholesky_solve_rmo(UTU, b, x, 4);

  ck_assert(ksl_allclose(4, x, xB));
}
END_TEST

START_TEST(test_linalg_symmetricMatrixInverse_rmo) {
  double A[16] = {4, 3,          2,          1,
                  3, 21.0 / 4.0, 7.0 / 2.0,  7.0 / 4.0,
                  2, 7.0 / 2.0,  13.0 / 3.0, 13.0 / 6.0,
                  1, 7.0 / 4.0,  13.0 / 6.0, 25.0 / 12.0};
  double AI[16] = {7.0 / 16.0, -0.25, 0,    0,          -0.25,     5.0 / 9.0,
                   -1.0 / 3.0, 0,     0,    -1.0 / 3.0, 3.0 / 4.0, -0.5,
                   0,          0,     -0.5, 1};

  ksl_linalg_symmetricMatrixInverse_rmo(A, 4);

  ck_assert(ksl_allclose(16, A, AI));
}
END_TEST

START_TEST(test_linalg_lu_full_cmo) {
  double AO[30] = {3, 2, -4, 8, 6,  1, 2, 0, -4, 1,  -3, -6, 0, 12, 3,
                   3, 3, -3, 0, -3, 0, 1, 1, 4,  -2, 2,  7,  3, 4,  -2};
  //  3,  1, -3, 3, 0,  2,
  //  2,  2, -6, 3, 1,  7,
  // -4,  0,  0,-3, 1,  3,
  //  8, -4, 12, 0, 4,  4,
  //  6, 1, 3, -3, -2, -2};
  double LU[30] = {12,         -0.5,      0,  0.25, -0.25, 4,  9,    1.0 / 3.0,
                   -1.0 / 3.0, 1.0 / 3.0, 8,  6,    -6,    -1, -0.5, 0,
                   3,          -4,        -6, 0,    4,     3,  0,    -2,
                   0,          -4,        0,  0,    2,     0};
  //   12,          4,    8,  0,  4, -4,
  // -0.5,          9,    6,  3,  3,  0,
  //    0,  1.0 / 3.0,   -6, -4,  0,  0,
  // 0.25, -1.0 / 3.0,   -1, -6, -2,  2,
  //-0.25,  1.0 / 3.0, -0.5,  0,  0,  0

  int pr[5];
  int prB[5] = {3, 1, 2, 4, 0};
  int pc[6];
  int pcB[6] = {2, 5, 0, 3, 4, 1};
  double eps = 1e-12;

  ksl_linalg_lu_full_cmo(5, 6, AO, eps, pr, pc);

  ck_assert(ksl_allclose(30, AO, LU));
  ck_assert(ksl_allclosei(5, pr, prB));
  ck_assert(ksl_allclosei(6, pc, pcB));
}
END_TEST

START_TEST(test_linalg_lu_full_specified_cmo) {
  double AO[30] = {3, 2, -4, 8, 6,  1, 2, 0, -4, 1,  -3, -6, 0, 12, 3,
                   3, 3, -3, 0, -3, 0, 1, 1, 4,  -2, 2,  7,  3, 4,  -2};
  //  3,  1, -3, 3, 0,  2,
  //  2,  2, -6, 3, 1,  7,
  // -4,  0,  0,-3, 1,  3,
  //  8, -4, 12, 0, 4,  4,
  //  6, 1, 3, -3, -2, -2};
  double LU[30] = {12,         -0.5,      0,  0.25, -0.25, 4,  9,    1.0 / 3.0,
                   -1.0 / 3.0, 1.0 / 3.0, 8,  6,    -6,    -1, -0.5, 0,
                   3,          -4,        -6, 0,    -4,    0,  0,    2,
                   0,          4,         3,  0,    -2,    0};
  //   12,          4,    8,  0, -4,  4,
  // -0.5,          9,    6,  3,  0,  3,
  //    0,  1.0 / 3.0,   -6, -4,  0,  0,
  // 0.25, -1.0 / 3.0,   -1, -6,  2, -2,
  //-0.25,  1.0 / 3.0, -0.5,  0,  0,  0
  int pr[5];
  int prB[5] = {3, 1, 2, 4, 0};
  int pc[6];
  int pcB[6] = {2, 5, 0, 3, 1, 4};
  double eps = 1e-12;

  ksl_linalg_lu_full_specified_cmo(5, 6, AO, eps, pr, pc, 4);

  ck_assert(ksl_allclose(30, AO, LU));
  ck_assert(ksl_allclosei(5, pr, prB));
  ck_assert(ksl_allclosei(6, pc, pcB));
}

END_TEST

START_TEST(test_linalg_lu_cmo) {
  double AO[24] = {12, -6, 0,  3,  4, 7, 3, -2, 8,  2, -4, 6,
                   0,  3,  -3, -3, 4, 1, 1, -2, -4, 2, 0,  1};
  // 12,  4,  8,  0,  4, -4,
  // -6,  7,  2,  3,  1,  2,
  //  0,  3, -4, -3,  1,  0,
  //  3, -2,  6, -3, -2,  1
  double LU[24] = {12, -0.5, 0,  0.25, 4,  9, 1.0 / 3.0, -1.0 / 3.0,
                   8,  6,    -6, -1,   0,  3, -4,        -6,
                   4,  3,    0,  -2,   -4, 0, 0,         2};
  //   12,          4,  8,  0,  4, -4,
  // -0.5,          9,  6,  3,  3,  0,
  //    0,  1.0 / 3.0, -6, -4,  0,  0,
  // 0.25, -1.0 / 3.0, -1, -6, -2,  2

  ksl_linalg_lu_cmo(4, 6, AO);

  ck_assert(ksl_allclose(24, AO, LU));
}
END_TEST

START_TEST(test_linalg_lu_setBMatrix_cmo) {
  double LU[30] = {12,         -0.5,      0,  0.25, -0.25, 4,  9,    1.0 / 3.0,
                   -1.0 / 3.0, 1.0 / 3.0, 8,  6,    -6,    -1, -0.5, 0,
                   3,          -4,        -6, 0,    4,     3,  0,    -2,
                   0,          -4,        0,  0,    2,     0};
  //   12,          4,    8,  0,  4, -4,
  // -0.5,          9,    6,  3,  3,  0,
  //    0,  1.0 / 3.0,   -6, -4,  0,  0,
  // 0.25, -1.0 / 3.0,   -1, -6, -2,  2,
  //-0.25,  1.0 / 3.0, -0.5,  0,  0,  0
  double LUB[30] = {12,          -0.5,
                    0,           0.25,
                    -0.25,       4,
                    9,           1.0 / 3.0,
                    -1.0 / 3.0,  1.0 / 3.0,
                    8,           6,
                    -6,          -1,
                    -0.5,        0,
                    3,           -4,
                    -6,          0,
                    29.0 / 81.0, 10.0 / 27.0,
                    -2.0 / 9.0,  1.0 / 3.0,
                    0,           -38.0 / 81.0,
                    -1.0 / 27.0, 2.0 / 9.0,
                    -1.0 / 3.0,  0};
  //    12,          4,    8,  0, 29.0 / 81.0, -38.0 / 81.0,
  //  -0.5,          9,    6,  3,  0.0 / 27.0,  -1.0 / 27.0,
  //     0,  1.0 / 3.0,   -6, -4,  -2.0 / 9.0,    2.0 / 9.0,
  //  0.25, -1.0 / 3.0,   -1, -6,   1.0 / 3.0,   -1.0 / 3.0,
  // -0.25,  1.0 / 3.0, -0.5,  0,           0,            0

  ksl_linalg_lu_setBMatrix_cmo(5, 6, 4, LU);

  ck_assert(ksl_allclose(30, LU, LUB));
}
END_TEST

START_TEST(test_linalg_lu_setCMatrix_cmo) {

  double LU[30] = {12,         -0.5,      0,  0.25, -0.25, 4,  9,    1.0 / 3.0,
                   -1.0 / 3.0, 1.0 / 3.0, 8,  6,    -6,    -1, -0.5, 0,
                   3,          -4,        -6, 0,    4,     3,  0,    -2,
                   0,          -4,        0,  0,    2,     0};
  //   12,          4,    8,  0,  4, -4,
  // -0.5,          9,    6,  3,  3,  0,
  //    0,  1.0 / 3.0,   -6, -4,  0,  0,
  // 0.25, -1.0 / 3.0,   -1, -6, -2,  2,
  //-0.25,  1.0 / 3.0, -0.5,  0,  0,  0
  double LUC[30] = {12, -0.5, 0,  0.25, 0,    4,  9, 1.0 / 3.0, -1.0 / 3.0, 0.5,
                    8,  6,    -6, -1,   -0.5, 0,  3, -4,        -6,         0,
                    4,  3,    0,  -2,   0,    -4, 0, 0,         2,          0};

  ksl_linalg_lu_setCMatrix_cmo(5, 6, 4, LU);

  ck_assert(ksl_allclose(30, LU, LUC));
}
END_TEST

START_TEST(test_linalg_ldlt_cmo) {
  double A[16] = {4, 3,          2,          1,
                  3, 21.0 / 4.0, 7.0 / 2.0,  7.0 / 4.0,
                  2, 7.0 / 2.0,  13.0 / 3.0, 13.0 / 6.0,
                  1, 7.0 / 4.0,  13.0 / 6.0, 25.0 / 12.0};
  double LDLT[16] = {4,         0.75,      0.5, 0.25,      0.75, 3,
                     2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 2,    0.5,
                     0.25,      1.0 / 3.0, 0.5, 1};

  ksl_linalg_ldlt_cmo(A, 4);

  ck_assert(ksl_allclose(16, A, LDLT));
}
END_TEST

START_TEST(test_linalg_cholesky_cmo) {
  double A[16] = {4, 3,          2,          1,
                  3, 21.0 / 4.0, 7.0 / 2.0,  7.0 / 4.0,
                  2, 7.0 / 2.0,  13.0 / 3.0, 13.0 / 6.0,
                  1, 7.0 / 4.0,  13.0 / 6.0, 25.0 / 12.0};
  double UTU[16] = {2,
                    1.5,
                    1,
                    0.5,
                    1.5,
                    sqrt(3),
                    2.0 * sqrt(3) / 3.0,
                    sqrt(3) / 3.0,
                    1,
                    2.0 * sqrt(3) / 3.0,
                    sqrt(2),
                    sqrt(2) / 2.0,
                    0.5,
                    sqrt(3) / 3.0,
                    sqrt(2) / 2.0,
                    1};

  ksl_linalg_cholesky_cmo(A, 4);

  ck_assert(ksl_allclose(16, A, UTU));
}
END_TEST

START_TEST(test_linalg_ldlt_forwardElimination_cmo) {
  double LDLT[16] = {4,         0.75,      0.5, 0.25,      0.75, 3,
                     2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 2,    0.5,
                     0.25,      1.0 / 3.0, 0.5, 1};
  double b[4] = {8, 18, 12, 2};
  double y[4];
  double yB[4] = {8, 12, 0, -4};

  ksl_linalg_ldlt_forwardElimination_cmo(LDLT, b, y, 4);

  ck_assert(ksl_allclose(4, y, yB));
}
END_TEST

START_TEST(test_linalg_ldlt_backwardSubstitution_cmo) {
  double LDLT[16] = {4,         0.75,      0.5, 0.25,      0.75, 3,
                     2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 2,    0.5,
                     0.25,      1.0 / 3.0, 0.5, 1};
  double y[4] = {2, 4, 0, -4};
  double x[4];
  double xB[4] = {-1, 4, 2, -4};

  ksl_linalg_ldlt_backwardSubstitution_cmo(LDLT, y, x, 4);

  ck_assert(ksl_allclose(4, x, xB));
}
END_TEST

START_TEST(test_linalg_cholesky_forwardElimination_cmo) {
  double UTU[16] = {2,
                    1.5,
                    1,
                    0.5,
                    1.5,
                    sqrt(3),
                    2.0 * sqrt(3) / 3.0,
                    sqrt(3) / 3.0,
                    1,
                    2.0 * sqrt(3) / 3.0,
                    sqrt(2),
                    sqrt(2) / 2.0,
                    0.5,
                    sqrt(3) / 3.0,
                    sqrt(2) / 2.0,
                    1};
  double b[4] = {8, 18, 12, 2};
  double y[4];
  double yB[4] = {4, 4 * sqrt(3), 0, -4};

  ksl_linalg_cholesky_forwardElimination_cmo(UTU, b, y, 4);

  ck_assert(ksl_allclose(4, y, yB));
}
END_TEST

START_TEST(test_linalg_cholesky_backwardSubstitution_cmo) {
  double UTU[16] = {2,
                    1.5,
                    1,
                    0.5,
                    1.5,
                    sqrt(3),
                    2.0 * sqrt(3) / 3.0,
                    sqrt(3) / 3.0,
                    1,
                    2.0 * sqrt(3) / 3.0,
                    sqrt(2),
                    sqrt(2) / 2.0,
                    0.5,
                    sqrt(3) / 3.0,
                    sqrt(2) / 2.0,
                    1};
  double y[4] = {4, 4 * sqrt(3), 0, -4};
  double x[4];
  double xB[4] = {-1, 4, 2, -4};

  ksl_linalg_cholesky_backwardSubstitution_cmo(UTU, y, x, 4);

  ck_assert(ksl_allclose(4, x, xB));
}
END_TEST

START_TEST(test_linalg_ldlt_solve_cmo) {
  double LDLT[16] = {4,         0.75,      0.5, 0.25,      0.75, 3,
                     2.0 / 3.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 2,    0.5,
                     0.25,      1.0 / 3.0, 0.5, 1};
  double b[4] = {8, 18, 12, 2};
  double x[4];
  double xB[4] = {-1, 4, 2, -4};

  ksl_linalg_ldlt_solve_cmo(LDLT, b, x, 4);

  ck_assert(ksl_allclose(4, x, xB));
}
END_TEST

START_TEST(test_linalg_cholesky_solve_cmo) {
  double UTU[16] = {2,
                    1.5,
                    1,
                    0.5,
                    1.5,
                    sqrt(3),
                    2.0 * sqrt(3) / 3.0,
                    sqrt(3) / 3.0,
                    1,
                    2.0 * sqrt(3) / 3.0,
                    sqrt(2),
                    sqrt(2) / 2.0,
                    0.5,
                    sqrt(3) / 3.0,
                    sqrt(2) / 2.0,
                    1};
  double b[4] = {8, 18, 12, 2};
  double x[4];
  double xB[4] = {-1, 4, 2, -4};

  ksl_linalg_cholesky_solve_cmo(UTU, b, x, 4);

  ck_assert(ksl_allclose(4, x, xB));
}
END_TEST

START_TEST(test_linalg_symmetricMatrixInverse_cmo) {
  double A[16] = {4, 3,          2,          1,
                  3, 21.0 / 4.0, 7.0 / 2.0,  7.0 / 4.0,
                  2, 7.0 / 2.0,  13.0 / 3.0, 13.0 / 6.0,
                  1, 7.0 / 4.0,  13.0 / 6.0, 25.0 / 12.0};
  double AI[16] = {7.0 / 16.0, -0.25, 0,    0,          -0.25,     5.0 / 9.0,
                   -1.0 / 3.0, 0,     0,    -1.0 / 3.0, 3.0 / 4.0, -0.5,
                   0,          0,     -0.5, 1};

  ksl_linalg_symmetricMatrixInverse_cmo(A, 4);

  ck_assert(ksl_allclose(16, A, AI));
}
END_TEST

Suite* linalg_suite(void) {
  Suite* s = suite_create("linalg");
  TCase* tc_core = tcase_create("core");
  tcase_add_test(tc_core, test_swap_unit_stride);
  tcase_add_test(tc_core, test_swap_nonunit_stride);
  tcase_add_test(tc_core, test_maxIndex);
  tcase_add_test(tc_core, test_linalg_lu_full_rmo);
  tcase_add_test(tc_core, test_linalg_lu_full_specified_rmo);
  tcase_add_test(tc_core, test_linalg_lu_rmo);
  tcase_add_test(tc_core, test_linalg_lu_setBMatrix_rmo);
  tcase_add_test(tc_core, test_linalg_lu_setCMatrix_rmo);
  tcase_add_test(tc_core, test_linalg_ldlt_rmo);
  tcase_add_test(tc_core, test_linalg_cholesky_rmo);
  tcase_add_test(tc_core, test_linalg_ldlt_forwardElimination_rmo);
  tcase_add_test(tc_core, test_linalg_ldlt_backwardSubstitution_rmo);
  tcase_add_test(tc_core, test_linalg_cholesky_forwardElimination_rmo);
  tcase_add_test(tc_core, test_linalg_cholesky_backwardSubstitution_rmo);
  tcase_add_test(tc_core, test_linalg_ldlt_solve_rmo);
  tcase_add_test(tc_core, test_linalg_cholesky_solve_rmo);
  tcase_add_test(tc_core, test_linalg_symmetricMatrixInverse_rmo);
  tcase_add_test(tc_core, test_linalg_lu_full_cmo);
  tcase_add_test(tc_core, test_linalg_lu_full_specified_cmo);
  tcase_add_test(tc_core, test_linalg_lu_cmo);
  tcase_add_test(tc_core, test_linalg_lu_setBMatrix_cmo);
  tcase_add_test(tc_core, test_linalg_lu_setCMatrix_cmo);
  tcase_add_test(tc_core, test_linalg_ldlt_cmo);
  tcase_add_test(tc_core, test_linalg_cholesky_cmo);
  tcase_add_test(tc_core, test_linalg_ldlt_forwardElimination_cmo);
  tcase_add_test(tc_core, test_linalg_ldlt_backwardSubstitution_cmo);
  tcase_add_test(tc_core, test_linalg_cholesky_forwardElimination_cmo);
  tcase_add_test(tc_core, test_linalg_cholesky_backwardSubstitution_cmo);
  tcase_add_test(tc_core, test_linalg_ldlt_solve_cmo);
  tcase_add_test(tc_core, test_linalg_cholesky_solve_cmo);
  tcase_add_test(tc_core, test_linalg_symmetricMatrixInverse_cmo);
  suite_add_tcase(s, tc_core);
  return s;
}

int main(void) {
  dbg = fopen("test_debug.txt", "wa");
  int number_failed;
  Suite* s = linalg_suite();
  SRunner* sr = srunner_create(s);
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  fclose(dbg);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
