#include <config.h>

#include <check.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "print.h"

START_TEST(test_matrix_SE3) {
  ksl_SE3_t* D = ksl_SE3_alloc(3);
  ck_assert_ptr_nonnull(D);

  D[1] = ksl_SE3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  double k = 1;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      ck_assert_double_eq(D[1].as_array[j][i], k++);
    }
  }
  free(D);
}
END_TEST

START_TEST(test_matrix_SE3f) {
  ksl_SE3f_t* D = ksl_SE3f_alloc(3);
  ck_assert_ptr_nonnull(D);

  D[1] =
    ksl_SE3f(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);

  float k = 1;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      ck_assert_float_eq(D[1].as_array[j][i], k++);
    }
  }
  free(D);
}
END_TEST

START_TEST(test_matrix_mat3x3) {
  ksl_mat3x3_t* D = ksl_mat3x3_alloc(3);
  ck_assert_ptr_nonnull(D);

  D[1] = ksl_mat3x3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

  double k = 1;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ck_assert_double_eq(D[1].as_array[j][i], k++);
    }
  }
  free(D);
}
END_TEST

START_TEST(test_matrix_mat3x3f) {
  ksl_mat3x3f_t* D = ksl_mat3x3f_alloc(3);
  ck_assert_ptr_nonnull(D);

  D[1] = ksl_mat3x3f(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

  float k = 1;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ck_assert_float_eq(D[1].as_array[j][i], k++);
    }
  }
  free(D);
}
END_TEST

START_TEST(test_matrix_mat4x4) {
  ksl_mat4x4_t* D = ksl_mat4x4_alloc(3);
  ck_assert_ptr_nonnull(D);

  D[1] = ksl_mat4x4(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
                    12.0, 13.0, 14.0, 15.0, 16.0);

  double k = 1;
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      ck_assert_double_eq(D[1].as_array[j][i], k++);
    }
  }
  free(D);
}
END_TEST

START_TEST(test_matrix_mat4x4f) {
  ksl_mat4x4f_t* D = ksl_mat4x4f_alloc(3);
  ck_assert_ptr_nonnull(D);

  D[1] = ksl_mat4x4f(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
                     12.0, 13.0, 14.0, 15.0, 16.0);

  float k = 1;
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      ck_assert_float_eq(D[1].as_array[j][i], k++);
    }
  }
  free(D);
}
END_TEST

START_TEST(test_matrix_SE3_alloc) {
  ksl_SE3_t* D = ksl_SE3_alloc(1);
  ck_assert_ptr_nonnull(D);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert_double_eq(D[0].as_array[j][i], 1.0);
      } else {
        ck_assert_double_eq(D[0].as_array[j][i], 0.0);
      }
    }
  }
  free(D);
}
END_TEST

START_TEST(test_matrix_SE3f_alloc) {
  const int n = 2;
  ksl_SE3f_t* Df = ksl_SE3f_alloc(n);
  ck_assert_ptr_nonnull(Df);
  for(int k = 0; k < n; k++) {
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 4; j++) {
        if(i == j) {
          ck_assert_double_eq(Df[k].as_array[j][i], 1.0);
        } else {
          ck_assert_double_eq(Df[k].as_array[j][i], 0.0);
        }
      }
    }
  }
  free(Df);
}
END_TEST

START_TEST(test_matrix_mat3x3_alloc) {
  const int n = 3;
  ksl_mat3x3_t* R = ksl_mat3x3_alloc(n);
  ck_assert_ptr_nonnull(R);
  for(int k = 0; k < n; k++) {
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        if(i == j) {
          ck_assert_double_eq(R[k].as_array[j][i], 1.0);
        } else {
          ck_assert_double_eq(R[k].as_array[j][i], 0.0);
        }
      }
    }
  }
  free(R);
}
END_TEST

START_TEST(test_matrix_mat3x3f_alloc) {
  const int n = 4;
  ksl_mat3x3f_t* Rf = ksl_mat3x3f_alloc(n);
  ck_assert_ptr_nonnull(Rf);
  for(int k = 0; k < n; k++) {
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        if(i == j) {
          ck_assert_double_eq(Rf[k].as_array[j][i], 1.0);
        } else {
          ck_assert_double_eq(Rf[k].as_array[j][i], 0.0);
        }
      }
    }
  }
  free(Rf);
}
END_TEST

START_TEST(test_matrix_mat4x4_alloc) {
  const int n = 5;
  ksl_mat4x4_t* M = ksl_mat4x4_alloc(n);
  ck_assert_ptr_nonnull(M);
  for(int k = 0; k < n; k++) {
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        if(i == j) {
          ck_assert_double_eq(M[k].as_array[j][i], 1.0);
        } else {
          ck_assert_double_eq(M[k].as_array[j][i], 0.0);
        }
      }
    }
  }
  free(M);
}
END_TEST

START_TEST(test_matrix_mat4x4f_alloc) {
  const int n = 6;
  ksl_mat4x4f_t* Mf = ksl_mat4x4f_alloc(n);
  ck_assert_ptr_nonnull(Mf);
  for(int k = 0; k < n; k++) {
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        if(i == j) {
          ck_assert_double_eq(Mf[k].as_array[j][i], 1.0);
        } else {
          ck_assert_double_eq(Mf[k].as_array[j][i], 0.0);
        }
      }
    }
  }
  free(Mf);
}
END_TEST

START_TEST(test_matrix_SE3toMat4x4) {
  ksl_SE3_t d = {
    {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}};
  ksl_mat4x4_t m;
  ksl_SE3_toMat4x4(&d, &m);
  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < 3; i++) {
      ck_assert_double_eq(m.as_array[j][i], d.as_array[j][i]);
    }
    if(j != 3) {
      ck_assert_double_eq(m.as_array[j][3], 0.0);
    } else {
      ck_assert_double_eq(m.as_array[j][3], 1.0);
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3ftoMat4x4f) {
  ksl_SE3f_t d = {
    {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}};
  ksl_mat4x4f_t m;
  ksl_SE3f_toMat4x4f(&d, &m);
  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < 3; i++) {
      ck_assert_float_eq(m.as_array[j][i], d.as_array[j][i]);
    }
    if(j != 3) {
      ck_assert_float_eq(m.as_array[j][3], 0.0);
    } else {
      ck_assert_float_eq(m.as_array[j][3], 1.0);
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3toMat4x4f) {
  ksl_SE3_t d = {
    {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}};
  ksl_mat4x4f_t m;
  ksl_SE3_toMat4x4f(&d, &m);
  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < 3; i++) {
      ck_assert_float_eq(m.as_array[j][i], (float)d.as_array[j][i]);
    }
    if(j != 3) {
      ck_assert_float_eq(m.as_array[j][3], 0.0);
    } else {
      ck_assert_float_eq(m.as_array[j][3], 1.0);
    }
  }
}
END_TEST

START_TEST(test_matrix_dc) {
  double dc[2];
  ksl_dc(0.2, dc);
  ck_assert_double_eq(dc[0], 0.19866933079506122);
  ck_assert_double_eq(dc[1], 0.9800665778412416);
  ksl_dc(5.2, dc);
  ck_assert_double_eq(dc[0], -0.8834546557201531);
  ck_assert_double_eq(dc[1], 0.4685166713003771);
}
END_TEST

START_TEST(test_matrix_SE3_setIdentity) {
  ksl_SE3_t d = {
    {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}};
  ksl_SE3_setIdentity(&d);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert_double_eq(d.as_array[j][i], 1.0);
      } else {
        ck_assert_double_eq(d.as_array[j][i], 0.0);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3f_setIdentity) {
  ksl_SE3f_t d = {
    {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}};
  ksl_SE3f_setIdentity(&d);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert_float_eq(d.as_array[j][i], 1.0);
      } else {
        ck_assert_float_eq(d.as_array[j][i], 0.0);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_setIdentity) {
  ksl_mat3x3_t r = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}};
  ksl_mat3x3_setIdentity(&r);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(i == j) {
        ck_assert_double_eq(r.as_array[j][i], 1.0);
      } else {
        ck_assert_double_eq(r.as_array[j][i], 0.0);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3f_setIdentity) {
  ksl_mat3x3f_t rf = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}};
  ksl_mat3x3f_setIdentity(&rf);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(i == j) {
        ck_assert_float_eq(rf.as_array[j][i], 1.0);
      } else {
        ck_assert_float_eq(rf.as_array[j][i], 0.0);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_mat4x4_setIdentity) {
  ksl_mat4x4_t m = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
                     12.0, 13.0, 14.0, 15.0, 16.0}};
  ksl_mat4x4_setIdentity(&m);
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert_double_eq(m.as_array[j][i], 1.0);
      } else {
        ck_assert_double_eq(m.as_array[j][i], 0.0);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_mat4x4f_setIdentity) {
  ksl_mat4x4f_t m = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
                      12.0, 13.0, 14.0, 15.0, 16.0}};
  ksl_mat4x4f_setIdentity(&m);
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert_float_eq(m.as_array[j][i], 1.0);
      } else {
        ck_assert_float_eq(m.as_array[j][i], 0.0);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_setAndGet) {
  ksl_mat3x3_t m;
  double k = 0.0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ksl_mat3x3_set(&m, i, j, k);
      double result = ksl_mat3x3_get(&m, i, j);
      ck_assert_double_eq(result, k);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3f_setAndGet) {
  ksl_mat3x3f_t m;
  float k = 0.0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ksl_mat3x3f_set(&m, i, j, k);
      float result = ksl_mat3x3f_get(&m, i, j);
      ck_assert_float_eq(result, k);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat4x4_setAndGet) {
  ksl_mat4x4_t m;
  double k = 0.0;
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      ksl_mat4x4_set(&m, i, j, k);
      double result = ksl_mat4x4_get(&m, i, j);
      ck_assert_double_eq(result, k);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat4x4f_setAndGet) {
  ksl_mat4x4f_t m;
  float k = 0.0;
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      ksl_mat4x4f_set(&m, i, j, k);
      float result = ksl_mat4x4f_get(&m, i, j);
      ck_assert_float_eq(result, k);
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3_setAndGet) {
  ksl_SE3_t m;
  double k = 0.0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      ksl_SE3_set(&m, i, j, k);
      double result = ksl_SE3_get(&m, i, j);
      ck_assert_double_eq(result, k);
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3f_setAndGet) {
  ksl_SE3f_t m;
  float k = 0.0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      ksl_SE3f_set(&m, i, j, k);
      float result = ksl_SE3f_get(&m, i, j);
      ck_assert_float_eq(result, k);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_setFromVectors) {
  ksl_mat3x3_t m;
  ksl_vec3_t v[3] = {{{1.0, 2.0, 3.0}}, {{4.0, 5.0, 6.0}}, {{7.0, 8.0, 9.0}}};
  ksl_mat3x3_setFromVectors(&m, &v[0], &v[1], &v[2]);

  for(int j = 0; j < 3; j++) {
    for(int i = 0; i < 3; i++) {
      ck_assert_double_eq(m.as_array[j][i], v[j].at[i]);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3f_setFromVectors) {
  ksl_mat3x3f_t m;
  ksl_vec3f_t v[3] = {{{1.0, 2.0, 3.0}}, {{4.0, 5.0, 6.0}}, {{7.0, 8.0, 9.0}}};
  ksl_mat3x3f_setFromVectors(&m, &v[0], &v[1], &v[2]);

  for(int j = 0; j < 3; j++) {
    for(int i = 0; i < 3; i++) {
      ck_assert_float_eq(m.as_array[j][i], v[j].at[i]);
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3_getTranslation) {
  ksl_SE3_t d = {
    {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}};
  ksl_vec3_t v;
  ksl_SE3_getTranslation(&d, &v);

  for(int i = 0; i < 3; i++) {
    ck_assert_double_eq(d.t.at[i], v.at[i]);
  }
}
END_TEST

START_TEST(test_matrix_SE3f_getTranslation) {
  ksl_SE3f_t d = {
    {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0}};
  ksl_vec3f_t v;
  ksl_SE3f_getTranslation(&d, &v);

  for(int i = 0; i < 3; i++) {
    ck_assert_float_eq(d.t.at[i], v.at[i]);
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_copy) {
  ksl_mat3x3_t r1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}};
  ksl_mat3x3_t r2;
  ksl_mat3x3_copy(&r1, &r2);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq(r1.at[i], r2.at[i]);
  }
}
END_TEST

START_TEST(test_matrix_mat3x3f_copy) {
  ksl_mat3x3f_t r1 = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}};
  ksl_mat3x3f_t r2;
  ksl_mat3x3f_copy(&r1, &r2);

  for(int i = 0; i < 9; i++) {
    ck_assert_float_eq(r1.at[i], r2.at[i]);
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_invert) {
  ksl_mat3x3_t r1 = ksl_mat3x3(13.0, 8.0, 9.0, 3.0, 5.0, 1.0, 1.0, 0.5, 3.0);

  // FILE* f = fopen("matrix_inv.txt", "w");
  ksl_mat3x3_invert(&r1);

  // ksl_mat3x3_print(f, &r1);
  // fclose(f);

  ksl_mat3x3_t r2 =
    ksl_mat3x3(0.15591397849462368, -0.20967741935483872, -0.39784946236559143,
               -0.08602150537634409, 0.32258064516129037, 0.1505376344086022,
               -0.037634408602150546, 0.016129032258064523, 0.4408602150537635);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq_tol(r1.at[i], r2.at[i], 1e-9);
  }
}
END_TEST

START_TEST(test_matrix_mat3x3f_invert) {
  ksl_mat3x3f_t r1 = ksl_mat3x3f(13.0, 8.0, 9.0, 3.0, 5.0, 1.0, 1.0, 0.5, 3.0);

  // FILE* f = fopen("matrix_inv.txt", "w");
  ksl_mat3x3f_invert(&r1);

  // ksl_mat3x3_print(f, &r1);
  // fclose(f);

  ksl_mat3x3f_t r2 = ksl_mat3x3f(
    0.15591397849462368, -0.20967741935483872, -0.39784946236559143,
    -0.08602150537634409, 0.32258064516129037, 0.1505376344086022,
    -0.037634408602150546, 0.016129032258064523, 0.4408602150537635);

  for(int i = 0; i < 9; i++) {
    ck_assert_float_eq_tol(r1.at[i], r2.at[i], 1e-6);
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_inverted) {
  ksl_mat3x3_t r1 = ksl_mat3x3(13.0, 8.0, 9.0, 3.0, 5.0, 1.0, 1.0, 0.5, 3.0);
  ksl_mat3x3_t r2;
  ksl_mat3x3_inverted(&r1, &r2);
  ksl_mat3x3_invert(&r1);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq(r1.at[i], r2.at[i]);
  }
}
END_TEST

START_TEST(test_matrix_mat3x3f_inverted) {
  ksl_mat3x3f_t r1 = ksl_mat3x3f(13.0, 8.0, 9.0, 3.0, 5.0, 1.0, 1.0, 0.5, 3.0);
  ksl_mat3x3f_t r2;
  ksl_mat3x3f_inverted(&r1, &r2);
  ksl_mat3x3f_invert(&r1);

  for(int i = 0; i < 9; i++) {
    ck_assert(fabs(r1.at[i] - r2.at[i]) < 1e-8);
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_transpose) {
  ksl_mat3x3_t r1 = ksl_mat3x3(13.0, 8.0, 9.0, 3.0, 5.0, 1.0, 1.0, 0.5, 3.0);
  ksl_mat3x3_t r2;
  ksl_mat3x3_copy(&r1, &r2);
  ksl_mat3x3_transpose(&r1);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ck_assert_double_eq(r1.as_array[i][j], r2.as_array[j][i]);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3f_transpose) {
  ksl_mat3x3f_t r1 = ksl_mat3x3f(13.0, 8.0, 9.0, 3.0, 5.0, 1.0, 1.0, 0.5, 3.0);
  ksl_mat3x3f_t r2;
  ksl_mat3x3f_copy(&r1, &r2);
  ksl_mat3x3f_transpose(&r1);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ck_assert_float_eq(r1.as_array[i][j], r2.as_array[j][i]);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_transposed) {
  ksl_mat3x3_t r1 = ksl_mat3x3(13.0, 8.0, 9.0, 3.0, 5.0, 1.0, 1.0, 0.5, 3.0);
  ksl_mat3x3_t r2;
  ksl_mat3x3_transposed(&r1, &r2);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ck_assert_double_eq(r1.as_array[i][j], r2.as_array[j][i]);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3f_transposed) {
  ksl_mat3x3f_t r1 = ksl_mat3x3f(13.0, 8.0, 9.0, 3.0, 5.0, 1.0, 1.0, 0.5, 3.0);
  ksl_mat3x3f_t r2;
  ksl_mat3x3f_transposed(&r1, &r2);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      ck_assert_float_eq(r1.as_array[i][j], r2.as_array[j][i]);
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3_copy) {
  ksl_SE3_t d1 =
    ksl_SE3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);
  ksl_SE3_t d2;
  ksl_SE3_copy(&d1, &d2);
  for(int i = 0; i < 12; i++) {
    ck_assert_double_eq(d1.at[i], d2.at[i]);
  }
}
END_TEST

START_TEST(test_matrix_SE3f_copy) {
  ksl_SE3f_t d1 =
    ksl_SE3f(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);
  ksl_SE3f_t d2;
  ksl_SE3f_copy(&d1, &d2);
  for(int i = 0; i < 12; i++) {
    ck_assert_float_eq(d1.at[i], d2.at[i]);
  }
}
END_TEST

START_TEST(test_matrix_SE3_inverted) {
  ksl_SE3_t d1 =
    ksl_SE3(0.28953232855036204, 0.13809317325638165, 0.9471543201739556, 1.0,
            -0.9568081313741513, 0.014614837677269121, 0.29035255510494773, 2.0,
            0.026253199052874154, -0.99031140658868, 0.13636013904305067, 3.0);
  ksl_SE3_t d2;
  ksl_SE3_inverted(&d1, &d2);
  ksl_SE3_t d3;
  ksl_product_dd(&d1, &d2, &d3);

  /* verify that the result is identity */
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert_double_eq_tol(1.0, d3.as_array[j][i], 1e-9);
      } else {
        ck_assert_double_eq_tol(0.0, d3.as_array[j][i], 1e-9);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3f_inverted) {
  ksl_SE3f_t d1 = ksl_SE3f(
    0.28953232855036204, 0.13809317325638165, 0.9471543201739556, 1.0,
    -0.9568081313741513, 0.014614837677269121, 0.29035255510494773, 2.0,
    0.026253199052874154, -0.99031140658868, 0.13636013904305067, 3.0);
  ksl_SE3f_t d2;
  ksl_SE3f_inverted(&d1, &d2);
  ksl_SE3f_t d3;
  ksl_product_ddf(&d1, &d2, &d3);

  /* verify that the result is identity */
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert_float_eq_tol(1.0, d3.as_array[j][i], 1e-6);
      } else {
        ck_assert_float_eq_tol(0.0, d3.as_array[j][i], 1e-6);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3_invert) {
  ksl_SE3_t d1 =
    ksl_SE3(0.28953232855036204, 0.13809317325638165, 0.9471543201739556, 1.0,
            -0.9568081313741513, 0.014614837677269121, 0.29035255510494773, 2.0,
            0.026253199052874154, -0.99031140658868, 0.13636013904305067, 3.0);
  ksl_SE3_t d2;
  ksl_SE3_copy(&d1, &d2);
  ksl_SE3_invert(&d1);
  ksl_SE3_t d3;
  ksl_product_dd(&d1, &d2, &d3);

  /* verify that the result is identity */
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert_double_eq_tol(1.0, d3.as_array[j][i], 1e-9);
      } else {
        ck_assert_double_eq_tol(0.0, d3.as_array[j][i], 1e-9);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_SE3f_invert) {
  ksl_SE3f_t d1 = ksl_SE3f(
    0.28953232855036204, 0.13809317325638165, 0.9471543201739556, 1.0,
    -0.9568081313741513, 0.014614837677269121, 0.29035255510494773, 2.0,
    0.026253199052874154, -0.99031140658868, 0.13636013904305067, 3.0);
  ksl_SE3f_t d2;
  ksl_SE3f_copy(&d1, &d2);
  ksl_SE3f_invert(&d1);
  ksl_SE3f_t d3;
  ksl_product_ddf(&d1, &d2, &d3);

  /* verify that the result is identity */
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 4; j++) {
      if(i == j) {
        ck_assert(fabs(1.0 - d3.as_array[j][i]) < 1e-6);
      } else {
        ck_assert(fabs(d3.as_array[j][i]) < 1e-6);
      }
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_getEulerAngles) {
  ksl_mat3x3_t r1;
  ksl_mat3x3_t r2;

  for(int j = 0; j < 12; j++) {
    ksl_mat3x3_setIdentity(&r1);
    ksl_axis_enum_t axis = j;
    ksl_vec3i_t sequence = ksl_axis_getVector(axis);
    ksl_vec3_t angles_test = {{0.2, 0.3, 0.4}};
    double dc[2];
    for(int i = 0; i < 3; i++) {
      ksl_dc(angles_test.at[i], dc);
      switch(sequence.at[i]) {
        case 0: {
          ksl_product_drdrx(&r1, dc, &r2);
          break;
        }
        case 1: {
          ksl_product_drdry(&r1, dc, &r2);
          break;
        }
        case 2: {
          ksl_product_drdrz(&r1, dc, &r2);
          break;
        }
      }
      r1 = r2;
    }

    ksl_vec3_t angles;
    ksl_vec3_t angles_prev = {{0.21, 0.31, 0.41}};
    ksl_mat3x3_getEulerAngles(&r1, axis, &angles, &angles_prev);

    for(int i = 0; i < 3; i++) {
      ck_assert_msg(fabs(angles_test.at[i] - angles.at[i]) < 1e-9,
                    "failed for angle sequence: %d%d%d\n", sequence.x,
                    sequence.y, sequence.z);
    }

    ksl_mat3x3_getEulerAngles(&r1, axis, &angles);
    for(int i = 0; i < 3; i++) {
      ck_assert_msg(fabs(angles_test.at[i] - angles.at[i]) < 1e-9,
                    "failed for angle sequence with reference: %d%d%d\n",
                    sequence.x, sequence.y, sequence.z);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat3x3_setFromEulerAngles) {
  ksl_mat3x3_t r1;
  ksl_mat3x3_t r2;

  for(int j = 0; j < 12; j++) {
    ksl_mat3x3_setIdentity(&r1);
    ksl_axis_enum_t axis = j;
    ksl_vec3i_t sequence = ksl_axis_getVector(axis);
    ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
    double dc[2];
    for(int i = 0; i < 3; i++) {
      ksl_dc(angles.at[i], dc);
      switch(sequence.at[i]) {
        case 0: {
          ksl_product_drdrx(&r1, dc, &r2);
          break;
        }
        case 1: {
          ksl_product_drdry(&r1, dc, &r2);
          break;
        }
        case 2: {
          ksl_product_drdrz(&r1, dc, &r2);
          break;
        }
      }
      r1 = r2;
    }
    ksl_mat3x3_setFromEulerAngles(&r2, axis, &angles);
    for(int i = 0; i < 9; i++) {
      ck_assert_msg(fabs(r1.at[i] - r2.at[i]) < 1e-6,
                    "setFromEulerAngles failed for angle sequence: %d%d%d",
                    sequence.at[0], sequence.at[1], sequence.at[2]);
    }
  }
}
END_TEST

START_TEST(test_matrix_mat4x4_getTranslation) {
  ksl_mat4x4_t m = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
                     12.0, 13.0, 14.0, 15.0, 16.0}};
  ksl_vec3_t v;
  ksl_mat4x4_getTranslation(&m, &v);
  ck_assert_double_eq(v.x, m.m03);
  ck_assert_double_eq(v.y, m.m13);
  ck_assert_double_eq(v.z, m.m23);
}
END_TEST

START_TEST(test_matrix_mat4x4f_getTranslation) {
  ksl_mat4x4f_t m = {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
                      12.0, 13.0, 14.0, 15.0, 16.0}};
  ksl_vec3f_t v;
  ksl_mat4x4f_getTranslation(&m, &v);
  ck_assert_float_eq(v.x, m.m03);
  ck_assert_float_eq(v.y, m.m13);
  ck_assert_float_eq(v.z, m.m23);
}
END_TEST

START_TEST(test_matrix_mat3x3_determinant) {
  ksl_mat3x3_t m = {{3.0, 2.0, 1.0, 4.0, 5.0, 6.0, 7.0, 9.0, 8.0}};
  double det = ksl_mat3x3_determinant(&m);
  ck_assert_double_eq_tol(det, -21, 1e-9);
}
END_TEST

START_TEST(test_matrix_mat3x3f_determinant) {
  ksl_mat3x3f_t m = {{3.0, 2.0, 1.0, 4.0, 5.0, 6.0, 7.0, 9.0, 8.0}};
  float det = ksl_mat3x3f_determinant(&m);
  ck_assert_float_eq_tol(det, -21, 1e-6);
}
END_TEST

START_TEST(test_matrix_mat3x3_getAxisAngle) {
  ksl_vec3_t axis = {{1.0, 2.0, 3.0}};
  double angle = 0.5;
  ksl_mat3x3_t m;
  ksl_mat3x3_setFromAxisAngle(&m, &axis, angle);

  ksl_vec3_t axis_test;
  double angle_test;
  ksl_mat3x3_getAxisAngle(&m, &axis_test, &angle_test);
  ksl_vec3_normalize(&axis);

  ck_assert_double_eq_tol(axis.x, axis_test.x, 1e-9);
  ck_assert_double_eq_tol(axis.y, axis_test.y, 1e-9);
  ck_assert_double_eq_tol(axis.z, axis_test.z, 1e-9);
  ck_assert_double_eq_tol(angle, angle_test, 1e-9);
}
END_TEST

START_TEST(test_matrix_mat3x3f_getAxisAngle) {
  ksl_vec3f_t axis = {{1.0, 2.0, 3.0}};
  float angle = 0.5;
  ksl_mat3x3f_t m;
  ksl_mat3x3f_setFromAxisAngle(&m, &axis, angle);

  ksl_vec3f_t axis_test;
  float angle_test;
  ksl_mat3x3f_getAxisAngle(&m, &axis_test, &angle_test);
  ksl_vec3f_normalize(&axis);

  ck_assert_float_eq_tol(axis.x, axis_test.x, 1e-6);
  ck_assert_float_eq_tol(axis.y, axis_test.y, 1e-6);
  ck_assert_float_eq_tol(axis.z, axis_test.z, 1e-6);
  ck_assert_float_eq_tol(angle, angle_test, 1e-6);
}
END_TEST

/* matrix vector operations */
START_TEST(test_matrix_product_drv) {

  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);

  FILE* f = fopen("product_drv_test.txt", "w");
  ksl_mat3x3_print(f, &r);
  fclose(f);
  ksl_vec3_t v1 = {{1.0, 2.0, 3.0}};
  ksl_vec3_t v2;
  ksl_product_drv(&r, &v1, &v2);

  ck_assert_double_eq_tol(v2.x, 1.0224326923807563, 1e-9);
  ck_assert_double_eq_tol(v2.y, 1.6260200151342845, 1e-9);
  ck_assert_double_eq_tol(v2.z, 3.2110263623853568, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drvf) {

  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v1 = {{1.0, 2.0, 3.0}};

  ksl_vec3f_t v2;
  ksl_product_drvf(&r, &v1, &v2);
  ck_assert_float_eq_tol(v2.x, 1.0224326923807563, 1e-6);
  ck_assert_float_eq_tol(v2.y, 1.6260200151342845, 1e-6);
  ck_assert_float_eq_tol(v2.z, 3.2110263623853568, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drvinv) {

  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v1 = {{1.0, 2.0, 3.0}};
  ksl_vec3_t v2;
  ksl_product_drvinv(&r, &v1, &v2);
  ck_assert_double_eq_tol(v2.x, -1.0224326923807563, 1e-9);
  ck_assert_double_eq_tol(v2.y, -1.6260200151342845, 1e-9);
  ck_assert_double_eq_tol(v2.z, -3.2110263623853568, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drvinvf) {

  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v1 = {{1.0, 2.0, 3.0}};

  ksl_vec3f_t v2;
  ksl_product_drvinvf(&r, &v1, &v2);
  ck_assert_float_eq_tol(v2.x, -1.0224326923807563, 1e-6);
  ck_assert_float_eq_tol(v2.y, -1.6260200151342845, 1e-6);
  ck_assert_float_eq_tol(v2.z, -3.2110263623853568, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drinvv) {

  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v1 = {{1.0, 2.0, 3.0}};
  ksl_vec3_t v2;
  ksl_product_drinvv(&r, &v1, &v2);
  ck_assert_double_eq_tol(v2.x, 1.1831846399394617, 1e-9);
  ck_assert_double_eq_tol(v2.y, 2.274971321748124, 1e-9);
  ck_assert_double_eq_tol(v2.z, 2.7248081754565625, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drinvvf) {

  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v1 = {{1.0, 2.0, 3.0}};

  ksl_vec3f_t v2;
  ksl_product_drinvvf(&r, &v1, &v2);
  ck_assert_float_eq_tol(v2.x, 1.1831846399394617, 1e-6);
  ck_assert_float_eq_tol(v2.y, 2.274971321748124, 1e-6);
  ck_assert_float_eq_tol(v2.z, 2.7248081754565625, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drinvvinv) {

  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v1 = {{1.0, 2.0, 3.0}};
  ksl_vec3_t v2;
  ksl_product_drinvvinv(&r, &v1, &v2);
  ck_assert_double_eq_tol(v2.x, -1.1831846399394617, 1e-9);
  ck_assert_double_eq_tol(v2.y, -2.274971321748124, 1e-9);
  ck_assert_double_eq_tol(v2.z, -2.7248081754565625, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drinvvinvf) {

  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v1 = {{1.0, 2.0, 3.0}};

  ksl_vec3f_t v2;
  ksl_product_drinvvinvf(&r, &v1, &v2);
  ck_assert_float_eq_tol(v2.x, -1.1831846399394617, 1e-6);
  ck_assert_float_eq_tol(v2.y, -2.274971321748124, 1e-6);
  ck_assert_float_eq_tol(v2.z, -2.7248081754565625, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drvtx) {
  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v2;
  ksl_product_drvtx(&r, 2.0, &v2);
  ck_assert_double_eq_tol(v2.x, 1.759846352562514, 1e-9);
  ck_assert_double_eq_tol(v2.y, 0.8714642629237408, 1e-9);
  ck_assert_double_eq_tol(v2.z, -0.3788018661770242, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drvtxf) {
  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v2;
  ksl_product_drvtxf(&r, 2.0, &v2);
  ck_assert_float_eq_tol(v2.x, 1.759846352562514, 1e-6);
  ck_assert_float_eq_tol(v2.y, 0.8714642629237408, 1e-6);
  ck_assert_float_eq_tol(v2.z, -0.3788018661770242, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drvtxinv) {
  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v2;
  ksl_product_drvtxinv(&r, 2.0, &v2);
  ck_assert_double_eq_tol(v2.x, -1.759846352562514, 1e-9);
  ck_assert_double_eq_tol(v2.y, -0.8714642629237408, 1e-9);
  ck_assert_double_eq_tol(v2.z, 0.3788018661770242, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drvtxinvf) {
  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v2;
  ksl_product_drvtxinvf(&r, 2.0, &v2);
  ck_assert_float_eq_tol(v2.x, -1.759846352562514, 1e-6);
  ck_assert_float_eq_tol(v2.y, -0.8714642629237408, 1e-6);
  ck_assert_float_eq_tol(v2.z, 0.3788018661770242, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drvty) {
  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v2;
  ksl_product_drvty(&r, 2.0, &v2);
  ck_assert_double_eq_tol(v2.x, -0.7440511038845192, 1e-9);
  ck_assert_double_eq_tol(v2.y, 1.7596760666084763, 1e-9);
  ck_assert_double_eq_tol(v2.z, 0.5915472047212714, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drvtyf) {
  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v2;
  ksl_product_drvtyf(&r, 2.0, &v2);
  ck_assert_float_eq_tol(v2.x, -0.7440511038845192, 1e-6);
  ck_assert_float_eq_tol(v2.y, 1.7596760666084763, 1e-6);
  ck_assert_float_eq_tol(v2.z, 0.5915472047212714, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drvtyinv) {
  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v2;
  ksl_product_drvtyinv(&r, 2.0, &v2);
  ck_assert_double_eq_tol(v2.x, 0.7440511038845192, 1e-9);
  ck_assert_double_eq_tol(v2.y, -1.7596760666084763, 1e-9);
  ck_assert_double_eq_tol(v2.z, -0.5915472047212714, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drvtyinvf) {
  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v2;
  ksl_product_drvtyinvf(&r, 2.0, &v2);
  ck_assert_float_eq_tol(v2.x, 0.7440511038845192, 1e-6);
  ck_assert_float_eq_tol(v2.y, -1.7596760666084763, 1e-6);
  ck_assert_float_eq_tol(v2.z, -0.5915472047212714, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drvtz) {
  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v2;
  ksl_product_drvtz(&r, 2.0, &v2);
  ck_assert_double_eq_tol(v2.x, 0.591040413322679, 1e-9);
  ck_assert_double_eq_tol(v2.y, -0.3795921219573748, 1e-9);
  ck_assert_double_eq_tol(v2.z, 1.8725867271683985, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drvtzf) {
  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v2;
  ksl_product_drvtzf(&r, 2.0, &v2);
  ck_assert_float_eq_tol(v2.x, 0.591040413322679, 1e-6);
  ck_assert_float_eq_tol(v2.y, -0.3795921219573748, 1e-6);
  ck_assert_float_eq_tol(v2.z, 1.8725867271683985, 1e-6);
}
END_TEST

START_TEST(test_matrix_product_drvtzinv) {
  ksl_mat3x3_t r;
  ksl_mat3x3_setIdentity(&r);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3_t v2;
  ksl_product_drvtzinv(&r, 2.0, &v2);
  ck_assert_double_eq_tol(v2.x, -0.591040413322679, 1e-9);
  ck_assert_double_eq_tol(v2.y, 0.3795921219573748, 1e-9);
  ck_assert_double_eq_tol(v2.z, -1.8725867271683985, 1e-9);
}
END_TEST

START_TEST(test_matrix_product_drvtzinvf) {
  ksl_mat3x3f_t r;
  ksl_mat3x3f_setIdentity(&r);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r, KSL_AXIS_XYZ, &angles);
  ksl_vec3f_t v2;
  ksl_product_drvtzinvf(&r, 2.0, &v2);
  ck_assert_float_eq_tol(v2.x, -0.591040413322679, 1e-6);
  ck_assert_float_eq_tol(v2.y, 0.3795921219573748, 1e-6);
  ck_assert_float_eq_tol(v2.z, -1.8725867271683985, 1e-6);
}
END_TEST

/* matrix-matrix operations */
START_TEST(test_matrix_product_drdrx) {
  ksl_mat3x3_t r1;
  ksl_mat3x3_setIdentity(&r1);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3_t r2;
  double dc[2];
  ksl_dc(0.5, dc);
  ksl_product_drdrx(&r1, dc, &r2);

  ksl_mat3x3_t r3 =
    ksl_mat3x3(0.879923176281257, -0.1848032027151301, 0.4377019306666745,
               0.4357321314618704, 0.6811374365560571, -0.588378536431725,
               -0.1894009330885121, 0.7084487058270867, 0.6798733100785226);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq_tol(r2.at[i], r3.at[i], 1e-9);
  }
}
END_TEST

START_TEST(test_matrix_product_drdrxf) {
  ksl_mat3x3f_t r1;
  ksl_mat3x3f_setIdentity(&r1);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3f_t r2;
  float dc[2];
  ksl_dcf(0.5, dc);
  ksl_product_drdrxf(&r1, dc, &r2);

  ksl_mat3x3f_t r3 =
    ksl_mat3x3f(0.879923176281257, -0.1848032027151301, 0.4377019306666745,
                0.4357321314618704, 0.6811374365560571, -0.588378536431725,
                -0.1894009330885121, 0.7084487058270867, 0.6798733100785226);

  for(int i = 0; i < 9; i++) {
    ck_assert_float_eq_tol(r2.at[i], r3.at[i], 1e-6);
  }
}
END_TEST

START_TEST(test_matrix_product_drdrxinv) {
  ksl_mat3x3_t r1;
  ksl_mat3x3_setIdentity(&r1);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3_t r2;
  double dc[2];
  ksl_dc(0.5, dc);
  ksl_product_drdrxinv(&r1, dc, &r2);

  ksl_mat3x3_t r3 =
    ksl_mat3x3(0.879923176281257, -0.4681630712092062, 0.080984829437787,
               0.4357321314618704, 0.8631235940753837, 0.2552551095709692,
               -0.1894009330885121, -0.1893171944287045, 0.963476147311829);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq_tol(r2.at[i], r3.at[i], 1e-9);
  }
}
END_TEST

START_TEST(test_matrix_product_drdrxinvf) {
  ksl_mat3x3f_t r1;
  ksl_mat3x3f_setIdentity(&r1);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3f_t r2;
  float dc[2];
  ksl_dcf(0.5, dc);
  ksl_product_drdrxinvf(&r1, dc, &r2);

  ksl_mat3x3f_t r3 =
    ksl_mat3x3f(0.879923176281257, -0.4681630712092062, 0.080984829437787,
                0.4357321314618704, 0.8631235940753837, 0.2552551095709692,
                -0.1894009330885121, -0.1893171944287045, 0.963476147311829);

  for(int i = 0; i < 9; i++) {
    ck_assert_float_eq_tol(r2.at[i], r3.at[i], 1e-6);
  }
}
END_TEST

START_TEST(test_matrix_product_drdry) {
  ksl_mat3x3_t r1;
  ksl_mat3x3_setIdentity(&r1);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3_t r2;
  double dc[2];
  ksl_dc(0.5, dc);
  ksl_product_drdry(&r1, dc, &r2);

  ksl_mat3x3_t r3 =
    ksl_mat3x3(0.6305253010605816, -0.3720255519422596, 0.6812010227711935,
               0.4733839989859243, 0.8798380333042382, 0.0423393983828867,
               -0.6150979062121391, 0.2957736023606357, 0.7308710843370773);

  FILE* f = fopen("drdry.txt", "w");

  ksl_mat3x3_print(f, &r1, "r1: ");
  ksl_mat3x3_print(f, &r2, "r2: ");
  fclose(f);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq_tol(r2.at[i], r3.at[i], 1e-9);
  }
}
END_TEST

START_TEST(test_matrix_product_drdryf) {
  ksl_mat3x3f_t r1;
  ksl_mat3x3f_setIdentity(&r1);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3f_t r2;
  float dc[2];
  ksl_dcf(0.5, dc);
  ksl_product_drdryf(&r1, dc, &r2);

  ksl_mat3x3f_t r3 =
    ksl_mat3x3f(0.6305253010605816, -0.3720255519422596, 0.6812010227711935,
                0.4733839989859243, 0.8798380333042382, 0.0423393983828867,
                -0.6150979062121391, 0.2957736023606357, 0.7308710843370773);

  for(int i = 0; i < 9; i++) {
    ck_assert_float_eq_tol(r2.at[i], r3.at[i], 1e-6);
  }
}
END_TEST

START_TEST(test_matrix_product_drdryinv) {
  ksl_mat3x3_t r1;
  ksl_mat3x3_setIdentity(&r1);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3_t r2;
  double dc[2];
  ksl_dc(0.5, dc);
  ksl_product_drdryinv(&r1, dc, &r2);

  ksl_mat3x3_t r3 =
    ksl_mat3x3(0.9138851695546577, -0.3720255519422596, -0.1625142626667319,
               0.2913978414665975, 0.8798380333042382, -0.3754628252436425,
               0.2826679940436521, 0.2957736023606357, 0.9124783730532743);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq_tol(r2.at[i], r3.at[i], 1e-9);
  }
}
END_TEST

START_TEST(test_matrix_product_drdryinvf) {
  ksl_mat3x3f_t r1;
  ksl_mat3x3f_setIdentity(&r1);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3f_t r2;
  float dc[2];
  ksl_dcf(0.5, dc);
  ksl_product_drdryinvf(&r1, dc, &r2);

  ksl_mat3x3f_t r3 =
    ksl_mat3x3f(0.9138851695546577, -0.3720255519422596, -0.1625142626667319,
                0.2913978414665975, 0.8798380333042382, -0.3754628252436425,
                0.2826679940436521, 0.2957736023606357, 0.9124783730532743);

  for(int i = 0; i < 9; i++) {
    ck_assert_float_eq_tol(r2.at[i], r3.at[i], 1e-6);
  }
}
END_TEST

START_TEST(test_matrix_product_drdrz) {
  ksl_mat3x3_t r1;
  ksl_mat3x3_setIdentity(&r1);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3_t r2;
  double dc[2];
  ksl_dc(0.5, dc);
  ksl_product_drdrz(&r1, dc, &r2);

  ksl_mat3x3_t r3 =
    ksl_mat3x3(0.5938466846931759, -0.7483407796811308, 0.2955202066613395,
               0.804207743227608, 0.5632294035024559, -0.1897960609786874,
               -0.0244135374675904, 0.3503694000572896, 0.9362933635841992);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq_tol(r2.at[i], r3.at[i], 1e-9);
  }
}
END_TEST

START_TEST(test_matrix_product_drdrzf) {
  ksl_mat3x3f_t r1;
  ksl_mat3x3f_setIdentity(&r1);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3f_t r2;
  float dc[2];
  ksl_dcf(0.5, dc);
  ksl_product_drdrzf(&r1, dc, &r2);

  ksl_mat3x3f_t r3 =
    ksl_mat3x3f(0.5938466846931759, -0.7483407796811308, 0.2955202066613395,
                0.804207743227608, 0.5632294035024559, -0.1897960609786874,
                -0.0244135374675904, 0.3503694000572896, 0.9362933635841992);

  for(int i = 0; i < 9; i++) {
    ck_assert_float_eq_tol(r2.at[i], r3.at[i], 1e-6);
  }
}
END_TEST

START_TEST(test_matrix_product_drdrzinv) {
  ksl_mat3x3_t r1;
  ksl_mat3x3_setIdentity(&r1);
  ksl_vec3_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3_t r2;
  double dc[2];
  ksl_dc(0.5, dc);
  ksl_product_drdrzinv(&r1, dc, &r2);

  ksl_mat3x3_t r3 =
    ksl_mat3x3(0.9505637859220634, 0.0953745057567946, 0.2955202066613395,
               -0.0394259027750862, 0.981031627128985, -0.1897960609786874,
               -0.3080163747008967, 0.1687621113410926, 0.9362933635841992);

  for(int i = 0; i < 9; i++) {
    ck_assert_double_eq_tol(r2.at[i], r3.at[i], 1e-9);
  }
}
END_TEST

START_TEST(test_matrix_product_drdrzinvf) {
  ksl_mat3x3f_t r1;
  ksl_mat3x3f_setIdentity(&r1);
  ksl_vec3f_t angles = {{0.2, 0.3, 0.4}};
  ksl_mat3x3f_setFromEulerAngles(&r1, KSL_AXIS_XYZ, &angles);
  ksl_mat3x3f_t r2;
  float dc[2];
  ksl_dcf(0.5, dc);
  ksl_product_drdrzinvf(&r1, dc, &r2);

  ksl_mat3x3f_t r3 =
    ksl_mat3x3f(0.9505637859220634, 0.0953745057567946, 0.2955202066613395,
                -0.0394259027750862, 0.981031627128985, -0.1897960609786874,
                -0.3080163747008967, 0.1687621113410926, 0.9362933635841992);

  for(int i = 0; i < 9; i++) {
    ck_assert_float_eq_tol(r2.at[i], r3.at[i], 1e-6);
  }
}
END_TEST

//
// inline void ksl_product_drdr(const ksl_mat3x3_t* restrict r1i,
//                              const ksl_mat3x3_t* restrict r2i,
//                              ksl_mat3x3_t* restrict ro) {
//   ro->m00 = r1i->m00 * r2i->m00 + r1i->m01 * r2i->m10 + r1i->m02 * r2i->m20;
//   ro->m01 = r1i->m00 * r2i->m01 + r1i->m01 * r2i->m11 + r1i->m02 * r2i->m21;
//   ro->m02 = r1i->m00 * r2i->m02 + r1i->m01 * r2i->m12 + r1i->m02 * r2i->m22;
//   ro->m10 = r1i->m10 * r2i->m00 + r1i->m11 * r2i->m10 + r1i->m12 * r2i->m20;
//   ro->m11 = r1i->m10 * r2i->m01 + r1i->m11 * r2i->m11 + r1i->m12 * r2i->m21;
//   ro->m12 = r1i->m10 * r2i->m02 + r1i->m11 * r2i->m12 + r1i->m12 * r2i->m22;
//   ro->m20 = r1i->m20 * r2i->m00 + r1i->m21 * r2i->m10 + r1i->m22 * r2i->m20;
//   ro->m21 = r1i->m20 * r2i->m01 + r1i->m21 * r2i->m11 + r1i->m22 * r2i->m21;
//   ro->m22 = r1i->m20 * r2i->m02 + r1i->m21 * r2i->m12 + r1i->m22 * r2i->m22;
// }
//
// inline void ksl_product_drdrf(const ksl_mat3x3f_t* restrict r1i,
//                               const ksl_mat3x3f_t* restrict r2i,
//                               ksl_mat3x3f_t* restrict ro) {
//   ro->m00 = r1i->m00 * r2i->m00 + r1i->m01 * r2i->m10 + r1i->m02 * r2i->m20;
//   ro->m01 = r1i->m00 * r2i->m01 + r1i->m01 * r2i->m11 + r1i->m02 * r2i->m21;
//   ro->m02 = r1i->m00 * r2i->m02 + r1i->m01 * r2i->m12 + r1i->m02 * r2i->m22;
//   ro->m10 = r1i->m10 * r2i->m00 + r1i->m11 * r2i->m10 + r1i->m12 * r2i->m20;
//   ro->m11 = r1i->m10 * r2i->m01 + r1i->m11 * r2i->m11 + r1i->m12 * r2i->m21;
//   ro->m12 = r1i->m10 * r2i->m02 + r1i->m11 * r2i->m12 + r1i->m12 * r2i->m22;
//   ro->m20 = r1i->m20 * r2i->m00 + r1i->m21 * r2i->m10 + r1i->m22 * r2i->m20;
//   ro->m21 = r1i->m20 * r2i->m01 + r1i->m21 * r2i->m11 + r1i->m22 * r2i->m21;
//   ro->m22 = r1i->m20 * r2i->m02 + r1i->m21 * r2i->m12 + r1i->m22 * r2i->m22;
// }
//
// inline void ksl_product_drdrinv(const ksl_mat3x3_t* restrict r1i,
//                                 const ksl_mat3x3_t* restrict r2i,
//                                 ksl_mat3x3_t* restrict ro) {
//   ro->m00 = ksl_dot_vv(&r1i->v0, &r2i->v0);
//   ro->m01 = ksl_dot_vv(&r1i->v0, &r2i->v1);
//   ro->m02 = ksl_dot_vv(&r1i->v0, &r2i->v2);
//   ro->m10 = ksl_dot_vv(&r1i->v1, &r2i->v0);
//   ro->m11 = ksl_dot_vv(&r1i->v1, &r2i->v1);
//   ro->m12 = ksl_dot_vv(&r1i->v1, &r2i->v2);
//   ro->m20 = ksl_dot_vv(&r1i->v2, &r2i->v0);
//   ro->m21 = ksl_dot_vv(&r1i->v2, &r2i->v1);
//   ro->m22 = ksl_dot_vv(&r1i->v2, &r2i->v2);
// }
//
// inline void ksl_product_drdrinvf(const ksl_mat3x3f_t* restrict r1i,
//                                  const ksl_mat3x3f_t* restrict r2i,
//                                  ksl_mat3x3f_t* restrict ro) {
//   ro->m00 = ksl_dot_vvf(&r1i->v0, &r2i->v0);
//   ro->m01 = ksl_dot_vvf(&r1i->v0, &r2i->v1);
//   ro->m02 = ksl_dot_vvf(&r1i->v0, &r2i->v2);
//   ro->m10 = ksl_dot_vvf(&r1i->v1, &r2i->v0);
//   ro->m11 = ksl_dot_vvf(&r1i->v1, &r2i->v1);
//   ro->m12 = ksl_dot_vvf(&r1i->v1, &r2i->v2);
//   ro->m20 = ksl_dot_vvf(&r1i->v2, &r2i->v0);
//   ro->m21 = ksl_dot_vvf(&r1i->v2, &r2i->v1);
//   ro->m22 = ksl_dot_vvf(&r1i->v2, &r2i->v2);
// }
//
// inline void ksl_product_drinvdr(const ksl_mat3x3_t* restrict r1i,
//                                 const ksl_mat3x3_t* restrict r2i,
//                                 ksl_mat3x3_t* restrict ro) {
//   ro->m00 = r1i->m00 * r2i->m00 + r1i->m10 * r2i->m10 + r1i->m20 * r2i->m20;
//   ro->m01 = r1i->m00 * r2i->m01 + r1i->m10 * r2i->m11 + r1i->m20 * r2i->m21;
//   ro->m02 = r1i->m00 * r2i->m02 + r1i->m10 * r2i->m12 + r1i->m20 * r2i->m22;
//   ro->m10 = r1i->m01 * r2i->m00 + r1i->m11 * r2i->m10 + r1i->m21 * r2i->m20;
//   ro->m11 = r1i->m01 * r2i->m01 + r1i->m11 * r2i->m11 + r1i->m21 * r2i->m21;
//   ro->m12 = r1i->m01 * r2i->m02 + r1i->m11 * r2i->m12 + r1i->m21 * r2i->m22;
//   ro->m20 = r1i->m02 * r2i->m00 + r1i->m12 * r2i->m10 + r1i->m22 * r2i->m20;
//   ro->m21 = r1i->m02 * r2i->m01 + r1i->m12 * r2i->m11 + r1i->m22 * r2i->m21;
//   ro->m22 = r1i->m02 * r2i->m02 + r1i->m12 * r2i->m12 + r1i->m22 * r2i->m22;
// }
//
// inline void ksl_product_drinvdrf(const ksl_mat3x3f_t* restrict r1i,
//                                  const ksl_mat3x3f_t* restrict r2i,
//                                  ksl_mat3x3f_t* restrict ro) {
//   ro->m00 = r1i->m00 * r2i->m00 + r1i->m10 * r2i->m10 + r1i->m20 * r2i->m20;
//   ro->m01 = r1i->m00 * r2i->m01 + r1i->m10 * r2i->m11 + r1i->m20 * r2i->m21;
//   ro->m02 = r1i->m00 * r2i->m02 + r1i->m10 * r2i->m12 + r1i->m20 * r2i->m22;
//   ro->m10 = r1i->m01 * r2i->m00 + r1i->m11 * r2i->m10 + r1i->m21 * r2i->m20;
//   ro->m11 = r1i->m01 * r2i->m01 + r1i->m11 * r2i->m11 + r1i->m21 * r2i->m21;
//   ro->m12 = r1i->m01 * r2i->m02 + r1i->m11 * r2i->m12 + r1i->m21 * r2i->m22;
//   ro->m20 = r1i->m02 * r2i->m00 + r1i->m12 * r2i->m10 + r1i->m22 * r2i->m20;
//   ro->m21 = r1i->m02 * r2i->m01 + r1i->m12 * r2i->m11 + r1i->m22 * r2i->m21;
//   ro->m22 = r1i->m02 * r2i->m02 + r1i->m12 * r2i->m12 + r1i->m22 * r2i->m22;
// }
//
// inline void ksl_product_dv(const ksl_SE3_t* restrict Di,
//                            const ksl_vec3_t* restrict vi,
//                            ksl_vec3_t* restrict vo) {
//   ksl_product_drv(&(Di->R), vi, vo);
//   ksl_xpy_vv(&Di->t, vo);
// }
//
// inline void ksl_product_dinvv(const ksl_SE3_t* restrict Di,
//                               const ksl_vec3_t* restrict vi,
//                               ksl_vec3_t* restrict vo) {
//   ksl_vec3_t temp;
//   ksl_product_drinvv(&Di->R, &Di->t, &temp); /* R^-1 * t -> temp*/
//   ksl_product_drinvv(&Di->R, vi, vo);        /* R^-1 * vi -> vo */
//   ksl_nxpy_vv(&temp, vo);                    /* vo -= temp */
// }
//
// inline void ksl_product_ddrx(const ksl_SE3_t* restrict Di, const double
// dc[2],
//                              ksl_SE3_t* restrict Do) {
//   ksl_product_drdrx(&Di->R, dc, &Do->R);
//   ksl_vec3_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrxf(const ksl_SE3f_t* restrict Di, const float
// dc[2],
//                               ksl_SE3f_t* restrict Do) {
//   ksl_product_drdrxf(&Di->R, dc, &Do->R);
//   ksl_vec3f_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrxinv(const ksl_SE3_t* restrict Di,
//                                 const double dc[2], ksl_SE3_t* restrict Do) {
//   ksl_product_drdrxinv(&Di->R, dc, &Do->R);
//   ksl_vec3_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrxinvf(const ksl_SE3f_t* restrict Di,
//                                  const float dc[2], ksl_SE3f_t* restrict Do)
//                                  {
//   ksl_product_drdrxinvf(&Di->R, dc, &Do->R);
//   ksl_vec3f_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddry(const ksl_SE3_t* restrict Di, const double
// dc[2],
//                              ksl_SE3_t* restrict Do) {
//   ksl_product_drdry(&Di->R, dc, &Do->R);
//   ksl_vec3_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddryf(const ksl_SE3f_t* restrict Di, const float
// dc[2],
//                               ksl_SE3f_t* restrict Do) {
//   ksl_product_drdryf(&Di->R, dc, &Do->R);
//   ksl_vec3f_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddryinv(const ksl_SE3_t* restrict Di,
//                                 const double dc[2], ksl_SE3_t* restrict Do) {
//   ksl_product_drdryinv(&Di->R, dc, &Do->R);
//   ksl_vec3_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddryinvf(const ksl_SE3f_t* restrict Di,
//                                  const float dc[2], ksl_SE3f_t* restrict Do)
//                                  {
//   ksl_product_drdryinvf(&Di->R, dc, &Do->R);
//   ksl_vec3f_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrz(const ksl_SE3_t* restrict Di, const double
// dc[2],
//                              ksl_SE3_t* restrict Do) {
//   ksl_product_drdrz(&Di->R, dc, &Do->R);
//   ksl_vec3_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrzf(const ksl_SE3f_t* restrict Di, const float
// dc[2],
//                               ksl_SE3f_t* restrict Do) {
//   ksl_product_drdrzf(&Di->R, dc, &Do->R);
//   ksl_vec3f_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrzinv(const ksl_SE3_t* restrict Di,
//                                 const double dc[2], ksl_SE3_t* restrict Do) {
//   ksl_product_drdrzinv(&Di->R, dc, &Do->R);
//   ksl_vec3_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrzinvf(const ksl_SE3f_t* restrict Di,
//                                  const float dc[2], ksl_SE3f_t* restrict Do)
//                                  {
//   ksl_product_drdrzinvf(&Di->R, dc, &Do->R);
//   ksl_vec3f_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddt(const ksl_SE3_t* restrict Di,
//                             const ksl_vec3_t* restrict t,
//                             ksl_SE3_t* restrict Do) {
//   ksl_mat3x3_copy(&Di->R, &Do->R);
//   ksl_product_drv(&Di->R, t, &Do->t);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtf(const ksl_SE3f_t* restrict Di,
//                              const ksl_vec3f_t* restrict t,
//                              ksl_SE3f_t* restrict Do) {
//   ksl_mat3x3f_copy(&Di->R, &Do->R);
//   ksl_product_drvf(&Di->R, t, &Do->t);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtinv(const ksl_SE3_t* restrict Di,
//                                const ksl_vec3_t* restrict t,
//                                ksl_SE3_t* restrict Do) {
//   ksl_mat3x3_copy(&Di->R, &Do->R);
//   ksl_product_drvinv(&Di->R, t, &Do->t);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtinvf(const ksl_SE3f_t* restrict Di,
//                                 const ksl_vec3f_t* restrict t,
//                                 ksl_SE3f_t* restrict Do) {
//   ksl_mat3x3f_copy(&Di->R, &Do->R);
//   ksl_product_drvinvf(&Di->R, t, &Do->t);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_dinvdt(const ksl_SE3_t* restrict Di,
//                                const ksl_vec3_t* restrict t,
//                                ksl_SE3_t* restrict Do) {
//   ksl_subtract_vv(t, &Di->t, &Do->t);
//   ksl_mat3x3_transposed(&Di->R, &Do->R);
//   ksl_product_drinvv(&Di->R, &Do->t, &Do->t);
// }
//
// inline void ksl_product_dinvdtf(const ksl_SE3f_t* restrict Di,
//                                 const ksl_vec3f_t* restrict t,
//                                 ksl_SE3f_t* restrict Do) {
//   ksl_subtract_vvf(t, &Di->t, &Do->t);
//   ksl_mat3x3f_transposed(&Di->R, &Do->R);
//   ksl_product_drinvvf(&Di->R, &Do->t, &Do->t);
// }
//
// inline void ksl_product_ddtx(const ksl_SE3_t* restrict Di, const double ti,
//                              ksl_SE3_t* restrict Do) {
//   ksl_mat3x3_copy(&Di->R, &Do->R);
//   ksl_product_drvtx(&Di->R, ti, &Do->t);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtxf(const ksl_SE3f_t* restrict Di, const float ti,
//                               ksl_SE3f_t* restrict Do) {
//   ksl_mat3x3f_copy(&Di->R, &Do->R);
//   ksl_product_drvtxf(&Di->R, ti, &Do->t);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtxinv(const ksl_SE3_t* restrict Di, const double
// ti,
//                                 ksl_SE3_t* restrict Do) {
//   ksl_mat3x3_copy(&Di->R, &Do->R);
//   ksl_product_drvtxinv(&Di->R, ti, &Do->t);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtxinvf(const ksl_SE3f_t* restrict Di, const float
// ti,
//                                  ksl_SE3f_t* restrict Do) {
//   ksl_mat3x3f_copy(&Di->R, &Do->R);
//   ksl_product_drvtxinvf(&Di->R, ti, &Do->t);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddty(const ksl_SE3_t* restrict Di, const double ti,
//                              ksl_SE3_t* restrict Do) {
//   ksl_mat3x3_copy(&Di->R, &Do->R);
//   ksl_product_drvty(&Di->R, ti, &Do->t);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtyf(const ksl_SE3f_t* restrict Di, const float ti,
//                               ksl_SE3f_t* restrict Do) {
//   ksl_mat3x3f_copy(&Di->R, &Do->R);
//   ksl_product_drvtyf(&Di->R, ti, &Do->t);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtyinv(const ksl_SE3_t* restrict Di, const double
// ti,
//                                 ksl_SE3_t* restrict Do) {
//   ksl_mat3x3_copy(&Di->R, &Do->R);
//   ksl_product_drvtyinv(&Di->R, ti, &Do->t);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtyinvf(const ksl_SE3f_t* restrict Di, const float
// ti,
//                                  ksl_SE3f_t* restrict Do) {
//   ksl_mat3x3f_copy(&Di->R, &Do->R);
//   ksl_product_drvtyinvf(&Di->R, ti, &Do->t);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtz(const ksl_SE3_t* restrict Di, const double ti,
//                              ksl_SE3_t* restrict Do) {
//   ksl_mat3x3_copy(&Di->R, &Do->R);
//   ksl_product_drvtz(&Di->R, ti, &Do->t);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtzf(const ksl_SE3f_t* restrict Di, const float ti,
//                               ksl_SE3f_t* restrict Do) {
//   ksl_mat3x3f_copy(&Di->R, &Do->R);
//   ksl_product_drvtzf(&Di->R, ti, &Do->t);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtzinv(const ksl_SE3_t* restrict Di, const double
// ti,
//                                 ksl_SE3_t* restrict Do) {
//   ksl_mat3x3_copy(&Di->R, &Do->R);
//   ksl_product_drvtzinv(&Di->R, ti, &Do->t);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddtzinvf(const ksl_SE3f_t* restrict Di, const float
// ti,
//                                  ksl_SE3f_t* restrict Do) {
//   ksl_mat3x3f_copy(&Di->R, &Do->R);
//   ksl_product_drvtzinvf(&Di->R, ti, &Do->t);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddr(const ksl_SE3_t* restrict Di,
//                             const ksl_mat3x3_t* restrict Ri,
//                             ksl_SE3_t* restrict Do) {
//   ksl_product_drdr(&Di->R, Ri, &Do->R);
//   ksl_vec3_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrf(const ksl_SE3f_t* restrict Di,
//                              const ksl_mat3x3f_t* restrict Ri,
//                              ksl_SE3f_t* restrict Do) {
//   ksl_product_drdrf(&Di->R, Ri, &Do->R);
//   ksl_vec3f_copy(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrinv(const ksl_SE3_t* restrict Di,
//                                const ksl_mat3x3_t* restrict Ri,
//                                ksl_SE3_t* restrict Do) {
//   ksl_product_drdrinv(&Di->R, Ri, &Do->R);
//   ksl_xpy_vv(&Di->t, &Do->t);
// }
//
// inline void ksl_product_ddrinvf(const ksl_SE3f_t* restrict Di,
//                                 const ksl_mat3x3f_t* restrict Ri,
//                                 ksl_SE3f_t* restrict Do) {
//   ksl_product_drdrinvf(&Di->R, Ri, &Do->R);
//   ksl_xpy_vvf(&Di->t, &Do->t);
// }
//
// inline void ksl_product_dd(const ksl_SE3_t* restrict D1i,
//                            const ksl_SE3_t* restrict D2i,
//                            ksl_SE3_t* restrict Do) {
//   ksl_product_drdr(&D1i->R, &D2i->R, &Do->R);
//   ksl_product_drv(&D1i->R, &D2i->t, &Do->t);
//   ksl_xpy_vv(&D1i->t, &Do->t);
// }
//
// inline void ksl_product_ddf(const ksl_SE3f_t* restrict D1i,
//                             const ksl_SE3f_t* restrict D2i,
//                             ksl_SE3f_t* restrict Do) {
//   ksl_product_drdrf(&D1i->R, &D2i->R, &Do->R);
//   ksl_product_drvf(&D1i->R, &D2i->t, &Do->t);
//   ksl_xpy_vvf(&D1i->t, &Do->t);
// }
//
// inline void ksl_product_ddinv(const ksl_SE3_t* restrict D1i,
//                               const ksl_SE3_t* restrict D2i,
//                               ksl_SE3_t* restrict Do) {
//   ksl_product_drdrinv(&D1i->R, &D2i->R, &Do->R);
//   ksl_product_drvinv(&Do->R, &D2i->t, &Do->t);
//   ksl_xpy_vv(&D1i->t, &Do->t);
// }
//
// inline void ksl_product_ddinvf(const ksl_SE3f_t* D1i, const ksl_SE3f_t* D2i,
//                                ksl_SE3f_t* Do) {
//   ksl_product_drdrinvf(&D1i->R, &D2i->R, &Do->R);
//   ksl_product_drvinvf(&Do->R, &D2i->t, &Do->t);
//   ksl_xpy_vvf(&D1i->t, &Do->t);
// }

Suite* matrix_suite(void) {
  Suite* s = suite_create("matrix");
  TCase* tc_core = tcase_create("core");
  tcase_add_test(tc_core, test_matrix_SE3);
  tcase_add_test(tc_core, test_matrix_SE3f);
  tcase_add_test(tc_core, test_matrix_mat3x3);
  tcase_add_test(tc_core, test_matrix_mat3x3f);
  tcase_add_test(tc_core, test_matrix_mat4x4);
  tcase_add_test(tc_core, test_matrix_mat4x4f);
  tcase_add_test(tc_core, test_matrix_SE3_alloc);
  tcase_add_test(tc_core, test_matrix_SE3f_alloc);
  tcase_add_test(tc_core, test_matrix_mat3x3_alloc);
  tcase_add_test(tc_core, test_matrix_mat3x3f_alloc);
  tcase_add_test(tc_core, test_matrix_mat4x4_alloc);
  tcase_add_test(tc_core, test_matrix_mat4x4f_alloc);
  tcase_add_test(tc_core, test_matrix_SE3toMat4x4);
  tcase_add_test(tc_core, test_matrix_SE3ftoMat4x4f);
  tcase_add_test(tc_core, test_matrix_SE3toMat4x4f);
  tcase_add_test(tc_core, test_matrix_dc);
  tcase_add_test(tc_core, test_matrix_SE3_setIdentity);
  tcase_add_test(tc_core, test_matrix_SE3f_setIdentity);
  tcase_add_test(tc_core, test_matrix_mat3x3_setIdentity);
  tcase_add_test(tc_core, test_matrix_mat3x3f_setIdentity);
  tcase_add_test(tc_core, test_matrix_mat4x4_setIdentity);
  tcase_add_test(tc_core, test_matrix_mat4x4f_setIdentity);
  tcase_add_test(tc_core, test_matrix_SE3_setAndGet);
  tcase_add_test(tc_core, test_matrix_SE3f_setAndGet);
  tcase_add_test(tc_core, test_matrix_mat3x3_setAndGet);
  tcase_add_test(tc_core, test_matrix_mat3x3f_setAndGet);
  tcase_add_test(tc_core, test_matrix_mat4x4_setAndGet);
  tcase_add_test(tc_core, test_matrix_mat4x4f_setAndGet);
  tcase_add_test(tc_core, test_matrix_mat3x3_setFromVectors);
  tcase_add_test(tc_core, test_matrix_mat3x3f_setFromVectors);
  tcase_add_test(tc_core, test_matrix_SE3_getTranslation);
  tcase_add_test(tc_core, test_matrix_SE3f_getTranslation);
  tcase_add_test(tc_core, test_matrix_mat3x3_copy);
  tcase_add_test(tc_core, test_matrix_mat3x3f_copy);
  tcase_add_test(tc_core, test_matrix_mat3x3_invert);
  tcase_add_test(tc_core, test_matrix_mat3x3f_invert);
  tcase_add_test(tc_core, test_matrix_mat3x3_inverted);
  tcase_add_test(tc_core, test_matrix_mat3x3f_inverted);
  tcase_add_test(tc_core, test_matrix_mat3x3_transpose);
  tcase_add_test(tc_core, test_matrix_mat3x3f_transpose);
  tcase_add_test(tc_core, test_matrix_mat3x3_transposed);
  tcase_add_test(tc_core, test_matrix_mat3x3f_transposed);
  tcase_add_test(tc_core, test_matrix_SE3_copy);
  tcase_add_test(tc_core, test_matrix_SE3f_copy);
  tcase_add_test(tc_core, test_matrix_SE3_inverted);
  tcase_add_test(tc_core, test_matrix_SE3f_inverted);
  tcase_add_test(tc_core, test_matrix_SE3_invert);
  tcase_add_test(tc_core, test_matrix_SE3f_invert);
  tcase_add_test(tc_core, test_matrix_mat3x3_getEulerAngles);
  tcase_add_test(tc_core, test_matrix_mat3x3_setFromEulerAngles);
  tcase_add_test(tc_core, test_matrix_mat4x4_getTranslation);
  tcase_add_test(tc_core, test_matrix_mat4x4f_getTranslation);
  tcase_add_test(tc_core, test_matrix_mat3x3_determinant);
  tcase_add_test(tc_core, test_matrix_mat3x3f_determinant);
  tcase_add_test(tc_core, test_matrix_mat3x3_getAxisAngle);
  tcase_add_test(tc_core, test_matrix_mat3x3f_getAxisAngle);
  tcase_add_test(tc_core, test_matrix_product_drv);
  tcase_add_test(tc_core, test_matrix_product_drvf);
  tcase_add_test(tc_core, test_matrix_product_drvinv);
  tcase_add_test(tc_core, test_matrix_product_drvinvf);
  tcase_add_test(tc_core, test_matrix_product_drinvv);
  tcase_add_test(tc_core, test_matrix_product_drinvvf);
  tcase_add_test(tc_core, test_matrix_product_drinvvinv);
  tcase_add_test(tc_core, test_matrix_product_drinvvinvf);
  tcase_add_test(tc_core, test_matrix_product_drvtx);
  tcase_add_test(tc_core, test_matrix_product_drvtxf);
  tcase_add_test(tc_core, test_matrix_product_drvtxinv);
  tcase_add_test(tc_core, test_matrix_product_drvtxinvf);
  tcase_add_test(tc_core, test_matrix_product_drvty);
  tcase_add_test(tc_core, test_matrix_product_drvtyf);
  tcase_add_test(tc_core, test_matrix_product_drvtyinv);
  tcase_add_test(tc_core, test_matrix_product_drvtyinvf);
  tcase_add_test(tc_core, test_matrix_product_drvtz);
  tcase_add_test(tc_core, test_matrix_product_drvtzf);
  tcase_add_test(tc_core, test_matrix_product_drvtzinv);
  tcase_add_test(tc_core, test_matrix_product_drvtzinvf);
  tcase_add_test(tc_core, test_matrix_product_drdrx);
  tcase_add_test(tc_core, test_matrix_product_drdrxf);
  tcase_add_test(tc_core, test_matrix_product_drdrxinv);
  tcase_add_test(tc_core, test_matrix_product_drdrxinvf);
  tcase_add_test(tc_core, test_matrix_product_drdry);
  tcase_add_test(tc_core, test_matrix_product_drdryf);
  tcase_add_test(tc_core, test_matrix_product_drdryinv);
  tcase_add_test(tc_core, test_matrix_product_drdryinvf);
  tcase_add_test(tc_core, test_matrix_product_drdrz);
  tcase_add_test(tc_core, test_matrix_product_drdrzf);
  tcase_add_test(tc_core, test_matrix_product_drdrzinv);
  tcase_add_test(tc_core, test_matrix_product_drdrzinvf);
  suite_add_tcase(s, tc_core);
  return s;
}

int main(void) {
  int number_failed;
  Suite* s = matrix_suite();
  SRunner* sr = srunner_create(s);
  srunner_set_log(sr, "matrix_test.log");
  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
