#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "inertia.h"
#include "util.h"

/*!
@brief double precision spatial inertia matrix constructor
@param m [in] mass
@param t [in] vector from reference frame to body centroid
@param Ixx [in] centroidal inertia term, expressed in reference frame
coordinates
@param Iyy [in] centroidal inertia term, expressed in reference frame
coordinates
@param Izz [in] centroidal inertia term, expressed in reference frame
coordinates
@param Ixy [in] centroidal inertia term, expressed in reference frame
coordinates
@param Iyz [in] centroidal inertia term, expressed in reference frame
coordinates
@param Izx [in] centroidal inertia term, expressed in reference frame
coordinates
*/
ksl_inertia_t ksl_inertia(const double m, const ksl_vec3_t t, const double Ixx,
                          const double Iyy, const double Izz, const double Ixy,
                          const double Iyz, const double Izx) {
  ksl_inertia_t i;
  i.m = m;
  i.mt.x = m * t.x;
  i.mt.y = m * t.y;
  i.mt.z = m * t.z;
  i.Ixx = Ixx;
  i.Iyy = Iyy;
  i.Izz = Izz;
  i.Ixy = Ixy;
  i.Iyz = Iyz;
  i.Izx = Izx;
  return i;
}

/*!
@brief single precision spatial inertia matrix constructor
@param m [in] mass
@param t [in] vector from reference frame to body centroid
@param Ixx [in] centroidal inertia term, expressed in reference frame
coordinates
@param Iyy [in] centroidal inertia term, expressed in reference frame
coordinates
@param Izz [in] centroidal inertia term, expressed in reference frame
coordinates
@param Ixy [in] centroidal inertia term, expressed in reference frame
coordinates
@param Iyz [in] centroidal inertia term, expressed in reference frame
coordinates
@param Izx [in] centroidal inertia term, expressed in reference frame
coordinates
*/
ksl_inertiaf_t ksl_inertiaf(const float m, const ksl_vec3f_t t, const float Ixx,
                            const float Iyy, const float Izz, const float Ixy,
                            const float Iyz, const float Izx) {
  ksl_inertiaf_t i;
  i.m = m;
  i.mt.x = m * t.x;
  i.mt.y = m * t.y;
  i.mt.z = m * t.z;
  i.Ixx = Ixx;
  i.Iyy = Iyy;
  i.Izz = Izz;
  i.Ixy = Ixy;
  i.Iyz = Iyz;
  i.Izx = Izx;
  return i;
}

inline ksl_inertia_t* ksl_inertia_alloc(int n) {
  return calloc(n, sizeof(ksl_inertia_t));
}

inline ksl_inertiaf_t* ksl_inertiaf_alloc(int n) {
  return calloc(n, sizeof(ksl_inertiaf_t));
}

inline void ksl_inertia_copy(const ksl_inertia_t* restrict inertia_i,
                             ksl_inertia_t* restrict inertia_o) {
  memcpy(inertia_o, inertia_i, sizeof(ksl_inertia_t));
}

inline void ksl_inertiaf_copy(const ksl_inertiaf_t* restrict inertia_i,
                              ksl_inertiaf_t* restrict inertia_o) {
  memcpy(inertia_o, inertia_i, sizeof(ksl_inertiaf_t));
}

/*!
@brief rotate a double precision centroidal inertia matrix

  I_ff = CoAd(D(R_{fc})) * Icc * Ad(D(R_{fc}^{-1}))

*/
inline void ksl_inertia_rotated(const ksl_inertia_t* restrict inertia_i,
                                const ksl_mat3x3_t* restrict r,
                                ksl_inertia_t* restrict inertia_o) {

  ksl_vec3_t temp;

  temp.at[0] = r->m00 * inertia_i->Ixx;
  temp.at[0] += r->m01 * inertia_i->Ixy;
  temp.at[0] += r->m02 * inertia_i->Izx;
  temp.at[1] = r->m00 * inertia_i->Ixy;
  temp.at[1] += r->m01 * inertia_i->Iyy;
  temp.at[1] += r->m02 * inertia_i->Iyz;
  temp.at[2] = r->m00 * inertia_i->Izx;
  temp.at[2] += r->m01 * inertia_i->Iyz;
  temp.at[2] += r->m02 * inertia_i->Izz;

  inertia_o->Ixx = temp.x * r->m00 + temp.y * r->m01 + temp.z * r->m02;
  inertia_o->Ixy = temp.x * r->m10 + temp.y * r->m11 + temp.z * r->m12;
  inertia_o->Izx = temp.x * r->m20 + temp.y * r->m21 + temp.z * r->m22;

  temp.at[0] = r->m10 * inertia_i->Ixx;
  temp.at[0] += r->m11 * inertia_i->Ixy;
  temp.at[0] += r->m12 * inertia_i->Izx;
  temp.at[1] = r->m10 * inertia_i->Ixy;
  temp.at[1] += r->m11 * inertia_i->Iyy;
  temp.at[1] += r->m12 * inertia_i->Iyz;
  temp.at[2] = r->m10 * inertia_i->Izx;
  temp.at[2] += r->m11 * inertia_i->Iyz;
  temp.at[2] += r->m12 * inertia_i->Izz;

  inertia_o->Iyy = temp.x * r->m10 + temp.y * r->m11 + temp.z * r->m12;
  inertia_o->Iyz = temp.x * r->m20 + temp.y * r->m21 + temp.z * r->m22;

  temp.at[0] = r->m20 * inertia_i->Ixx;
  temp.at[0] += r->m21 * inertia_i->Ixy;
  temp.at[0] += r->m22 * inertia_i->Izx;
  temp.at[1] = r->m20 * inertia_i->Ixy;
  temp.at[1] += r->m21 * inertia_i->Iyy;
  temp.at[1] += r->m22 * inertia_i->Iyz;
  temp.at[2] = r->m20 * inertia_i->Izx;
  temp.at[2] += r->m21 * inertia_i->Iyz;
  temp.at[2] += r->m22 * inertia_i->Izz;

  inertia_o->Izz = temp.x * r->m20 + temp.y * r->m21 + temp.z * r->m22;
  inertia_o->m = inertia_i->m;

  /* inertia_o->mt = Rfb * m_c * t_bc */
  ksl_product_drv(r, &inertia_i->mt, &inertia_o->mt);
}

/*!
@brief rotate a single precision centroidal inertia matrix

  I_ff = CoAd(D(R_{fc})) * Icc * Ad(D(R_{fc}^{-1}))

*/
inline void ksl_inertiaf_rotated(const ksl_inertiaf_t* restrict inertia_i,
                                 const ksl_mat3x3f_t* restrict r,
                                 ksl_inertiaf_t* restrict inertia_o) {

  ksl_vec3f_t temp;

  temp.at[0] = r->m00 * inertia_i->Ixx;
  temp.at[0] += r->m01 * inertia_i->Ixy;
  temp.at[0] += r->m02 * inertia_i->Izx;
  temp.at[1] = r->m00 * inertia_i->Ixy;
  temp.at[1] += r->m01 * inertia_i->Iyy;
  temp.at[1] += r->m02 * inertia_i->Iyz;
  temp.at[2] = r->m00 * inertia_i->Izx;
  temp.at[2] += r->m01 * inertia_i->Iyz;
  temp.at[2] += r->m02 * inertia_i->Izz;

  inertia_o->Ixx = temp.x * r->m00 + temp.y * r->m01 + temp.z * r->m02;
  inertia_o->Ixy = temp.x * r->m10 + temp.y * r->m11 + temp.z * r->m12;
  inertia_o->Izx = temp.x * r->m20 + temp.y * r->m21 + temp.z * r->m22;

  temp.at[0] = r->m10 * inertia_i->Ixx;
  temp.at[0] += r->m11 * inertia_i->Ixy;
  temp.at[0] += r->m12 * inertia_i->Izx;
  temp.at[1] = r->m10 * inertia_i->Ixy;
  temp.at[1] += r->m11 * inertia_i->Iyy;
  temp.at[1] += r->m12 * inertia_i->Iyz;
  temp.at[2] = r->m10 * inertia_i->Izx;
  temp.at[2] += r->m11 * inertia_i->Iyz;
  temp.at[2] += r->m12 * inertia_i->Izz;

  inertia_o->Iyy = temp.x * r->m10 + temp.y * r->m11 + temp.z * r->m12;
  inertia_o->Iyz = temp.x * r->m20 + temp.y * r->m21 + temp.z * r->m22;

  temp.at[0] = r->m20 * inertia_i->Ixx;
  temp.at[0] += r->m21 * inertia_i->Ixy;
  temp.at[0] += r->m22 * inertia_i->Izx;
  temp.at[1] = r->m20 * inertia_i->Ixy;
  temp.at[1] += r->m21 * inertia_i->Iyy;
  temp.at[1] += r->m22 * inertia_i->Iyz;
  temp.at[2] = r->m20 * inertia_i->Izx;
  temp.at[2] += r->m21 * inertia_i->Iyz;
  temp.at[2] += r->m22 * inertia_i->Izz;

  inertia_o->Izz = temp.x * r->m20 + temp.y * r->m21 + temp.z * r->m22;
  inertia_o->m = inertia_i->m;

  /* inertia_o->mt = Rfb * m_c * t_bc */
  ksl_product_drvf(r, &inertia_i->mt, &inertia_o->mt);
}

/*!
@brief rotate a double precision centroidal inertia matrix in place

  I_ff = CoAd(D(R_{fc})) * Icc * Ad(D(R_{fc}^{-1}))

*/
inline void ksl_inertia_rotate(ksl_inertia_t* restrict inertia_i,
                               const ksl_mat3x3_t* restrict r) {
  ksl_inertia_t inertia_o;
  ksl_inertia_rotated(inertia_i, r, &inertia_o);
  ksl_inertia_copy(&inertia_o, inertia_i);
}

/*!
@brief rotate a single precision centroidal inertia matrix in place

  I_ff = CoAd(D(R_{fc})) * Icc * Ad(D(R_{fc}^{-1}))

*/
inline void ksl_inertiaf_rotate(ksl_inertiaf_t* restrict inertia_i,
                                const ksl_mat3x3f_t* restrict r) {
  ksl_inertiaf_t inertia_o;
  ksl_inertiaf_rotated(inertia_i, r, &inertia_o);
  ksl_inertiaf_copy(&inertia_o, inertia_i);
}

/*!
@brief translate a double precision centroidal inertia matrix

  I_ff = CoAd(D(t_{fc})) * Icc * Ad(D(t_{fc}^{-1}))
*/
void ksl_inertia_translated(const ksl_inertia_t* restrict inertia_i,
                            const ksl_vec3_t* restrict t,
                            ksl_inertia_t* restrict inertia_o) {
  ksl_axpy_vv(inertia_i->m, t, &inertia_o->mt);

  ksl_vec3_t temp;
  memcpy(&temp, &inertia_o->mt, sizeof(ksl_vec3_t));
  assert(inertia_i->m > 0);
  ksl_vec3_scale(&temp, 1 / inertia_i->m);

  inertia_o->Ixy -= inertia_o->mt.x * temp.y;
  inertia_o->Iyz -= inertia_o->mt.y * temp.z;
  inertia_o->Izx -= inertia_o->mt.z * temp.x;

  ksl_vec3_t temp2;
  temp2.x = inertia_o->mt.x * temp.x;
  temp2.y = inertia_o->mt.y * temp.y;
  temp2.z = inertia_o->mt.z * temp.z;

  inertia_o->Ixx += temp2.y + temp2.z;
  inertia_o->Iyy += temp2.x + temp2.z;
  inertia_o->Izz += temp2.x + temp2.y;
}

/*!
@brief translate a single precision centroidal inertia matrix

  I_ff = CoAd(D(t_{fc})) * Icc * Ad(D(t_{fc}^{-1}))
*/
void ksl_inertiaf_translated(const ksl_inertiaf_t* restrict inertia_i,
                             const ksl_vec3f_t* restrict t,
                             ksl_inertiaf_t* restrict inertia_o) {
  ksl_axpy_vvf(inertia_i->m, t, &inertia_o->mt);

  ksl_vec3f_t temp;
  memcpy(&temp, &inertia_o->mt, sizeof(ksl_vec3_t));
  assert(inertia_i->m > 0);
  ksl_vec3f_scale(&temp, 1 / inertia_i->m);

  inertia_o->Ixy -= inertia_o->mt.x * temp.y;
  inertia_o->Iyz -= inertia_o->mt.y * temp.z;
  inertia_o->Izx -= inertia_o->mt.z * temp.x;

  ksl_vec3_t temp2;
  temp2.x = inertia_o->mt.x * temp.x;
  temp2.y = inertia_o->mt.y * temp.y;
  temp2.z = inertia_o->mt.z * temp.z;

  inertia_o->Ixx += temp2.y + temp2.z;
  inertia_o->Iyy += temp2.x + temp2.z;
  inertia_o->Izz += temp2.x + temp2.y;
}

/*!
@brief translate a double precision centroidal inertia matrix in place

  I_ff = CoAd(D(t_{fc})) * Icc * Ad(D(t_{fc}^{-1}))
*/
inline void ksl_inertia_translate(ksl_inertia_t* restrict inertia_i,
                                  const ksl_vec3_t* restrict t) {
  ksl_inertia_t inertia_o;
  ksl_inertia_translated(inertia_i, t, &inertia_o);
  ksl_inertia_copy(&inertia_o, inertia_i);
}

/*!
@brief translate a single precision centroidal inertia matrix in place

  I_ff = CoAd(D(t_{fc})) * Icc * Ad(D(t_{fc}^{-1}))
*/
inline void ksl_inertiaf_translate(ksl_inertiaf_t* restrict inertia_i,
                                   const ksl_vec3f_t* restrict t) {
  ksl_inertiaf_t inertia_o;
  ksl_inertiaf_translated(inertia_i, t, &inertia_o);
  ksl_inertiaf_copy(&inertia_o, inertia_i);
}

/*!
@brief spatial transform a double precision centroidal inertia matrix

  I_ff = CoAd(D_{fc}) * Icc * Ad(D_{fc}^{-1})
*/
inline void ksl_inertia_transformed(const ksl_inertia_t* restrict inertia_i,
                                    const ksl_SE3_t* restrict d,
                                    ksl_inertia_t* restrict inertia_o) {
  ksl_inertia_rotated(inertia_i, &d->R, inertia_o);
  ksl_axpy_vv(inertia_i->m, &d->t, &inertia_o->mt);

  ksl_vec3_t t;
  ksl_vec3_copy(&inertia_o->mt, &t);
  ksl_vec3_scale(&t, 1 / inertia_i->m);

  inertia_o->Ixy -= inertia_o->mt.x * t.y;
  inertia_o->Iyz -= inertia_o->mt.y * t.z;
  inertia_o->Izx -= inertia_o->mt.z * t.x;

  ksl_vec3_t temp2;
  temp2.x = inertia_o->mt.x * t.x;
  temp2.y = inertia_o->mt.y * t.y;
  temp2.z = inertia_o->mt.z * t.z;

  inertia_o->Ixx += temp2.y + temp2.z;
  inertia_o->Iyy += temp2.x + temp2.z;
  inertia_o->Izz += temp2.x + temp2.y;
}

/*!
@brief spatial transform a single precision centroidal inertia matrix

  I_ff = CoAd(D_{fc}) * Icc * Ad(D_{fc}^{-1})
*/
inline void ksl_inertiaf_transformed(const ksl_inertiaf_t* restrict inertia_i,
                                     const ksl_SE3f_t* restrict d,
                                     ksl_inertiaf_t* restrict inertia_o) {
  ksl_inertiaf_rotated(inertia_i, &d->R, inertia_o);
  ksl_axpy_vvf(inertia_i->m, &d->t, &inertia_o->mt);

  ksl_vec3f_t t;
  ksl_vec3f_copy(&inertia_o->mt, &t);
  ksl_vec3f_scale(&t, 1 / inertia_i->m);

  inertia_o->Ixy -= inertia_o->mt.x * t.y;
  inertia_o->Iyz -= inertia_o->mt.y * t.z;
  inertia_o->Izx -= inertia_o->mt.z * t.x;

  ksl_vec3_t temp2;
  temp2.x = inertia_o->mt.x * t.x;
  temp2.y = inertia_o->mt.y * t.y;
  temp2.z = inertia_o->mt.z * t.z;

  inertia_o->Ixx += temp2.y + temp2.z;
  inertia_o->Iyy += temp2.x + temp2.z;
  inertia_o->Izz += temp2.x + temp2.y;
}

/*!
@brief translate a double precision centroidal inertia matrix in place

  I_ff = CoAd(D(t_{fc})) * Icc * Ad(D(t_{fc}^{-1}))
*/
inline void ksl_inertia_transform(ksl_inertia_t* restrict inertia_i,
                                  const ksl_SE3_t* restrict d) {
  ksl_inertia_t inertia_o;
  ksl_inertia_transformed(inertia_i, d, &inertia_o);
  ksl_inertia_copy(&inertia_o, inertia_i);
}

/*!
@brief translate a single precision centroidal inertia matrix in place

  I_ff = CoAd(D(t_{fc})) * Icc * Ad(D(t_{fc}^{-1}))
*/
inline void ksl_inertiaf_transform(ksl_inertiaf_t* restrict inertia_i,
                                   const ksl_SE3f_t* restrict d) {
  ksl_inertiaf_t inertia_o;
  ksl_inertiaf_transformed(inertia_i, d, &inertia_o);
  ksl_inertiaf_copy(&inertia_o, inertia_i);
}

/*!
@brief merge double precision child inertia (j) with parent inertia (i)
@param parent inertia [in/out] parent inertia to be merged
@param t_ic_i [in/out] vector from parent body reference frame i to body center
of mass frame c, expressed in body frame i
@param child inertia [in]
@param t_jc_j [in] vector from child body reference frame j to body center of
mass frame c, expressed in body frame j
@param D_ij [in] transformation from parent reference frame i to child reference
frame j
*/
void ksl_inertia_merge(ksl_inertia_t* restrict inertia_i,
                       ksl_vec3_t* restrict t_ic_i,
                       ksl_inertia_t* restrict inertia_j,
                       const ksl_vec3_t* restrict t_jc_j,
                       const ksl_SE3_t* restrict D_ij) {

  /* locate the combined center of mass */
  ksl_vec3_t* mt_ic_i = &inertia_i->mt;
  ksl_vec3_t t_jc_i = {{0.0, 0.0, 0.0}};
  ksl_vec3_t mt_jc_i = {{0.0, 0.0, 0.0}};
  ksl_vec3_t t_icbar = {{0.0, 0.0, 0.0}};
  ksl_vec3_t mt_icbar = {{0.0, 0.0, 0.0}};

  double m_i = inertia_i->m;
  double m_j = inertia_j->m;
  double m_ij = m_i + m_j;

  /* 1. Locate the combined center of mass t_icbar^i */

  /* a. compute R^ij * t_jc^j (second half of Eq. 1.6)
  store results in t_jc_i */
  ksl_product_drv(&D_ij->R, t_jc_j, &t_jc_i);

  /* b. compute t_ij^i + t_jc^i -> t_jc^i (eq. 1.6) */
  ksl_xpy_vv(&D_ij->t, &t_jc_i);

  /* c. compute mt_jc_i (eq 1.6) */
  ksl_axpy_vv(m_j, &t_jc_i, &mt_jc_i);

  /* d. compute quantity m_{i} * t_{ic}^i + m_{j} * t_{ic} (eq. 1.6) */
  ksl_add_vv(mt_ic_i, &mt_jc_i, &mt_icbar);

  /* e. compute t_icbar^i =
    m_{i} * t_{ic}^i + m_{j} * t_{ic} / (m_i + m_j) (eq 1.7) */
  assert(m_ij > 0.0);
  ksl_product_av(1.0 / m_ij, &mt_icbar, &t_icbar);

  /* 2. compute combined moment of inertia referenced to the combined center
  of mass */

  /* a. find vector t_ccbar^i from center of mass to new combined centroid
  (cbar) (eq. 1.10) */
  ksl_vec3_t t_ccbar = {{0.0, 0.0, 0.0}};
  ksl_subtract_vv(&t_icbar, t_ic_i, &t_ccbar);

  /* b. transform centroidal inertia matrix of parent body to new centroid
  (eq. 1.11) */
  ksl_inertia_t* M_i_cc_ii = inertia_i;
  ksl_inertia_t M_i_ccbar;
  ksl_inertia_translated(M_i_cc_ii, &t_ccbar, &M_i_ccbar);

  /* c. transform child body centroidal moment of inertia to new centroid
  (eq. 1.12 (a)) */
  ksl_inertia_t* M_j_cc_jj = inertia_j;
  ksl_inertia_t M_j_ccbar_temp;
  ksl_inertia_rotated(M_j_cc_jj, &D_ij->R, &M_j_ccbar_temp);

  /* equation 1.13 */
  ksl_vec3_t temp = {{0.0, 0.0, 0.0}};
  ksl_product_drv(&D_ij->R, t_jc_j, &temp);
  ksl_xpy_vv(&D_ij->t, &temp);
  ksl_subtract_vv(&t_icbar, &temp, &t_ccbar);

  ksl_inertia_t M_j_ccbar;
  ksl_inertia_translated(&M_j_ccbar_temp, &t_ccbar, &M_j_ccbar);

  memset(inertia_j, 0, sizeof(ksl_inertia_t));

  inertia_i->m = m_ij;
  inertia_i->Ixx = M_i_ccbar.Ixx + M_j_ccbar.Ixx;
  inertia_i->Iyy = M_i_ccbar.Iyy + M_j_ccbar.Iyy;
  inertia_i->Izz = M_i_ccbar.Izz + M_j_ccbar.Izz;
  inertia_i->Ixy = M_i_ccbar.Ixy + M_j_ccbar.Ixy;
  inertia_i->Iyz = M_i_ccbar.Iyz + M_j_ccbar.Iyz;
  inertia_i->Izx = M_i_ccbar.Izx + M_j_ccbar.Izx;
  ksl_product_av(m_ij, &t_icbar, &inertia_i->mt);
}

/*!
@brief merge single precision child inertia (j) with parent inertia (i)
@param parent inertia [in/out] parent inertia to be merged
@param t_ic_i [in/out] vector from parent body reference frame i to body center
of mass frame c, expressed in body frame i
@param child inertia [in]
@param t_jc_j [in] vector from child body reference frame j to body center of
mass frame c, expressed in body frame j
@param D_ij [in] transformation from parent reference frame i to child reference
frame j
*/
void ksl_inertiaf_merge(ksl_inertiaf_t* restrict inertia_i,
                        ksl_vec3f_t* restrict t_ic_i,
                        ksl_inertiaf_t* restrict inertia_j,
                        const ksl_vec3f_t* restrict t_jc_j,
                        const ksl_SE3f_t* restrict D_ij) {

  /* locate the combined center of mass */
  ksl_vec3f_t* mt_ic_i = &inertia_i->mt;
  ksl_vec3f_t t_jc_i = {{0.0, 0.0, 0.0}};
  ksl_vec3f_t mt_jc_i = {{0.0, 0.0, 0.0}};
  ksl_vec3f_t t_icbar = {{0.0, 0.0, 0.0}};
  ksl_vec3f_t mt_icbar = {{0.0, 0.0, 0.0}};

  float m_i = inertia_i->m;
  float m_j = inertia_j->m;
  float m_ij = m_i + m_j;

  /* 1. Locate the combined center of mass t_icbar^i */

  /* a. compute R^ij * t_jc^j (second half of Eq. 1.6)
  store results in t_jc_i */
  ksl_product_drvf(&D_ij->R, t_jc_j, &t_jc_i);

  /* b. compute t_ij^i + t_jc^i -> t_jc^i (eq. 1.6) */
  ksl_xpy_vvf(&D_ij->t, &t_jc_i);

  /* c. compute mt_jc_i (eq 1.6) */
  ksl_axpy_vvf(m_j, &t_jc_i, &mt_jc_i);

  /* d. compute quantity m_{i} * t_{ic}^i + m_{j} * t_{ic} (eq. 1.6) */
  ksl_add_vvf(mt_ic_i, &mt_jc_i, &mt_icbar);

  /* e. compute t_icbar^i =
    m_{i} * t_{ic}^i + m_{j} * t_{ic} / (m_i + m_j) (eq 1.7) */
  assert(m_ij > 0.0);
  ksl_product_avf(1.0 / m_ij, &mt_icbar, &t_icbar);

  /* 2. compute combined moment of inertia referenced to the combined center
  of mass */

  /* a. find vector t_ccbar^i from center of mass to new combined centroid
  (cbar) (eq. 1.10) */
  ksl_vec3f_t t_ccbar = {{0.0, 0.0, 0.0}};
  ksl_subtract_vvf(&t_icbar, t_ic_i, &t_ccbar);

  /* b. transform centroidal inertia matrix of parent body to new centroid
  (eq. 1.11) */
  ksl_inertiaf_t* M_i_cc_ii = inertia_i;
  ksl_inertiaf_t M_i_ccbar;
  ksl_inertiaf_translated(M_i_cc_ii, &t_ccbar, &M_i_ccbar);

  /* c. transform child body centroidal moment of inertia to new centroid
  (eq. 1.12 (a)) */
  ksl_inertiaf_t* M_j_cc_jj = inertia_j;
  ksl_inertiaf_t M_j_ccbar_temp;
  ksl_inertiaf_rotated(M_j_cc_jj, &D_ij->R, &M_j_ccbar_temp);

  /* equation 1.13 */
  ksl_vec3f_t temp = {{0.0, 0.0, 0.0}};
  ksl_product_drvf(&D_ij->R, t_jc_j, &temp);
  ksl_xpy_vvf(&D_ij->t, &temp);
  ksl_subtract_vvf(&t_icbar, &temp, &t_ccbar);

  ksl_inertiaf_t M_j_ccbar;
  ksl_inertiaf_translated(&M_j_ccbar_temp, &t_ccbar, &M_j_ccbar);

  memset(inertia_j, 0, sizeof(ksl_inertiaf_t));

  inertia_i->m = m_ij;
  inertia_i->Ixx = M_i_ccbar.Ixx + M_j_ccbar.Ixx;
  inertia_i->Iyy = M_i_ccbar.Iyy + M_j_ccbar.Iyy;
  inertia_i->Izz = M_i_ccbar.Izz + M_j_ccbar.Izz;
  inertia_i->Ixy = M_i_ccbar.Ixy + M_j_ccbar.Ixy;
  inertia_i->Iyz = M_i_ccbar.Iyz + M_j_ccbar.Iyz;
  inertia_i->Izx = M_i_ccbar.Izx + M_j_ccbar.Izx;
  ksl_product_avf(m_ij, &t_icbar, &inertia_i->mt);
}

/*!
@brief factor double precision inertia matrix

@todo finish this

inertia matrix is overwritten with its factors

4 cases:
  * rank 2+, full rank, positive definite, or rank 2, positive semidefinite
  * rank 1, positive semidefinite
  * rank 0, do nothing

  @return rank of inertia matrix
*/
int ksl_inertia_factor(ksl_inertia_t* restrict inertia) {

  ksl_mat3x3_t M2;
  M2.as_array[0][0] = inertia->Ixx;
  M2.as_array[1][1] = inertia->Iyy;
  M2.as_array[2][2] = inertia->Izz;
  M2.as_array[0][1] = M2.as_array[1][0] = inertia->Ixy;
  M2.as_array[0][2] = M2.as_array[2][0] = inertia->Izx;
  M2.as_array[1][2] = M2.as_array[2][1] = inertia->Iyz;

  ksl_mat3x3_t M2WRTReference;
  ksl_mat3x3_copy(&M2, &M2WRTReference);
  double f[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  ksl_vec3i_t indices = {{0.0, 0.0, 0.0}};

  /* @todo this needs to be updated to use ksl functions */

  // if(body->frames->referenceFrame != body->frames->origin) {
  //   double temp[3][3];
  //   kinematics_ProductRInverseR(
  //     temp, body->frames->referenceFrame->rotationWRTOrigin, M2);
  //   kinematics_ProductRR(M2WRTReference, temp,
  //                        body->frames->referenceFrame->rotationWRTOrigin);
  //   kinematics_SubtractVectors(
  //     temp[0], body->inertia->centerOfMass,
  //     body->frames->referenceFrame->translationWRTOrigin);
  //   kinematics_ProductRInverseV(body->inertia->translationWRTReference,
  //                               body->frames->referenceFrame->rotationWRTOrigin,
  //                               temp[0]);
  // } else {
  //   memcpy(M2WRTReference, &M2, sizeof(M2));
  // }
  // fprintf(f, "3. Processing body: %s\n", body->name);

  /* compute 3 determinants */
  double determinants[3];
  determinants[0] =
    M2WRTReference.as_array[1][1] * M2WRTReference.as_array[2][2] -
    M2WRTReference.as_array[1][2] * M2WRTReference.as_array[2][1];
  determinants[1] =
    M2WRTReference.as_array[2][2] * M2WRTReference.as_array[0][0] -
    M2WRTReference.as_array[2][0] * M2WRTReference.as_array[0][2];
  determinants[2] =
    M2WRTReference.as_array[0][0] * M2WRTReference.as_array[1][1] -
    M2WRTReference.as_array[0][1] * M2WRTReference.as_array[1][0];

  /* find largest value in determinants */
  int i = 1;
  int j = 2;
  int k = 0;

  double maxVal = determinants[0];

  if(determinants[1] > maxVal) {
    i = 2;
    j = 0;
    k = 1;
    maxVal = determinants[1];
  }

  if(determinants[2] > maxVal) {
    i = 0;
    j = 1;
    k = 2;
    maxVal = determinants[2];
  }

  if(M2WRTReference.as_array[j][j] > M2WRTReference.as_array[i][i]) {
    ksl_swapi(&i, &j);
  }
  // fprintf(f, "4. Processing body: %s\n", body->name);
  // if maxVal > epsilon, then we have a rank of 2 or 3
  if(maxVal > 1e-9) {
    // F[i,i]
    f[0] = sqrt(M2WRTReference.as_array[i][i]);

    // F[j,i]
    f[1] = M2WRTReference.as_array[j][i] / f[0];

    // F[k,i]
    f[2] = M2WRTReference.as_array[k][i] / f[0];

    // F[j,j]
    f[3] = sqrt(M2WRTReference.as_array[j][j] - f[1] * f[1]);

    /* F[k,j] */
    f[4] = (M2WRTReference.as_array[k][j] - f[2] * f[1]) / f[3];

    /* F[k,k] */
    f[5] = sqrt(M2WRTReference.as_array[k][k] - f[2] * f[2] - f[4] * f[4]);
    indices.x = i;
    indices.y = j;
    indices.z = k;

    /* assign function pointer for rank 2,3 case */
    return 2;
  } else {

    maxVal = 0;

    /* scan all entries to see if there is at least one value greater than
    epsilon */
    for(i = 0; i < 3; i++) {
      for(j = i; j < 3; j++) {
        maxVal = fmax(maxVal, fabs(M2WRTReference.as_array[i][j]));
      }
    }
    if(maxVal < 1e-9) {
      /*rank = 0 */
      return 0;
    } else {
      maxVal = 0;
      j = 0;
      // rank = 1
      // 1.find largest value on diagonal
      for(i = 0; i < 3; i++) {
        // fprintf(f, "M2 diagonal term: %f\n", M2WRTReference[i][i]);
        if(fabs(M2WRTReference.as_array[i][i]) > maxVal) {
          j = i;
          // fprintf(f, "\n\nj = %d\n\n", j);
          maxVal = fabs(M2WRTReference.as_array[i][i]);
        }
      }
      /* 2. select a column and normalize it
         3. compute sqrt of trace
         4. multipy sqrt of trace with unit vector */
      double sqrtMassDivideByNorm =
        sqrt((M2WRTReference.as_array[0][0] + M2WRTReference.as_array[1][1] +
              M2WRTReference.as_array[2][2]) /
             (M2WRTReference.as_array[0][j] * M2WRTReference.as_array[0][j] +
              M2WRTReference.as_array[1][j] * M2WRTReference.as_array[1][j] +
              M2WRTReference.as_array[2][j] * M2WRTReference.as_array[2][j]));

      f[0] = M2WRTReference.as_array[0][j] * sqrtMassDivideByNorm;
      f[1] = M2WRTReference.as_array[1][j] * sqrtMassDivideByNorm;
      f[2] = M2WRTReference.as_array[2][j] * sqrtMassDivideByNorm;

      /* return rank 1 */
      return 1;
    }
  }
  return 0;
}
