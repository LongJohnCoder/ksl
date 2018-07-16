/**
@file
@author Kristopher Wehage, Roger Wehage
@brief Utilities to initialize and operate on quaternions
@date 2018
@copyright Kristopher Wehage 2018

@remark
THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
*/

#ifndef _KSL_QUATERNION_H_
#define _KSL_QUATERNION_H_

#include "matrix.h"
#include "vector.h"

/*!
@brief double precision quaternion data structure
*/
typedef union ksl_quaternion_t {
  double at[4];
  struct {
    union {
      ksl_vec3_t r;
      struct {
        double x;
        double y;
        double z;
      };
    };
    double w;
  };
} ksl_quaternion_t;

/*!
@brief single precision quaternion data structure
*/
typedef union ksl_quaternionf_t {
  float at[4];
  struct {
    union {
      ksl_vec3f_t r;
      struct {
        float x;
        float y;
        float z;
      };
    };
    float w;
  };
} ksl_quaternionf_t;

ksl_quaternion_t ksl_quaternion(const double x, const double y, const double z,
                                const double w);

ksl_quaternionf_t ksl_quaternionf(const float x, const float y, const float z,
                                  const float w);

ksl_quaternion_t ksl_quaternion_unit(const double x, const double y,
                                     const double z, const double w);

ksl_quaternionf_t ksl_quaternionf_unit(const float x, const float y,
                                       const float z, const float w);

// ksl_quaternion_t ksl_quaternion_axis_angle();

ksl_quaternion_t* ksl_quaternion_alloc(int);

ksl_quaternionf_t* ksl_quaternionf_alloc(int);

int ksl_quaternion_maxIndex(const ksl_quaternion_t* q);

int ksl_quaternionf_maxIndex(const ksl_quaternionf_t* q);

double ksl_dot_qq(const ksl_quaternion_t* q1, const ksl_quaternion_t* q2);

float ksl_dot_qqf(const ksl_quaternionf_t* q1, const ksl_quaternionf_t* q2);

void ksl_quaternion_normalize(ksl_quaternion_t* qi);

void ksl_quaternionf_normalize(ksl_quaternionf_t* qi);

void ksl_axpy_qq(double a, const ksl_quaternion_t* x, ksl_quaternion_t* y);

void ksl_axpy_qqf(float a, const ksl_quaternionf_t* x, ksl_quaternionf_t* y);

void ksl_xpy_qq(const ksl_quaternion_t* x, ksl_quaternion_t* y);

void ksl_xpy_qqf(const ksl_quaternionf_t* x, ksl_quaternionf_t* y);

void ksl_nxpy_qq(const ksl_quaternion_t* x, ksl_quaternion_t* y);

void ksl_nxpy_qqf(const ksl_quaternionf_t* x, ksl_quaternionf_t* y);

void ksl_mat3x3_to_quaternion(const ksl_mat3x3_t* ri, ksl_quaternion_t* qo);

void ksl_mat3x3f_to_quaternionf(const ksl_mat3x3f_t* ri, ksl_quaternionf_t* qo);

void ksl_quaternion_to_mat3x3(ksl_quaternion_t* qi, ksl_mat3x3_t* ro);

void ksl_quaternionf_to_mat3x3f(ksl_quaternionf_t* qi, ksl_mat3x3f_t* ro);

void ksl_slerp(const ksl_quaternion_t* q1i, const ksl_quaternion_t* q2i,
               double* t, ksl_quaternion_t* qo);

void ksl_slerpf(const ksl_quaternionf_t* q1i, const ksl_quaternionf_t* q2i,
                float* t, ksl_quaternionf_t* qo);

void ksl_squad(const ksl_quaternion_t qi[4], double* t, ksl_quaternion_t* qo);

void ksl_squadf(const ksl_quaternionf_t qi[4], float* t, ksl_quaternionf_t* qo);

void ksl_nlerp(const ksl_quaternion_t* q1i, const ksl_quaternion_t* q2i,
               const double* t, ksl_quaternion_t* qo);

void ksl_nlerpf(const ksl_quaternionf_t* q1i, const ksl_quaternionf_t* q2i,
                const float* t, ksl_quaternionf_t* qo);

void ksl_dpi_p_to_dq(const ksl_vec3_t* dpii, const ksl_quaternion_t* qi,
                     ksl_quaternion_t* dqo);

void ksl_dpi_pf_to_dqf(const ksl_vec3f_t* dpii, const ksl_quaternionf_t* qi,
                       ksl_quaternionf_t* dqo);

void ksl_dpi_c_to_dq(const ksl_vec3_t* dpii, const ksl_quaternion_t* qi,
                     ksl_quaternion_t* dqo);

void ksl_dpi_cf_to_dqf(const ksl_vec3f_t* dpii, const ksl_quaternionf_t* qi,
                       ksl_quaternionf_t* dqo);

void ksl_dq_to_dpi_p(const ksl_quaternion_t* qi, const ksl_quaternion_t* dqi,
                     ksl_vec3_t* dpio);

void ksl_dqf_to_dpi_pf(const ksl_quaternionf_t* qi,
                       const ksl_quaternionf_t* dqi, ksl_vec3f_t* dpio);

void ksl_dq_to_dpi_c(const ksl_quaternion_t* qi, const ksl_quaternion_t* dqi,
                     ksl_vec3_t* dpio);

void ksl_dqf_to_dpi_cf(const ksl_quaternionf_t* qi,
                       const ksl_quaternionf_t* dqi, ksl_vec3f_t* dpio);

inline void ksl_omega_p_to_dquaternion(const ksl_vec3_t* wi,
                                       const ksl_quaternion_t* qi,
                                       ksl_quaternion_t* dqo);

inline void ksl_omega_pf_to_dquaternionf(const ksl_vec3f_t* wi,
                                         const ksl_quaternionf_t* qi,
                                         ksl_quaternionf_t* dqo);

inline void ksl_omega_c_to_dquaternion(const ksl_vec3_t* wi,
                                       const ksl_quaternion_t* qi,
                                       ksl_quaternion_t* dqo);

inline void ksl_omega_cf_to_dquaternionf(const ksl_vec3f_t* wi,
                                         const ksl_quaternionf_t* qi,
                                         ksl_quaternionf_t* dqo);

void ksl_product_qq(const ksl_quaternion_t* q1i, const ksl_quaternion_t* q2i,
                    ksl_quaternion_t* qo);

void ksl_product_qqf(const ksl_quaternionf_t* q1i, const ksl_quaternionf_t* q2i,
                     ksl_quaternionf_t* qo);

void ksl_product_qconjq(const ksl_quaternion_t* q1i,
                        const ksl_quaternion_t* q2i, ksl_quaternion_t* qo);

void ksl_product_qconjqf(const ksl_quaternionf_t* q1i,
                         const ksl_quaternionf_t* q2i, ksl_quaternionf_t* qo);

void ksl_product_qqconj(const ksl_quaternion_t* q1i,
                        const ksl_quaternion_t* q2i, ksl_quaternion_t* qo);

void ksl_product_qqconjf(const ksl_quaternionf_t* q1i,
                         const ksl_quaternionf_t* q2i, ksl_quaternionf_t* qo);

void ksl_product_qconjqconj(const ksl_quaternion_t* q1i,
                            const ksl_quaternion_t* q2i, ksl_quaternion_t* qo);

void ksl_product_qconjqconjf(const ksl_quaternionf_t* q1i,
                             const ksl_quaternionf_t* q2i,
                             ksl_quaternionf_t* qo);

void ksl_product_qxq(const double dci[2], const ksl_quaternion_t* qi,
                     ksl_quaternion_t* qo);

void ksl_product_qxqf(const float dci[2], const ksl_quaternionf_t* qi,
                      ksl_quaternionf_t* qo);

void ksl_product_qqx(const ksl_quaternion_t* qi, const double dci[2],
                     ksl_quaternion_t* qo);

void ksl_product_qqxf(const ksl_quaternionf_t* qi, const float dci[2],
                      ksl_quaternionf_t* qo);

void ksl_product_qxconjq(const double dci[2], const ksl_quaternion_t* qi,
                         ksl_quaternion_t* qo);

void ksl_product_qxconjqf(const float dci[2], const ksl_quaternionf_t* qi,
                          ksl_quaternionf_t* qo);

void ksl_product_qconjqx(const ksl_quaternion_t* qi, const double dci[2],
                         ksl_quaternion_t* qo);

void ksl_product_qconjqxf(const ksl_quaternionf_t* qi, const float dci[2],
                          ksl_quaternionf_t* qo);

void ksl_product_qxqconj(const double dci[2], const ksl_quaternion_t* qi,
                         ksl_quaternion_t* qo);

void ksl_product_qxqconjf(const float dci[2], const ksl_quaternionf_t* qi,
                          ksl_quaternionf_t* qo);

void ksl_product_qqxconj(const ksl_quaternion_t* qi, const double dci[2],
                         ksl_quaternion_t* qo);

void ksl_product_qqxconjf(const ksl_quaternionf_t* qi, const float dci[2],
                          ksl_quaternionf_t* qo);

void ksl_product_qxconjqconj(const double dci[2], const ksl_quaternion_t* qi,
                             ksl_quaternion_t* qo);

void ksl_product_qxconjqconjf(const float dci[2], const ksl_quaternionf_t* qi,
                              ksl_quaternionf_t* qo);

void ksl_product_qconjqxconj(const ksl_quaternion_t* qi, const double dci[2],
                             ksl_quaternion_t* qo);

void ksl_product_qconjqxconjf(const ksl_quaternionf_t* qi, const float dci[2],
                              ksl_quaternionf_t* qo);

void ksl_product_qyq(const double dci[2], const ksl_quaternion_t* qi,
                     ksl_quaternion_t* qo);

void ksl_product_qyqf(const float dci[2], const ksl_quaternionf_t* qi,
                      ksl_quaternionf_t* qo);

void ksl_product_qqy(const ksl_quaternion_t* qi, const double dci[2],
                     ksl_quaternion_t* qo);

void ksl_product_qqyf(const ksl_quaternionf_t* qi, const float dci[2],
                      ksl_quaternionf_t* qo);

void ksl_product_qyconjq(const double dci[2], const ksl_quaternion_t* qi,
                         ksl_quaternion_t* qo);

void ksl_product_qyconjqf(const float dci[2], const ksl_quaternionf_t* qi,
                          ksl_quaternionf_t* qo);

void ksl_product_qconjqy(const ksl_quaternion_t* qi, const double dci[2],
                         ksl_quaternion_t* qo);

void ksl_product_qconjqyf(const ksl_quaternionf_t* qi, const float dci[2],
                          ksl_quaternionf_t* qo);

void ksl_product_qyqconj(const double dci[2], const ksl_quaternion_t* qi,
                         ksl_quaternion_t* qo);

void ksl_product_qyqconjf(const float dci[2], const ksl_quaternionf_t* qi,
                          ksl_quaternionf_t* qo);

void ksl_product_qqyconj(const ksl_quaternion_t* qi, const double dci[2],
                         ksl_quaternion_t* qo);

void ksl_product_qqyconjf(const ksl_quaternionf_t* qi, const float dci[2],
                          ksl_quaternionf_t* qo);

void ksl_product_qyconjqconj(const double dci[2], const ksl_quaternion_t* qi,
                             ksl_quaternion_t* qo);

void ksl_product_qyconjqconjf(const float dci[2], const ksl_quaternionf_t* qi,
                              ksl_quaternionf_t* qo);

void ksl_product_qconjqyconj(const ksl_quaternion_t* qi, const double dci[2],
                             ksl_quaternion_t* qo);

void ksl_product_qconjqyconjf(const ksl_quaternionf_t* qi, const float dci[2],
                              ksl_quaternionf_t* qo);

void ksl_product_qzq(const double dci[2], const ksl_quaternion_t* qi,
                     ksl_quaternion_t* qo);

void ksl_product_qzqf(const float dci[2], const ksl_quaternionf_t* qi,
                      ksl_quaternionf_t* qo);

void ksl_product_qqz(const ksl_quaternion_t* qi, const double dci[2],
                     ksl_quaternion_t* qo);

void ksl_product_qqzf(const ksl_quaternionf_t* qi, const float dci[2],
                      ksl_quaternionf_t* qo);

void ksl_product_qzconjq(const double dci[2], const ksl_quaternion_t* qi,
                         ksl_quaternion_t* qo);

void ksl_product_qzconjqf(const float dci[2], const ksl_quaternionf_t* qi,
                          ksl_quaternionf_t* qo);

void ksl_product_qconjqz(const ksl_quaternion_t* qi, const double dci[2],
                         ksl_quaternion_t* qo);

void ksl_product_qconjqzf(const ksl_quaternionf_t* qi, const float dci[2],
                          ksl_quaternionf_t* qo);

void ksl_product_qzqconj(const double dci[2], const ksl_quaternion_t* qi,
                         ksl_quaternion_t* qo);

void ksl_product_qzqconjf(const float dci[2], const ksl_quaternionf_t* qi,
                          ksl_quaternionf_t* qo);

void ksl_product_qqzconj(const ksl_quaternion_t* qi, const double dci[2],
                         ksl_quaternion_t* qo);

void ksl_product_qqzconjf(const ksl_quaternionf_t* qi, const float dci[2],
                          ksl_quaternionf_t* qo);

void ksl_product_qzconjqconj(const double dci[2], const ksl_quaternion_t* qi,
                             ksl_quaternion_t* qo);

void ksl_product_qzconjqconjf(const float dci[2], const ksl_quaternionf_t* qi,
                              ksl_quaternionf_t* qo);

void ksl_product_qconjqzconj(const ksl_quaternion_t* qi, const double dci[2],
                             ksl_quaternion_t* qo);

void ksl_product_qconjqzconjf(const ksl_quaternionf_t* qi, const float dci[2],
                              ksl_quaternionf_t* qo);

void ksl_product_qv(const ksl_quaternion_t* qi, const ksl_vec3_t* vi,
                    ksl_quaternion_t* qo);

void ksl_product_qvf(const ksl_quaternionf_t* qi, const ksl_vec3f_t* vi,
                     ksl_quaternionf_t* qo);

void ksl_product_qconjv(const ksl_quaternion_t* qi, const ksl_vec3_t* vi,
                        ksl_quaternion_t* qo);

void ksl_product_qconjvf(const ksl_quaternionf_t* qi, const ksl_vec3f_t* vi,
                         ksl_quaternionf_t* qo);

void ksl_product_vq(const ksl_vec3_t* vi, const ksl_quaternion_t* qi,
                    ksl_quaternion_t* qo);

void ksl_product_vqf(const ksl_vec3f_t* vi, const ksl_quaternionf_t* qi,
                     ksl_quaternionf_t* qo);

void ksl_product_vqconj(const ksl_vec3_t* vi, const ksl_quaternion_t* qi,
                        ksl_quaternion_t* qo);

void ksl_product_vqconjf(const ksl_vec3f_t* vi, const ksl_quaternionf_t* qi,
                         ksl_quaternionf_t* qo);

void ksl_product_uqinverseq(ksl_quaternion_t* q1i, ksl_quaternion_t* q2i,
                            ksl_quaternion_t* qo);

#endif
