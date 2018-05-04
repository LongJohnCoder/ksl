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

ksl_quaternion_t* ksl_quaternion_alloc(int);

ksl_quaternionf_t* ksl_quaternionf_alloc(int);

double ksl_quaternion_dot_qq(const ksl_quaternion_t*, const ksl_quaternion_t*);

float ksl_quaternionf_dot_qq(const ksl_quaternionf_t*,
                             const ksl_quaternionf_t*);

void ksl_quaternion_normalize(ksl_quaternion_t* qi);

void ksl_quaternionf_normalize(ksl_quaternionf_t* qi);

void ksl_axpy_qq(double, const ksl_quaternion_t*, ksl_quaternion_t*);

void ksl_axpy_qqf(float, const ksl_quaternionf_t*, ksl_quaternionf_t*);

void ksl_mat3x3_toQuaternion(const ksl_mat3x3_t* ri, ksl_quaternion_t* qo);

void ksl_mat3x3f_toQuaternion(const ksl_mat3x3f_t* ri, ksl_quaternionf_t* qo);

void ksl_quaternion_toMat3x3(ksl_quaternion_t* qi, ksl_mat3x3_t* ro);

void ksl_quaternionf_toMat3x3f(ksl_quaternionf_t* qi, ksl_mat3x3f_t* ro);

void ksl_quaternion_slerp(const ksl_quaternion_t* q1i,
                          const ksl_quaternion_t* q2i, double* t,
                          ksl_quaternion_t* qo);

void ksl_quaternionf_slerp(const ksl_quaternionf_t* q1i,
                           const ksl_quaternionf_t* q2i, float* t,
                           ksl_quaternionf_t* qo);

void ksl_quaternion_squad(const ksl_quaternion_t qi[4], double* t,
                          ksl_quaternion_t* qo);

void ksl_quaternionf_squad(const ksl_quaternionf_t qi[4], float* t,
                           ksl_quaternionf_t* qo);

void ksl_quaternion_nlerp(const ksl_quaternion_t* q1i,
                          const ksl_quaternion_t* q2i, const double* t,
                          ksl_quaternion_t* qo);

void ksl_quaternionf_nlerp(const ksl_quaternionf_t* q1i,
                           const ksl_quaternionf_t* q2i, const float* t,
                           ksl_quaternionf_t* qo);

void ksl_add_dpiToq(const ksl_vec3_t* dpi, ksl_quaternion_t* qi);

void ksl_subtract_dpiFromq(const ksl_vec3_t* dpi, ksl_quaternion_t* qi);

void ksl_omega_to_dquaternion(const ksl_vec3_t* w, const ksl_quaternion_t* qi,
                              ksl_quaternion_t* dqo);

void ksl_quaternion_product_qq(const ksl_quaternion_t* q1i,
                               const ksl_quaternion_t* q2i,
                               ksl_quaternion_t* qo);

void ksl_quaternionf_product_qq(const ksl_quaternionf_t* q1i,
                                const ksl_quaternionf_t* q2i,
                                ksl_quaternionf_t* qo);

void ksl_quaternion_product_qconjq(const ksl_quaternion_t* q1i,
                                   const ksl_quaternion_t* q2i,
                                   ksl_quaternion_t* qo);

void ksl_quaternionf_product_qconjq(const ksl_quaternionf_t* q1i,
                                    const ksl_quaternionf_t* q2i,
                                    ksl_quaternionf_t* qo);

void ksl_quaternion_product_qqconj(const ksl_quaternion_t* q1i,
                                   const ksl_quaternion_t* q2i,
                                   ksl_quaternion_t* qo);

void ksl_quaternionf_product_qqconj(const ksl_quaternionf_t* q1i,
                                    const ksl_quaternionf_t* q2i,
                                    ksl_quaternionf_t* qo);

void ksl_quaternion_product_qconjqconj(const ksl_quaternion_t* q1i,
                                       const ksl_quaternion_t* q2i,
                                       ksl_quaternion_t* qo);

void ksl_quaternionf_product_qconjqconj(const ksl_quaternionf_t* q1i,
                                        const ksl_quaternionf_t* q2i,
                                        ksl_quaternionf_t* qo);

void ksl_quaternion_product_uqinverseq(ksl_quaternion_t* q1i,
                                       ksl_quaternion_t* q2i,
                                       ksl_quaternion_t* qo);

void ksl_quaternion_product_qv(const ksl_quaternion_t* q1i,
                               const ksl_vec3_t* vi, ksl_quaternion_t* qo);

void ksl_quaternionf_product_qv(const ksl_quaternionf_t* q1i,
                                const ksl_vec3f_t* vi, ksl_quaternionf_t* qo);

#endif
