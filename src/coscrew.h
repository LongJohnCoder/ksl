/**
@file
@author Kristopher Wehage, Roger Wehage
@brief Utilities to initialize and operate on coscrews
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

#ifndef _KSL_COSCREW_H_
#define _KSL_COSCREW_H_

#include "matrix.h"

/*!
@brief double precision coscrew (i.e. a linear operator on a screw), consisting
of a linear and angular vector pair in function space.
*/
typedef union ksl_coscrew_t {
  double at[6];
  struct {
    ksl_vec3_t lin; /*!< bound linear vector, e.g. force or linear
                         momentum */
    ksl_vec3_t ang; /*!< free angular vector, e.g. moment/torque
                         or angular momentum*/
  };
  struct {
    double m0, m1, m2, m3, m4, m5;
  };
} ksl_coscrew_t;

/*!
@brief single precision coscrew (i.e. a linear operator on a screw), consisting
of a linear and angular vector pair in function space.
*/
typedef union ksl_coscrewf_t {
  float at[6];
  struct {
    ksl_vec3f_t lin; /*!< bound linear vector, e.g. force or linear
                         momentum */
    ksl_vec3f_t ang; /*!< free angular vector, e.g. moment/torque
                         or angular momentum*/
  };
  struct {
    float m0, m1, m2, m3, m4, m5;
  };
} ksl_coscrewf_t;

ksl_coscrew_t ksl_coscrew(const double m0, const double m1, const double m2,
                          const double m3, const double m4, const double m5);

ksl_coscrewf_t ksl_coscrewf(const float m0, const float m1, const float m2,
                            const float m3, const float m4, const float m5);

ksl_coscrew_t* ksl_coscrew_alloc(int);

ksl_coscrewf_t* ksl_coscrewf_alloc(int);

double ksl_coscrew_norm(const ksl_coscrew_t* v);

float ksl_coscrewf_norm(const ksl_coscrewf_t* v);

void ksl_coscrew_normalize(ksl_coscrew_t* v);

void ksl_coscrewf_normalize(ksl_coscrewf_t* v);

void ksl_coscrew_scale(ksl_coscrew_t* c, const double k);

void ksl_coscrewf_scale(ksl_coscrewf_t* c, const float k);

void ksl_coscrew_copy(const ksl_coscrew_t* ci, ksl_coscrew_t* co);

void ksl_coscrewf_copy(const ksl_coscrewf_t* ci, ksl_coscrewf_t* co);

void ksl_coscrew_invert(ksl_coscrew_t* ci);

void ksl_coscrewf_invert(ksl_coscrewf_t* ci);

void ksl_coscrew_inverted(const ksl_coscrew_t* ci, ksl_coscrew_t* co);

void ksl_coscrewf_inverted(const ksl_coscrewf_t* ci, ksl_coscrewf_t* co);

void ksl_axpy_cc(const double, const ksl_coscrew_t*, ksl_coscrew_t*);

void ksl_axpy_ccf(const float, const ksl_coscrewf_t*, ksl_coscrewf_t*);

void ksl_xpy_cc(const ksl_coscrew_t*, ksl_coscrew_t*);

void ksl_xpy_ccf(const ksl_coscrewf_t*, ksl_coscrewf_t*);

void ksl_nxpy_cc(const ksl_coscrew_t*, ksl_coscrew_t*);

void ksl_nxpy_ccf(const ksl_coscrewf_t*, ksl_coscrewf_t*);

void ksl_product_ac(const double k, const ksl_coscrew_t* ci, ksl_coscrew_t* co);

void ksl_product_acf(const float k, const ksl_coscrewf_t* ci,
                     ksl_coscrewf_t* co);

void ksl_add_cc(const ksl_coscrew_t* c1i, const ksl_coscrew_t* c2i,
                ksl_coscrew_t* co);

void ksl_add_ccf(const ksl_coscrewf_t* c1i, const ksl_coscrewf_t* c2i,
                 ksl_coscrewf_t* co);

void ksl_subtract_cc(const ksl_coscrew_t* c1i, const ksl_coscrew_t* c2i,
                     ksl_coscrew_t* co);

void ksl_subtract_ccf(const ksl_coscrewf_t* c1i, const ksl_coscrewf_t* c2i,
                      ksl_coscrewf_t* co);

void ksl_add_cct(const ksl_coscrew_t* ci1, const ksl_coscrew_t* ci2,
                 ksl_coscrew_t* co);

void ksl_add_cctf(const ksl_coscrewf_t* ci1, const ksl_coscrewf_t* ci2,
                  ksl_coscrewf_t* co);

void ksl_hctx(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hctxf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hcty(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hctyf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hctz(const ksl_SE3_t* ri, ksl_coscrew_t* co);

void ksl_hctzf(const ksl_SE3f_t* ri, ksl_coscrewf_t* co);

void ksl_hcrx(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hcrxf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hcry(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hcryf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hcrz(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hcrzf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hctxinv(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hctxinvf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hctyinv(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hctyinvf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hctzinv(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hctzinvf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hcrxinv(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hcrxinvf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hcryinv(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hcryinvf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_hcrzinv(const ksl_SE3_t* Di, ksl_coscrew_t* co);

void ksl_hcrzinvf(const ksl_SE3f_t* Di, ksl_coscrewf_t* co);

void ksl_cross_sc(const ksl_screw_t* s1i, const ksl_coscrew_t* c2i,
                  ksl_coscrew_t* co);

void ksl_cross_scf(const ksl_screwf_t* s1i, const ksl_coscrewf_t* c2i,
                   ksl_coscrewf_t* co);

void ksl_cross_sca(const ksl_screw_t* s1i, const ksl_coscrew_t* c2i,
                   ksl_coscrew_t* co);

void ksl_cross_scaf(const ksl_screwf_t* s1i, const ksl_coscrewf_t* c2i,
                    ksl_coscrewf_t* co);

void ksl_product_CoAdrc(const ksl_mat3x3_t* ri, const ksl_coscrew_t* ci,
                        ksl_coscrew_t* co);

void ksl_product_CoAdrcf(const ksl_mat3x3f_t* ri, const ksl_coscrewf_t* ci,
                         ksl_coscrewf_t* co);

void ksl_product_CoAdrcinv(const ksl_mat3x3_t* ri, const ksl_coscrew_t* ci,
                           ksl_coscrew_t* co);

void ksl_product_CoAdrcinvf(const ksl_mat3x3f_t* ri, const ksl_coscrewf_t* ci,
                            ksl_coscrewf_t* co);

void ksl_product_CoAdrinvc(const ksl_mat3x3_t* ri, const ksl_coscrew_t* ci,
                           ksl_coscrew_t* co);

void ksl_product_CoAdrinvcf(const ksl_mat3x3f_t* ri, const ksl_coscrewf_t* ci,
                            ksl_coscrewf_t* co);

void ksl_product_CoAdtc(const ksl_vec3_t* ti, const ksl_coscrew_t* ci,
                        ksl_coscrew_t* co);

void ksl_product_CoAdtcf(const ksl_vec3f_t* ti, const ksl_coscrewf_t* ci,
                         ksl_coscrewf_t* co);

void ksl_product_CoAdtinvc(const ksl_vec3_t* ti, const ksl_coscrew_t* ci,
                           ksl_coscrew_t* co);

void ksl_product_CoAdtinvcf(const ksl_vec3f_t* ti, const ksl_coscrewf_t* ci,
                            ksl_coscrewf_t* co);

void ksl_product_CoAdtcinv(const ksl_vec3_t* ti, const ksl_coscrew_t* ci,
                           ksl_coscrew_t* co);

void ksl_product_CoAdtcinvf(const ksl_vec3f_t* ti, const ksl_coscrewf_t* ci,
                            ksl_coscrewf_t* co);

void ksl_product_CoAdc(const ksl_SE3_t* Di, const ksl_coscrew_t* ci,
                       ksl_coscrew_t* co);

void ksl_product_CoAdcf(const ksl_SE3f_t* Di, const ksl_coscrewf_t* ci,
                        ksl_coscrewf_t* co);

void ksl_product_CoAdcinv(const ksl_SE3_t* Di, const ksl_coscrew_t* ci,
                          ksl_coscrew_t* co);

void ksl_product_CoAdcinvf(const ksl_SE3f_t* Di, const ksl_coscrewf_t* ci,
                           ksl_coscrewf_t* co);

void ksl_product_CoAdinvc(const ksl_SE3_t* Di, const ksl_coscrew_t* ci,
                          ksl_coscrew_t* co);

void ksl_product_CoAdinvcf(const ksl_SE3f_t* Di, const ksl_coscrewf_t* ci,
                           ksl_coscrewf_t* co);

#endif
