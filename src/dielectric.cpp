/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    virtual Color3f eval(const BSDFQueryRecord &) const override {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    virtual float pdf(const BSDFQueryRecord &) const override {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override {
        // x determines reflection or refraction
        // scaled x, y is used to map to the true color
        bRec.measure = EDiscrete;
        float cosTheta = Frame::cosTheta(bRec.wi);
        float Fr = fresnel(cosTheta, m_extIOR, m_intIOR);

        if (sample.x() < Fr) {    // reflection
            /* Relative index of refraction: no change */
            bRec.eta = 1.f;

            /* Reflection in local coordinates */
            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                bRec.wi.z()
            );

            return 1;
        } else {    // refraction
            /* Relative index of refraction: eta2 / eta1 */
            Normal3f n;
            if (cosTheta < 0) {
                bRec.eta = m_extIOR / m_intIOR;
                n = Normal3f(0.f, 0.f, -1.f);
            } else if (cosTheta >= 0) {
                bRec.eta = m_intIOR / m_extIOR;
                n = Normal3f(0.f, 0.f, 1.f);
            }

            /* Refraction in local coordinates */
            float eta_inv = 1.f / bRec.eta;
            bRec.wo = -eta_inv * (bRec.wi - bRec.wi.dot(n) * n) - n * sqrtf(1 - eta_inv * eta_inv * (1 - bRec.wi.dot(n) * bRec.wi.dot(n)));
            bRec.wo = bRec.wo.normalized();

            return eta_inv * eta_inv;
        }
    }

    virtual std::string toString() const override {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
