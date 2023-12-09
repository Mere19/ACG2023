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
#include <nori/warp.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

class DisneyDiffuse : public BSDF {
public:
    DisneyDiffuse(const PropertyList &propList) {
        m_subsurface = propList.getFloat("subsurface", 0.5f);
        m_roughness = propList.getFloat("roughness", 0.7f);

        /* base color */
        if(propList.has("albedo")) {
            PropertyList l;
            l.setColor("value", propList.getColor("albedo"));
            m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));
        }
    }

    /// Add texture for the albedo
    virtual void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case ETexture:
                if(obj->getIdName() == "albedo") {
                    if (m_albedo)
                        throw NoriException("There is already an albedo defined!");
                    m_albedo = static_cast<Texture<Color3f> *>(obj);
                }
                else {
                    throw NoriException("The name of this texture does not match any field!");
                }
                break;

            default:
                throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    /// @brief evaluate base diffusion
    /// @param bRec 
    /// @return evaluated color given bRec
    Color3f evalBase(const BSDFQueryRecord &bRec) const {
        /* compute FD90 */
        float FD90 = 0.5f + 2 * m_roughness * abs(bRec.wo.dot((bRec.wi + bRec.wo).normalized()));

        /* compute FD */
        float FD_wi = 1 + (FD90 - 1) * (1 - pow(abs(bRec.wi.z()), 5));

        return m_albedo->eval(bRec.uv) * INV_PI * FD_wi * FD90 * abs(bRec.wo.z());
    }

    /// @brief evaluate subsurface diffusion
    /// @param bRec 
    /// @return evaluated color given bRec
    Color3f evalSubsurface(const BSDFQueryRecord &bRec) const {
        /* compute FSS90 */
        float FSS90 = m_roughness * pow(abs(bRec.wo.dot((bRec.wi + bRec.wo).normalized())), 2);

        /* compute FSS_wi */
        float FSS_wi = (1 + (FSS90 - 1) * pow(1 - abs(bRec.wi.z()), 5));

        /* compute FSS_wo */
        float FSS_wo = (1 + (FSS90 - 1) * pow(1 - abs(bRec.wo.z()), 5));

        Color3f x = 1.25f * m_albedo->eval(bRec.uv) * INV_PI;
        Color3f y = FSS_wi * FSS_wo * (1.f / (abs(bRec.wi.z()) + abs(bRec.wo.z())) - 0.5f) + 0.5f;
        Color3f z = abs(bRec.wo.z());

        return x * y * z;
    }

    /// Evaluate the BRDF for the given pair of directions
    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Color3f base = evalBase(bRec);
        Color3f subsurface = evalSubsurface(bRec);

        return (1.f - m_subsurface) * base + m_subsurface * subsurface;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;


        /* Importance sampling density wrt. solid angles:
           cos(theta) / pi.

           Note that the directions in 'bRec' are in local coordinates,
           so Frame::cosTheta() actually just returns the 'z' component.
        */
        return INV_PI * Frame::cosTheta(bRec.wo);
    }

    /// Sample the BRDF
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override {
        if (Frame::cosTheta(bRec.wi) <= 0)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */
        bRec.wo = Warp::squareToCosineHemisphere(sample);

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        return eval(bRec) / pdf(bRec) * Frame::cosTheta(bRec.wo);
    }

    virtual std::string toString() const override {
        return tfm::format(
            "DisneyDiffuse[\n"
            "  alpha = %f,\n"
            "]",
            m_roughness
        );
    }
private:
    float m_subsurface;
    float m_roughness;
    Texture<Color3f> * m_albedo;
};

NORI_REGISTER_CLASS(DisneyDiffuse, "disney_diffuse");
NORI_NAMESPACE_END