#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <nori/common.h>

#include <Eigen/Geometry>
#include <Eigen/LU>

NORI_NAMESPACE_BEGIN

class DisneySheen : public BSDF {
public:
    DisneySheen(const PropertyList &propList) {
        m_sheentint = propList.getFloat("sheen", 0.01f);

        /* base color */
        if(propList.has("albedo")) {
            PropertyList l;
            l.setColor("value", propList.getColor("albedo"));
            m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));
        }
    }

    /// Evaluate the Fresnel term Fm
    virtual Color3f eval(const BSDFQueryRecord& bRec) const override {
        Color3f basecolor = m_albedo->eval(bRec.uv);
        float luminance = basecolor.getLuminance();

        Color3f Ctint = luminance > 0 ? basecolor / luminance : Color3f(1.f);
        Color3f Csheen = (1 - m_sheentint) + m_sheentint * Ctint;
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        return Csheen * pow((1.f - abs(wh.dot(bRec.wo))), 5) * Frame::cosTheta(bRec.wo);
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
            "Sheen[\n"
            "  sheen = %f\n"
            "]",
            m_sheentint
        );
    }
private:
    Texture<Color3f> * m_albedo;
    float m_sheentint;
};

NORI_REGISTER_CLASS(DisneySheen, "disney_sheen");
NORI_NAMESPACE_END
