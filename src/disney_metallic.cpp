#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

class DisneyMetallic : public BSDF {
public:
    DisneyMetallic(const PropertyList &propList) {
        m_roughness = propList.getFloat("roughness");
        m_anisotropic = propList.getFloat("anisotropic");

        m_aspect = sqrt(1.f - 0.9 * m_anisotropic);
        m_alphax = max(0.0001, m_roughness * m_roughness / aspect);
        m_alphay = max(0.0001, m_roughness * m_roughness * aspect);

        /* base color */
        if(propList.has("albedo")) {
            PropertyList l;
            l.setColor("value", propList.getColor("albedo"));
            m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));
        }
    }

    /// Evaluate the Fresnel term Fm
    Color3f Fresnel(const BSDFQueryRecord& bRec, const Vector3f& wh) const {
        Color3f basecolor = m_albedo->eval(bRec.uv);

        return basecolor + (1.f - basecolor) * float(pow(1 - abs(wh.dot(bRec.wo)), 5));
    }

    /// Evaluate the Trowbridge-Reitz distribution Dm
    /// ref: pbrb - microfacet
    float TrowbridgeReitz(const Vector3f& wh) const {
        // parallel, no sampling
        if (Frame::cosTheta(wh) == 0.f) {
            return 0.f;
        }

        float tan2Theta = Frame::tanTheta(wh) * Frame::tanTheta;
        float cos2Theta = Frame::cosTheta(wh) * Frame::cosTheta(wh);
        float cos4Theta = cos2Theta * cos2Theta;

        float e = (Frame::cosPhi2(wh) / (m_alphax * m_alphax) +
                    Frame::sinPhi2(wh) / (m_alphay * m_alphay)) * tan2Theta;
        
        return INV_PI / (m_alphax * m_alphay * cos4Theta * (1 + e) * (1 + e));
    }

    /// Evaluate the Lambda function
    /// ref: pbrb - microfacet
    float Lambda(const Vector3f& w) const {
        // parallel, no microfacet seen
        if (Frame::cosTheta(w) == 0.f) {
            return 0.f;
        }

        float absTanTheta = abs(Frame::tanTheta(w));
        float cosPhi = Frame::cosPhi(w);
        float sinPhi = Frame::sinPhi(w);

        float alpha2 = cosPhi * cosPhi * m_alphax * m_alphax + sinPhi * sinPhi * m_alphay * m_alphay;
        float alpha2Tan2Theta = alpha2 * absTanTheta * absTanTheta;

        return (-1.f + sqrt(1.f + alpha2Tan2Theta)) / 2.f;
    }

    /// Evaluate Smith's shadowing-masking function G
    float SmithG(const Vector3f &w) const {
        return 1.f / (1.f + Lambda(w));
    }

    /// Evaluate the BRDF for the given pair of directions
    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        /* no reflection from the back */
        if (Frame::cosTheta(bRec.wo) <= 0.f || bRec.measure != ESolidAngle)
            return 0;

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        Color3f Fm = Fresnel(bRec, wh);
        float Dm = TrowbridgeReitz(wh);
        float Gm = SmithG(bRec.wi) * SmithG(bRec.wo);
        float cosTheta = Frame::cosTheta(wi);

        return Fm * Dm * Gm / (4 * cosTheta);
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    /// ref: https://cseweb.ucsd.edu/~tzli/cse272/wi2023/homework1.pdf
    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        if (Frame::cosTheta(bRec.wo) <= 0.f || bRec.measure != ESolidAngle)
            return 0;

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float D = TrowbridgeReitz(wh);
        float Gi = SmithG(bRec.wi);
        float Jh = 1.f / (4.f * abs(Frame::cosTheta(bRec.wi)));

        return D * Gi * Jh;
    }

    /// VNDF sample the BRDF
    /// ref: https://jcgt.org/published/0007/04/01/
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const override {
        bRec.eta = 1.f;
        bRec.measure = ESolidAngle;
        
        if (Frame::cosTheta(bRec.wi) <= 0.f)
            return 0;

        // transform the view direction to the hemisphere configuration
        Vector3f Vh = Vector3f(m_alphax * bRec.wi.x(), m_alphay * bRec.wi.y(), bRec.wi.z()).normalized();

        // orthonormal basis
        float lensq = Vh.x() * Vh.x() + Vh.y() * Vh.y();
        Vector3f T1 = lensq > 0 ? Vector3f(-Vh.y(), Vh.x(), 0) / sqrt(lensq) : Vector3f(1, 0, 0);
        Vector3f T2 = Vh.cross(T1);

        // parameterization
        float r = sqrt(_sample.x());
        float phi = 2.f * M_PI * _sample.y();
        float t1 = r * cos(phi);
        float t2 = r * sin(phi);
        float s = 0.5 * (1.0 + Vh.z());
        t2 = (1.0 - s) * sqrt(1.0 - t1 * t1) + s * t2;

        // reprojection onto hemisphere
        Vector3f Nh = t1 * T1 + t2 * T2 + sqrt(max(0.0, 1.0 - t1*t1 - t2*t2)) * Vh;

        // transforming normal back to ellipsoid configuration
        Vector3f Ne = Vector3f(m_alphax * Nh.x(), m_alphay * Nh.y(), max(0.0, Nh.z())).normalized();
        bRec.wo = ((2.f * bRec.wi.dot(Ne) * Ne) - bRec.wi).normalized();

        /* discard invalid samples */
        if (pdf(bRec) > 0) {
            return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);   // cos forshortening term
        } else {
            return 0;
        }
    }

    virtual std::string toString() const override {
        return tfm::format(
            "Metallic[\n"
            "  roughness = %f,\n"
            "  anisotropic = %f\n"
            "]",
            m_roughness,
            m_anisotropic,
        );
    }
private:
    Texture<Color3f> * m_albedo;
    float m_roughness;
    float m_anisotropic;
    float m_aspect;
    float m_alphax;
    float m_alphay;
};

NORI_REGISTER_CLASS(DisneyMetallic, "disney_metallic");
NORI_NAMESPACE_END
