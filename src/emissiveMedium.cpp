#include <nori/medium.h>
#include <nori/shape.h>
#include <nori/volume.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class EmissiveMedium : public Medium {
public:

    EmissiveMedium(const PropertyList& propList) {
        m_scale = propList.getFloat("scale", 1.f);

        m_maxDensity = m_scale;
        m_invMaxDensity = 1.0f / m_scale;
        m_isRGB = propList.getBoolean("isRGB", true);
        if (m_isRGB) {
            m_scale = propList.getFloat("scale", 1.f);
        }
        else {
            m_scale = propList.getFloat("scale", 1.f);
        }

    }

    // sample the freepath based on delta tracking
    bool sample_freepath(MediumQueryRecord& mRec, Sampler* sampler) const override {
        Vector3f direction = -mRec.wi;
        Ray3f ray(mRec.ref, direction, 0, mRec.tMax);
        float pdf_failure = 1.0f;
        float pdf_success = 1.0f;
        Color3f transmittance(1.0f);
        Color3f Le(0.0f);
        float mint, maxt;
        if (!m_shape->getBoundingBox().rayIntersect(ray, mint, maxt))
            return false;
        mint = std::max(mint, ray.mint);
        maxt = std::min(maxt, ray.maxt);
        float t = mint, densityP, albedoAtT = 0;
        while (true) {
            t += - log(1 - sampler->next1D()) * m_invMaxDensity;
            if (t >= maxt)
                break;

            Point3f p = ray(t);
            densityP = m_extinction->lookupFloat(p) * m_scale;
            transmittance = exp(-densityP * t);

            // combination of line and point integral
            //if ((sampler->next1D()) >= transmittance.x()) {
            //    Le += m_radiance->lookupFloat(p) * (1 - albedoAtT) * (1 - exp(-densityAtT * t));
            //}
            //else {
            //    Le += m_radiance->lookupFloat(p) * (1 - albedoAtT);
            //}
            if (densityP * m_invMaxDensity > sampler->next1D()) {
                mRec.p = p;
                mRec.ret = m_albedo->lookupRGB(p);
                mRec.radiance = eval_radiance(p);
                return true;
            }
        }
        mRec.p = mRec.ref + mRec.tMax * direction;
        mRec.ret = 1.f;
        return false;
    }


    /* evaluate the transmittance between mRec.ref and mRec.p*/
    Color3f Tr(const MediumQueryRecord& mRec, Sampler* sampler) const override {
        Ray3f ray = Ray3f(mRec.ref, (mRec.p - mRec.ref).normalized(), 0, (mRec.p - mRec.ref).norm());
        float mint, maxt;
        if (!m_shape->getBoundingBox().rayIntersect(ray, mint, maxt))
            return Color3f(1.f);
        mint = std::max(mint, ray.mint);
        maxt = std::min(maxt, ray.maxt);
        int nSamples = 2;
        float result = 0;
            float t = mint;
            while (true) {
                t += -log(1 - sampler->next1D()) * m_invMaxDensity;
                if (t >= maxt) {
                    result += 1;
                    break;
                }
                Point3f p = ray(t);
                float density = m_extinction->lookupFloat(p) * m_scale;

                if (density * m_invMaxDensity > sampler->next1D())
                    break;
            }
        return Color3f(result / nSamples);
    }

    virtual Color3f eval_radiance(Point3f& p) const override{
        return Color3f(m_radiance->lookupFloat(p));
    }

    /**
     * \brief sample a point inside the medium shape based on shape.sampleVolume()
     , return the radiance at that point
        only assume a given mRec.ref 
     * */
    virtual Color3f sample_radiance(MediumQueryRecord& mRec, Sampler* sampler) const override{ 
        Point3f sample = Point3f(sampler->next1D(), sampler->next1D(), sampler->next1D());
        ShapeQueryRecord sRec(mRec.ref);
        if (!m_shape)
            throw NoriException(
                "There is no shape attached to this medium!");
        m_shape->sampleVolume(sRec, sample);
        mRec.p = sRec.p;
        mRec.wi = (mRec.p - mRec.ref).normalized();
        mRec.shadowRay = Ray3f(mRec.ref, mRec.wi, Epsilon, (mRec.p - mRec.ref).norm());
        mRec.radiance_pdf = sRec.pdf;
        return eval_radiance(mRec.p);
    }

    void addChild(NoriObject* obj) {
        switch (obj->getClassType()) {
        case EPhaseFunction:
            if (m_phase)
                throw NoriException(
                    "Medium: tried to register multiple Phase functions!");
            m_phase = static_cast<PhaseFunction*>(obj);
            break;
        case EVolume:
            if (obj->getIdName() == "albedo") {
                m_albedo = static_cast<Volume*>(obj);
                m_albedo->setMedium(static_cast<Medium*>(this));
            }
            else if (obj->getIdName() == "sigma_t") {
                m_extinction = static_cast<Volume*>(obj);
                m_extinction->setMedium(static_cast<Medium*>(this));
            }
            else if (obj->getIdName() == "radiance") {
                m_radiance = static_cast<Volume*>(obj);
                m_radiance->setMedium(static_cast<Medium*>(this));
            }
            break;
        default:
            throw NoriException("Medium::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }

    virtual bool isHeterogeneous() const override { return true; }
    virtual bool isEmissive() const { return true; }
    virtual void volGrid() {
        if (m_extinction->isPerlin()) {
            m_extinction->gridGeneration();
        }
        if (m_albedo->isPerlin()) {
            m_albedo->gridGeneration();
        }
        if (m_radiance->isPerlin()) {
            m_radiance->gridGeneration();
        }

    }
    std::string toString() const override {
        return tfm::format(
            "HeterogeneousMedium[\n"
            "  scale = %s,\n"
            "  sigmaT = %s,\n"
            "  albedo = %s,\n"

            "]",
            m_scale,
            m_extinction->toString(),
            m_albedo->toString());
    }


protected:

    Volume* m_extinction = nullptr;
    Volume* m_albedo = nullptr;
    Volume* m_radiance = nullptr;
    float m_scale;
    float m_maxDensity;
    float m_invMaxDensity;
    bool m_isRGB;
};

NORI_REGISTER_CLASS(EmissiveMedium, "emissive");
NORI_NAMESPACE_END