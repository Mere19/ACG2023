#include <nori/medium.h>
#include <nori/shape.h>
#include <nori/volume.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class Heterogeneous : public Medium {
public:

    Heterogeneous(const PropertyList& propList) {
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

    bool sample_freepath(MediumQueryRecord& mRec, Sampler* sampler) const override{        
        // sample distance based on delta tracking
        Vector3f direction = -mRec.wi;
        Ray3f ray = Ray3f(mRec.ref, direction, 0, mRec.tMax);
        Color3f transmittance = Color3f(1);
        float mint, maxt;
        if (!m_shape->getBoundingBox().rayIntersect(ray, mint, maxt))
            return false;
        mint = std::max(mint, ray.mint);
        maxt = std::min(maxt, ray.maxt);
        float t = mint, densityP = 0;
        
        while (true) {
            t += - log(1 - sampler->next1D()) * m_invMaxDensity;
            if (t >= maxt)
                break;
            Point3f p = ray(t);
            densityP = m_extinction->lookupFloat(p) * m_scale;
            if (densityP * m_invMaxDensity > sampler->next1D()) {
                mRec.p = p;
                mRec.ret = m_albedo->lookupRGB(p);
                return true;
            }
        }
        mRec.p = mRec.ref + mRec.tMax * direction;
        mRec.ret = 1.f;
        return false;
    }

    Color3f eval(const MediumQueryRecord& mRec) const {
        return m_albedo->lookupFloat(mRec.ref);
    }



    void addChild(NoriObject* obj) {
        switch (obj->getClassType()) {
        case EPhaseFunction:
            if (m_phase)
                throw NoriException(
                    "Already have a phase function");
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
            break;
        default:
            throw NoriException("Invalid",
                classTypeName(obj->getClassType()));
        }
    }

    virtual bool isHeterogeneous() const override{ return true; }
    virtual void volGrid() {
        if (m_extinction->isPerlin()) {
            m_extinction->gridGeneration();
        }
        if (m_albedo->isPerlin()) {
            m_albedo->gridGeneration();
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
    float m_scale;
    float m_maxDensity;
    float m_invMaxDensity;
    Transform m_worldToGrid;
    bool m_isRGB;
};

NORI_REGISTER_CLASS(Heterogeneous, "heterogeneous");
NORI_NAMESPACE_END