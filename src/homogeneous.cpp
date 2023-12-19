#include <nori/medium.h>
#include <nori/shape.h>
NORI_NAMESPACE_BEGIN
using namespace std;

class HomogeneousMedium : public Medium {
public:
    HomogeneousMedium(const PropertyList& propList) {
        m_isRGB = propList.getBoolean("isRGB", true);
        if (m_isRGB) {
            m_albedo = propList.getColor("albedo", Color3f(0.5f));
            m_sigma_t = propList.getColor("sigma_t", Color3f(0.5f));
            m_scale = propList.getFloat("scale", 1.f);
            m_sigma_t *= m_scale;
            m_sigma_s = m_albedo * m_sigma_t;
            m_sigma_a = m_sigma_t - m_sigma_s;
        }
        else {
            m_albedo = propList.getColor("albedo", Color3f(0.5f));
            m_sigmaT_f = propList.getFloat("sigma_t", 0.5f);
            m_scale = propList.getFloat("scale", 1.f);
            m_sigmaT_f *= m_scale;
            m_sigmaS_f = m_albedo_f * m_sigmaT_f;
            m_sigmaA_f = m_sigmaT_f - m_sigmaS_f;
        }

    }

    bool sample_freepath(MediumQueryRecord& mRec, Sampler* sampler) const {
        if (m_isRGB) {
            // EBalance sampling for channel-variant extinction, pick a random channel each time
            int channel = std::min((int)(sampler->next1D() * 3.f), 2);
            float density = m_sigma_t[channel];
            if (density < Epsilon) {
                mRec.p = mRec.ref + mRec.tMax * -mRec.wi;
                mRec.ret = 1;
                return false;
            }
            float t = -log(1.f - sampler->next1D()) / density;
            float pdf_failure = 0;
            float pdf_success = 0;
            float sampled_distance = t < mRec.tMax ? t : mRec.tMax;
            bool valid = t < mRec.tMax ? true : false;
            for (int i = 0; i < 3; ++i) {
                float tmp = exp(-m_sigma_t[i] * sampled_distance);
                pdf_failure += tmp;
                pdf_success += m_sigma_t[i] * tmp;
            }
            pdf_success /= 3.f;
            pdf_failure /= 3.f;
            Color3f transmittance = (-m_sigma_t * sampled_distance).exp();
            if (valid) {
                mRec.p = mRec.ref + t * -mRec.wi;
                if (mRec.p == mRec.ref) return false;
                mRec.ret = m_sigma_s * transmittance / pdf_success;
            }
            else {
                mRec.p = mRec.ref + mRec.tMax * -mRec.wi;
                mRec.ret = transmittance / pdf_failure;
            }
            return valid;
        }else {
            // scala value for density function
            bool deltaTrack = false;
            if (deltaTrack) {
                //perform delta tracking
                Ray3f ray(mRec.ref, -mRec.wi, 0, mRec.tMax);
                float pdf_failure = 1.0f;
                float pdf_success = 1.0f;
                Color3f transmittance = Color3f(1.0f);
                float mint, maxt;
                if (!m_shape->getBoundingBox().rayIntersect(ray, mint, maxt))
                    return false;
                mint = std::max(mint, ray.mint);
                maxt = std::min(maxt, ray.maxt);
                float t = mint;
                float densityP = 0;
                while (true) {
                    t += -log(1 - sampler->next1D()) * 1.0f / m_sigmaT_f;
                    if (t >= mRec.tMax)
                        break;
                    Point3f p = ray(t);
                    densityP = m_sigmaT_f;
                    if (densityP * 1.0f / m_sigmaT_f > sampler->next1D()) {
                        mRec.p = p;
                        mRec.ret = Color3f(m_sigmaS_f / m_sigmaT_f);
                        return true;
                    }
                }
                mRec.p = mRec.ref + mRec.tMax * (- mRec.wi);
                mRec.ret = Color3f(1.f);
                return false;
            }else {
                // naive update of constant albedo and extinction
                float t = -log(1.f - sampler->next1D()) / m_sigmaT_f;

                float sampled_distance = t < mRec.tMax ? t : mRec.tMax;
                bool valid = t < mRec.tMax ? true : false;
                if (valid) {
                    mRec.p = mRec.ref + t * (- mRec.wi);
                    if (mRec.p == mRec.ref) return false;
                    mRec.ret = Color3f(m_albedo);
                }
                else {
                    mRec.p = mRec.ref + mRec.tMax * (- mRec.wi);
                    mRec.ret = Color3f(1.0);
                }
                return valid;
            }
        }
    }

    Color3f eval(const MediumQueryRecord& mRec) const {
        return m_albedo;
    }

    bool isHomogeneous() const override { return true; }

    Color3f& getSigmaA() { return m_sigma_a; }
    Color3f& getSigmaS() { return m_sigma_s; }
    Color3f& getSigmaT() { return m_sigma_t; }
    Color3f& getAlbedo() { return m_albedo; }

    void addChild(NoriObject* obj) {
        switch (obj->getClassType()) {
        case EPhaseFunction:
            if (m_phase)
                throw NoriException(
                    "Already have a phase function");
            m_phase = static_cast<PhaseFunction*>(obj);
            break;

        default:
            throw NoriException("Invalid",
                classTypeName(obj->getClassType()));
        }
    }
    std::string toString() const override {
        return tfm::format(
            "HomogeneousMedium[\n"
            "  sigmaA = %s,\n"
            "  sigmaS = %s,\n"
            "  sigmaT = %s,\n"
            "  albedo = %s,\n"
            "]",
            m_sigma_a,
            m_sigma_s,
            m_sigma_t,
            m_albedo);
    }

protected:
    Color3f m_sigma_a;
    Color3f m_sigma_s;
    Color3f m_sigma_t;
    Color3f m_albedo;
    float m_sigmaT_f;
    float m_albedo_f;
    float m_sigmaS_f;
    float m_sigmaA_f;
    float m_scale;
    bool m_isRGB;
};

NORI_REGISTER_CLASS(HomogeneousMedium, "homogeneous");
NORI_NAMESPACE_END