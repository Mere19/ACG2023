#include <nori/phasefunction.h>
#include <nori/warp.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class HenyeyGreenstein : public PhaseFunction {
public:
    HenyeyGreenstein(const PropertyList& props) {
        m_g = props.getFloat("g", 0.f);
    }
    float sample(PhaseFunctionQueryRecord& pRec, const Point2f& sample) const override{

        float cosTheta = 0.f;
        if (std::abs(m_g) < Epsilon) {
            cosTheta = 1 - 2 * sample.x();
        }
        else {
            float value = (1 - m_g * m_g) / (1 - m_g + 2 * m_g * sample.x());
            cosTheta = (1 + m_g * m_g - value * value) / (2 * m_g);
        }

        float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);
        float phi = 2.0 * M_PI * sample.y();
        float sinPhi = std::sin(phi);
        float cosPhi = std::cos(phi);

        pRec.wo = Frame(-pRec.wi).toWorld(Vector3f(
            sinTheta * cosPhi,
            sinTheta * sinPhi,
            cosTheta
        ));

        return 1.0f;
    }
    float pdf(PhaseFunctionQueryRecord& pRec) const override{
        float value = 1.0f + m_g * m_g + 2.0f * m_g * pRec.wi.dot(pRec.wo);
        return INV_FOURPI * (1 - m_g * m_g) / (value * std::sqrt(value));
    }
    float eval(PhaseFunctionQueryRecord& pRec) const override{
        float value = 1.0f + m_g * m_g + 2.0f * m_g * pRec.wi.dot(pRec.wo);
        return INV_FOURPI * (1 - m_g * m_g) / (value * std::sqrt(value));
    }
    std::string toString() const override {
        return tfm::format("HenyeyGreenstein[]");
    }
protected:
    float m_g;
};

NORI_REGISTER_CLASS(HenyeyGreenstein, "hg");
NORI_NAMESPACE_END