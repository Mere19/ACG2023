#include<nori/volume.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class ConstantFloatVolume : public Volume {
public:
    ConstantFloatVolume(const PropertyList& props) {
        m_value = props.getFloat("value");
    }
    virtual float lookupFloat(const Point3f& p) const override {
        return m_value;
    };
    virtual Color3f lookupRGB(const Point3f& p) const override {
        //not used
        return Color3f(0.0f);
    };
    std::string toString() const override {
        return tfm::format("ConstantFloatVolume[]");
    }
protected:
    float m_value;
};

NORI_REGISTER_CLASS(ConstantFloatVolume, "constFVol");
NORI_NAMESPACE_END