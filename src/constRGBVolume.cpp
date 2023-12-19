#include<nori/volume.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class ConstantRGBVolume : public Volume {
public:
    ConstantRGBVolume(const PropertyList& props) {
        m_value = props.getColor("value");
    }
    virtual float lookupFloat(const Point3f& p) const override {
        //not used
        return 0.0f;
    };
    virtual Color3f lookupRGB(const Point3f& p) const override {
        return m_value;
    };
    std::string toString() const override {
        return tfm::format("ConstantFloatVolume[]");
    }
protected:
    Color3f m_value;
};

NORI_REGISTER_CLASS(ConstantRGBVolume, "constRGBVol");
NORI_NAMESPACE_END